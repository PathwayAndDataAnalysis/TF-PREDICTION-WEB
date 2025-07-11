import shutil
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import os

import psutil
from werkzeug.utils import secure_filename
import pandas as pd
from flask import current_app, jsonify
import scanpy as sc
import numpy as np
from . import get_all_users_data, save_all_users_data, get_file_path
from filelock import FileLock

executor = ThreadPoolExecutor(max_workers=4)  # Adjust as needed

MAX_FILE_SIZE = 500 * 1024 * 1024 * 1024  # 500 GB
MIN_DISK_SPACE = 10 * 1024 * 1024 * 1024  # 10 GB

# Define constants for column names
CLUSTER_COL = "Cluster"
UMAP1_COL = "X_umap1"
UMAP2_COL = "X_umap2"
PCA1_COL = "X_pca1"
PCA2_COL = "X_pca2"


def run_in_background(fn, *args, **kwargs):
    app = current_app._get_current_object()

    def wrapped(*args, **kwargs):
        with app.app_context():
            try:
                return fn(*args, **kwargs)
            except Exception as e:
                current_app.logger.error(f"Background job {fn.__name__} failed: {e}")
                raise e

    current_app.logger.info(
        f"[UTILS] Submitting background job: {fn.__name__} with args={args} kwargs={kwargs}"
    )
    future = executor.submit(wrapped, *args, **kwargs)
    current_app.logger.info(f"[UTILS] Background job {fn.__name__} submitted. Future: {future}")
    return future


def _get_summary_stats(data_series):
    """
    Calculates summary statistics for a pandas Series, handling NaNs.
    """
    stats = data_series.describe()
    if stats.get("count", 0.0) == 0.0:
        return {
            "mean": 0.0,
            "sd": 0.0,
            "min": 0.0,
            "max": 0.0
        }
    return {
        "mean": round(float(stats.get("mean", 0)), 2),
        "sd": round(float(stats.get("std", 0)), 2),
        "min": round(float(stats.get("min", 0)), 2),
        "max": round(float(stats.get("max", 0)), 2),
    }


def calculate_and_save_qc_metrics(user_id, filename, file_path):
    """
    Reads an uploaded file, calculates QC metrics using Scanpy,
    and saves the histogram data to users.json.
    """
    try:
        if not file_path:
            current_app.logger.error(f"Error: Could not find path for {filename} for user {user_id}")
            return

        BINS = 100
        qc_results = {}

        # Read data into AnnData object
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad(file_path)

            # Get only the non-zero expression values, which is a much smaller array
            non_zero_expr = adata.X.data

            if non_zero_expr.size > 0:
                # Calculate histogram only on the non-zero values
                expr_counts, expr_bins = np.histogram(non_zero_expr, bins=BINS)

                # Calculate the number of zero values separately
                total_elements = adata.n_obs * adata.n_vars
                num_zeros = total_elements - len(non_zero_expr)

                # Add the count of zeros to the first bin of the histogram
                expr_counts[0] += num_zeros

                expr_stats = {
                    "mean": round(float(non_zero_expr.sum() / total_elements), 3),  # Mean across all elements
                    "sd": round(float(np.std(non_zero_expr)), 2),
                    "min": round(float(np.min(non_zero_expr)), 2),
                    "max": round(float(np.max(non_zero_expr)), 2),
                }
            else:
                expr_counts, expr_bins = [], []
                expr_stats = {"mean": 0.0, "sd": 0.0, "min": 0.0, "max": 0.0}

            qc_results["gene_expression"] = {
                "counts": expr_counts.tolist(),
                "bins": expr_bins.tolist(),
                **expr_stats
            }

        else:  # For .tsv, .csv
            # Infer delimiter based on file extension
            delimiter = infer_delimiter(file_path)
            adata = sc.AnnData(pd.read_csv(file_path, index_col=0, delimiter=delimiter))
            gene_exp_df = pd.read_csv(file_path, index_col=0, delimiter=delimiter)

            expr_flat = gene_exp_df.values.flatten()
            expr_flat = expr_flat[np.isfinite(expr_flat)]  # remove NaN/Inf
            if expr_flat.size > 0:
                expr_counts, expr_bins = np.histogram(expr_flat, bins=BINS)
                expr_stats = {
                    "mean": round(float(np.mean(expr_flat)), 2),
                    "sd": round(float(np.std(expr_flat)), 2),
                    "min": round(float(np.min(expr_flat)), 2),
                    "max": round(float(np.max(expr_flat)), 2),
                }
            else:
                expr_counts, expr_bins = [], []
                expr_stats = {"mean": 0.0, "sd": 0.0, "min": 0.0, "max": 0.0}
            qc_results["gene_expression"] = {
                "counts": expr_counts.tolist(),
                "bins": expr_bins.tolist(),
                **expr_stats
            }

        # Identify mitochondrial genes (works for both human 'MT-' and mouse 'mt-')
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Generate histogram data (lightweight for JSON storage)
        metrics = {
            "n_genes_by_counts": adata.obs['n_genes_by_counts'],  # Genes per cell
            "n_cells_by_counts": adata.var['n_cells_by_counts'],  # Cells per gene
            "pct_counts_mt": adata.obs['pct_counts_mt'],  # MT percentage
        }

        for name, data in metrics.items():
            valid_data = data.dropna()
            if not valid_data.empty:
                counts, bins = np.histogram(valid_data, bins=BINS)
                qc_results[name] = {
                    "counts": counts.tolist(),
                    "bins": bins.tolist(),
                    **_get_summary_stats(valid_data)
                }
            else:
                qc_results[name] = {"counts": [], "bins": []}

        qc_status = "completed"

    except Exception as e:
        current_app.logger.error(f"Error calculating QC for {filename}: {e}")
        qc_results = {}
        qc_status = "failed"

    # --- Safely update users.json ---
    users_json_path = "instance/users.json"
    lock_path = users_json_path + ".lock"

    with FileLock(lock_path):
        all_users_data = get_all_users_data()
        user_node = all_users_data.get(user_id)
        if user_node and 'files' in user_node:
            for i, file_info in enumerate(user_node['files']):
                if file_info['filename'] == filename:
                    # Add qc metrics and status to the file's record
                    all_users_data[user_id]['files'][i]['qc_metrics'] = qc_results
                    all_users_data[user_id]['files'][i]['qc_status'] = qc_status
                    break
        save_all_users_data(all_users_data)
        current_app.logger.info(
            f"[UTILS] QC metrics for {filename} saved successfully for user {user_id}."
        )


def update_analysis_status(
    user_id,
    analysis_id,
    status,
    umap_csv_path=None,
    metadata_cols=None,
    tfs=None,
    pvalues_path=None,
    activation_path=None,
    bh_reject_path=None,
    fdr_level=None,
    p_value_threshold=None,
    p_val_threshold_path=None,
    z_scores_path=None,
    error=None,
):
    try:
        current_app.logger.info(
            f"[UTILS] Updating analysis status: user_id={user_id}, analysis_id={analysis_id}, status={status}, error={error}"
        )

        all_users_data = get_all_users_data()
        user_node = all_users_data.get(user_id, {})

        if not user_node:
            current_app.logger.error(f"[UTILS] User {user_id} not found in users data")
            return

        found = False
        for analysis in user_node.get("analyses", []):
            if analysis["id"] == analysis_id:
                if status:
                    analysis["status"] = status
                if metadata_cols:
                    analysis["metadata_cols"] = metadata_cols
                if tfs:
                    analysis["tfs"] = tfs
                if pvalues_path:
                    analysis["pvalues_path"] = pvalues_path
                if activation_path:
                    analysis["activation_path"] = activation_path
                if bh_reject_path:
                    analysis["bh_reject_path"] = bh_reject_path
                if p_value_threshold or fdr_level:
                    if p_value_threshold:
                        analysis.get("inputs")["p_value_threshold"] = p_value_threshold
                        if "fdr_level" in analysis.get("inputs", {}):
                            del analysis["inputs"]["fdr_level"]
                    else:
                        analysis["inputs"]["fdr_level"] = fdr_level
                        if "p_value_threshold" in analysis.get("inputs", {}):
                            del analysis["inputs"]["p_value_threshold"]
                if p_val_threshold_path:
                    analysis["p_val_threshold_path"] = p_val_threshold_path
                if z_scores_path:
                    analysis["z_scores_path"] = z_scores_path
                if umap_csv_path:
                    layout = analysis.get("inputs", {}).get("layout", {})
                    layout["layout_filepath"] = umap_csv_path
                if error:
                    analysis["error"] = error
                    analysis["error_timestamp"] = datetime.now().isoformat()
                found = True
                break

        if found:
            save_all_users_data(all_users_data)
            current_app.logger.info(
                f"[UTILS] Analysis status updated and saved for analysis_id={analysis_id}."
            )
        else:
            current_app.logger.warning(
                f"[UTILS] Analysis with id={analysis_id} not found for user_id={user_id}. No update performed."
            )

    except Exception as e:
        current_app.logger.error(f"[UTILS] Error updating analysis status: {e}")
        # Don't raise the exception to prevent cascading failures


def infer_delimiter(filepath):
    ext = filepath.split('.')[-1]
    if ext == 'tsv':
        return '\t'
    elif ext == 'csv':
        return ','
    else:
        current_app.logger.warning(
            f"[UTILS] Unsupported file extension '{ext}' for delimiter inference. Defaulting to comma."
        )
        return ','


def check_system_resources():
    """Check if system has enough resources for processing."""
    try:
        # Check disk space
        disk_usage = shutil.disk_usage(current_app.config["UPLOAD_FOLDER"])
        if disk_usage.free < MIN_DISK_SPACE:
            raise ValueError(f"Insufficient disk space. Need at least {MIN_DISK_SPACE / (1024 ** 3):.1f} GB free.")

        # Check memory
        memory = psutil.virtual_memory()
        if memory.available < 2 * 1024 * 1024 * 1024:  # 2 GB
            raise ValueError("Insufficient memory available for processing.")

        return True
    except Exception as e:
        current_app.logger.error(f"System resource check failed: {e}")
        return False


def get_user_specific_data_path(username):
    """Returns the path to a user's dedicated data directory."""
    # current_app.config['UPLOAD_FOLDER'] is the root for all user uploads
    base_upload_folder = current_app.config["UPLOAD_FOLDER"]
    # Sanitize username for directory creation, though secure_filename is usually for files
    safe_username_dir = secure_filename(username)  # Ensures username is safe for path
    path = os.path.join(base_upload_folder, safe_username_dir)
    try:
        os.makedirs(path, exist_ok=True)
    except OSError as e:
        current_app.logger.error(f"Could not create user data path {path}: {e}")
        return None  # Indicate failure
    return path


def get_user_analysis_path(username, analysis_id):
    """Returns the path to a specific analysis directory for a user."""
    user_data_path = get_user_specific_data_path(username)
    if not user_data_path:
        return None
    analysis_base_dir_name = "analyses"  # Name of the subfolder for all analyses
    analysis_path = os.path.join(
        user_data_path, analysis_base_dir_name, secure_filename(analysis_id)
    )
    try:
        os.makedirs(analysis_path, exist_ok=True)  # Create if it doesn't exist
    except OSError as e:
        current_app.logger.error(f"Could not create analysis path {analysis_path}: {e}")
        return None
    return analysis_path


def allowed_file(filename):
    return (
            "." in filename
            and filename.rsplit(".", 1)[1].lower()
            in current_app.config["ALLOWED_EXTENSIONS"]
    )


def estimate_fdr_for_gene(p_values_df, gene_name, p_value_cutoff):
    if gene_name not in p_values_df.columns:
        return np.nan

    p_gene = p_values_df[[gene_name]].dropna()
    m = p_gene.shape[0]  # Total tests for this gene
    if m == 0:
        return np.nan

    discoveries = (p_gene < p_value_cutoff).sum().iloc[0]

    if discoveries > 0:
        fdr = (p_value_cutoff * m) / discoveries
        return fdr
    else:
        return np.nan  # Or 0.0 if you prefer


def generate_scatter_plot_response(analysis_to_view, plot_type=None):
    try:
        layout_filepath = (analysis_to_view.get("inputs", {}).get("layout", {}).get("layout_filepath", ""))

        if not os.path.exists(layout_filepath):
            current_app.logger.error(f"Layout file '{layout_filepath}' not found!")
            return jsonify({"error": "Layout file not found."}), 404

        # Determine which coordinates to load
        if plot_type == 'umap_plot':
            column_names = [UMAP1_COL, UMAP2_COL]
            plot_title = "UMAP Plot"
        elif plot_type == 'pca_plot':
            plot_title = "PCA Plot"
            column_names = [PCA1_COL, PCA2_COL]
        else:
            return jsonify({"error": "Invalid plot type specified."}), 400

        plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))

        traces = []
        for cluster_label, group_df in plot_df.groupby(CLUSTER_COL):
            traces.append(
                {
                    "cluster": str(cluster_label),
                    "x": group_df[column_names[0]].tolist(),
                    "y": group_df[column_names[1]].tolist(),
                    "mode": "markers",
                    "type": "scattergl",
                    "name": str(cluster_label),
                }
            )

        layout = {
            "title": plot_title,
            "xaxis": {"title": column_names[0]},
            "yaxis": {"title": column_names[1]},
        }

        graph_data = {
            "data": traces,
            "layout": layout,
            "metadata_cols": analysis_to_view.get("metadata_cols", []),
            "tfs": analysis_to_view.get("tfs", [])
        }
        return jsonify(graph_data), 200

    except Exception as e:
        current_app.logger.error(
            f"Unexpected error for analysis_id '{analysis_to_view['id']}': {e}",
            exc_info=True,
        )
        return jsonify({"error": "An unexpected error occurred while preparing the plot."}), 500


def get_layout_and_metadata_dfs(analysis, user_id):
    layout_source = analysis.get("inputs", {}).get("layout", {}).get("source")

    layout_filepath = analysis.get("inputs").get("layout").get("layout_filepath")
    plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))

    if layout_source == "FILE":
        metadata_filepath = analysis.get("inputs").get("gene_expression").get("metadata_filepath", {})
        metadata_df = pd.read_csv(get_file_path(metadata_filepath, user_id), index_col=0,
                                  sep=infer_delimiter(metadata_filepath))
    else:
        h5ad_file = analysis.get("inputs").get("gene_expression").get("h5ad_filepath")
        adata = sc.read_h5ad(get_file_path(h5ad_file, user_id))
        metadata_df = pd.DataFrame(adata.obs)

    return plot_df, metadata_df


def get_layout_and_gene_exp_levels_df(analysis, gene_name):
    try:
        layout_filepath = (analysis.get("inputs").get("layout").get("layout_filepath"))
        z_score_filepath = analysis.get("z_scores_path", "")

        if not os.path.exists(layout_filepath):
            current_app.logger.error(f"Layout file '{layout_filepath}' not found for analysis_id '{analysis['id']}'.")
            return jsonify({"error": "Layout file not found. Delete this analysis and create new analysis"}), 404

        if not os.path.exists(z_score_filepath):
            current_app.logger.error(f"Z-scores file '{z_score_filepath}' not found.")
            return jsonify({"error": "Z-scores file not found. Delete this analysis and create new analysis"}), 404

        plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))
        gene_exp_levels_df = pd.read_parquet(z_score_filepath, use_threads=True, columns=[gene_name])

        plot_df = plot_df.merge(gene_exp_levels_df, left_index=True, right_index=True)
        plot_df = plot_df.dropna(subset=[gene_name])

        return plot_df
    except Exception as e:
        current_app.logger.error(f"Error reading layout/z-score file for analysis_id '{analysis['id']}': {e}")
        return jsonify({"error": "Failed to read layout/z-score file. Invalid format."}), 500


def generate_colored_traces(
        plot_df, plot_type="umap_plot", cluster_col="Cluster", tf_activity=None
):
    if plot_type.lower() == "umap_plot":
        x_col, y_col = UMAP1_COL, UMAP2_COL
        title = (
            f"UMAP Plot Colored by {cluster_col if not tf_activity else tf_activity}"
        )
    elif plot_type.lower() == "pca_plot":
        x_col, y_col = PCA1_COL, PCA2_COL
        title = f"PCA Plot Colored by {cluster_col if not tf_activity else tf_activity}"
    else:
        current_app.logger.warning(f"Unknown plot type '{plot_type}'")
        return jsonify({"error": f"Unknown plot type: {plot_type}"}), 400

    traces = []
    if tf_activity:
        plot_df[tf_activity] = plot_df[tf_activity].replace(
            {1: "Active", -1: "Inactive", 0: "Insignificant", np.nan: "Not_Enough_Data"}
        )

        unique_clusters = {
            "Active": "red",
            "Inactive": "blue",
            "Insignificant": "gray",
            "Not_Enough_Data": "yellow",
        }
        for cluster in unique_clusters.keys():
            traces.append(
                {
                    "cluster": cluster,
                    "x": plot_df[plot_df[tf_activity] == cluster][x_col].tolist(),
                    "y": plot_df[plot_df[tf_activity] == cluster][y_col].tolist(),
                    "mode": "markers",
                    "type": "scattergl",
                    "name": cluster,
                    "marker": {
                        "color": unique_clusters[cluster],
                    },
                }
            )
    else:
        for cluster in plot_df[cluster_col].unique().tolist():
            traces.append(
                {
                    "cluster": cluster,
                    "x": plot_df[plot_df["Cluster"] == cluster][x_col].tolist(),
                    "y": plot_df[plot_df["Cluster"] == cluster][y_col].tolist(),
                    "mode": "markers",
                    "type": "scattergl",
                    "name": cluster,
                }
            )
    return traces, title, x_col, y_col
