from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

import pandas as pd
from flask import current_app
import scanpy as sc
import numpy as np
from . import get_file_path, get_all_users_data, save_all_users_data
from filelock import FileLock

executor = ThreadPoolExecutor(max_workers=4)  # Adjust as needed


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
    artifact_path=None,
    metadata_cols=None,
    tfs=None,
    bh_reject_path=None,
    fdr_level=None,
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
                if bh_reject_path:
                    analysis["bh_reject_path"] = bh_reject_path
                if fdr_level:
                    analysis.get("inputs", {})["fdr_level"] = fdr_level
                if artifact_path:
                    analysis["artifact_path"] = artifact_path
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
