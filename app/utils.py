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


def calculate_and_save_qc_metrics(user_id, filename, file_path):
    """
    Reads an uploaded file, calculates QC metrics using Scanpy,
    and saves the histogram data to users.json.
    """
    try:
        if not file_path:
            current_app.logger.error(f"Error: Could not find path for {filename} for user {user_id}")
            return

        # Read data into AnnData object
        if filename.endswith('.h5ad'):
            adata = sc.read_h5ad(file_path)
        else:  # For .tsv, .csv
            adata = sc.AnnData(pd.read_csv(file_path, index_col=0))

        # Identify mitochondrial genes (works for both human 'MT-' and mouse 'mt-')
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        # Generate histogram data (lightweight for JSON storage)
        qc_results = {}
        metrics = {
            "n_genes_by_counts": adata.obs['n_genes_by_counts'],  # Genes per cell
            "n_cells_by_counts": adata.var['n_cells_by_counts'],  # Cells per gene
            "pct_counts_mt": adata.obs['pct_counts_mt'],  # MT percentage
        }

        for name, data in metrics.items():
            counts, bins = np.histogram(data.dropna(), bins=300)
            qc_results[name] = {
                "counts": counts.tolist(),
                "bins": bins.tolist()
            }

        mean_n_genes_by_counts = adata.obs['n_genes_by_counts'].mean()
        sd_n_genes_by_counts = adata.obs['n_genes_by_counts'].std()
        min_n_genes_by_counts = adata.obs['n_genes_by_counts'].min()
        max_n_genes_by_counts = adata.obs['n_genes_by_counts'].max()

        mean_n_cells_by_counts = adata.var['n_cells_by_counts'].mean()
        sd_n_cells_by_counts = adata.var['n_cells_by_counts'].std()
        min_n_cells_by_counts = adata.var['n_cells_by_counts'].min()
        max_n_cells_by_counts = adata.var['n_cells_by_counts'].max()

        mean_pct_counts_mt = adata.obs['pct_counts_mt'].mean()
        sd_pct_counts_mt = adata.obs['pct_counts_mt'].std()
        min_pct_counts_mt = adata.obs['pct_counts_mt'].min()
        max_pct_counts_mt = adata.obs['pct_counts_mt'].max()

        data_summary = {
            "mean_n_genes_by_counts": round(float(mean_n_genes_by_counts), 2),
            "sd_n_genes_by_counts": round(float(sd_n_genes_by_counts), 2),
            "min_n_genes_by_counts": round(float(min_n_genes_by_counts), 2),
            "max_n_genes_by_counts": round(float(max_n_genes_by_counts), 2),
            "mean_n_cells_by_counts": round(float(mean_n_cells_by_counts), 2),
            "sd_n_cells_by_counts": round(float(sd_n_cells_by_counts), 2),
            "min_n_cells_by_counts": round(float(min_n_cells_by_counts), 2),
            "max_n_cells_by_counts": round(float(max_n_cells_by_counts), 2),
            "mean_pct_counts_mt": round(float(mean_pct_counts_mt), 2),
            "sd_pct_counts_mt": round(float(sd_pct_counts_mt), 2),
            "min_pct_counts_mt": round(float(min_pct_counts_mt), 2),
            "max_pct_counts_mt": round(float(max_pct_counts_mt), 2)
        }

        qc_status = "completed"
        qc_results['data_summary'] = data_summary

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
                if fdr_level:
                    analysis.get("inputs", {})["fdr_level"] = fdr_level
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


# def infer_delimiter(filepath):
#     with open(filepath, 'r') as f:
#         first_line = f.readline()
#         return '\t' if '\t' in first_line else ','

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
