from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from flask import current_app
from app import get_all_users_data, save_all_users_data

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
