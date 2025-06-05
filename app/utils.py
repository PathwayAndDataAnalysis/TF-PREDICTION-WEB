from concurrent.futures import ThreadPoolExecutor

from flask import current_app

from app import get_all_users_data, save_all_users_data

executor = ThreadPoolExecutor(max_workers=4)  # Adjust as needed


def run_in_background(fn, *args, **kwargs):
    app = current_app._get_current_object()  # Get the actual Flask app object

    def wrapped(*args, **kwargs):
        with app.app_context():
            return fn(*args, **kwargs)

    current_app.logger.info(
        f"[UTILS] Submitting background job: {fn.__name__} with args={args} kwargs={kwargs}"
    )
    future = executor.submit(wrapped, *args, **kwargs)
    current_app.logger.info(
        f"[UTILS] Background job {fn.__name__} submitted. Future: {future}"
    )
    return future


def update_analysis_status(
        user_id, analysis_id, status, umap_csv_path=None, metadata_cols=None, tfs=None, pvalues_path=None, bh_reject_path=None, error=None
):
    if metadata_cols is None:
        metadata_cols = []
    current_app.logger.info(
        f"[UTILS] Updating analysis status: user_id={user_id}, analysis_id={analysis_id}, status={status}, umap_csv_path={umap_csv_path}, metadata_cols={metadata_cols}, error={error}"
    )
    all_users_data = get_all_users_data()
    user_node = all_users_data.get(user_id, {})
    found = False
    for analysis in user_node.get("analyses", []):
        if analysis["id"] == analysis_id:
            analysis["status"] = status
            analysis["metadata_cols"] = metadata_cols
            if tfs:
                analysis["tfs"] = tfs
            if pvalues_path:
                analysis["pvalues_path"] = pvalues_path
            if bh_reject_path:
                analysis["bh_reject_path"] = bh_reject_path
            if umap_csv_path:
                analysis["umap_csv"] = umap_csv_path
            if error:
                analysis["error"] = error
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
