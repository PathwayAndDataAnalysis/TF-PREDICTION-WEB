import os
import uuid
from datetime import datetime, timezone  # For timestamps

import numpy as np
import pandas as pd
import scanpy as sc
from flask import (
    Blueprint,
    render_template,
    request,
    redirect,
    url_for,
    flash,
    current_app,
    jsonify, send_from_directory,
)
from flask_login import login_user, logout_user, login_required, current_user
from werkzeug.security import generate_password_hash, check_password_hash
from werkzeug.utils import secure_filename
import psutil
import shutil
from werkzeug.exceptions import RequestEntityTooLarge
from . import (
    User,
    get_all_users_data,
    save_all_users_data,
    find_analysis_by_id,
    get_file_path,
)
from .benjamini_hotchberg import bh_fdr_correction
from .tf_analysis import run_tf_analysis
from .umap_pipeline import run_umap_pipeline
from .utils import run_in_background, update_analysis_status, infer_delimiter, calculate_and_save_qc_metrics

# Create Blueprints
auth_bp = Blueprint("auth", __name__, url_prefix="/auth")
main_bp = Blueprint("main_routes", __name__)  # No prefix for main app routes like index

MAX_FILE_SIZE = 500 * 1024 * 1024 * 1024  # 500 GB
MIN_DISK_SPACE = 10 * 1024 * 1024 * 1024  # 10 GB

# Define constants for column names
CLUSTER_COL = "Cluster"
UMAP1_COL = "X_umap1"
UMAP2_COL = "X_umap2"
PCA1_COL = "X_pca1"
PCA2_COL = "X_pca2"


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


# --- Helper functions (get_user_specific_data_path, allowed_file as before) ---
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


# --- Jinja Filter for Datetime Formatting ---
@main_bp.app_template_filter()
def format_datetime(value, format="%Y-%m-%d %H:%M"):
    """Formats a datetime string or object."""
    if isinstance(value, str):
        try:
            # Attempt to parse if it's a common ISO format string
            dt_obj = datetime.fromisoformat(value.replace("Z", "+00:00"))
            return dt_obj.strftime(format)
        except ValueError:
            return value  # Return original if parsing fails
    if isinstance(value, datetime):
        return value.strftime(format)
    return value

# This filter will format bytes into KB, MB, GB, etc.
@main_bp.app_template_filter()
def format_filesize(value):
    if value is None:
        return "N/A"
    try:
        size_bytes = int(value)
        if size_bytes < 1024:
            return f"{size_bytes} B"
        elif size_bytes < 1024**2:
            return f"{size_bytes/1024:.1f} KB"
        elif size_bytes < 1024**3:
            return f"{size_bytes/1024**2:.1f} MB"
        else:
            return f"{size_bytes/1024**3:.2f} GB"
    except (ValueError, TypeError):
        return value


# --- Main Routes (index, data operations) ---
@main_bp.route("/")
def index():
    if not current_user.is_authenticated:
        return redirect(url_for("auth.login"))

    all_users_data = get_all_users_data()
    current_user_node = all_users_data.get(current_user.id, {})
    user_files = current_user_node.get("files", [])
    user_analyses = current_user_node.get("analyses", [])

    # Create a mapping from secured filename to original filename for easy lookup
    filename_map = {f['filename']: f.get('original_filename', f['filename']) for f in user_files}

    for analysis in user_analyses:
        inputs = analysis.get('inputs', {})
        display_inputs = {}

        # Enrich Gene Expression input
        ge_input = inputs.get('gene_expression', {})
        if ge_input.get('h5ad_filepath'):
            # The h5ad_filepath is already the secured name
            secured_name = ge_input['h5ad_filepath'].split("/")[-1]
            display_inputs['gene_expression'] = filename_map.get(secured_name, secured_name)
        elif ge_input.get('gene_exp_filepath'):
            secured_name = ge_input['gene_exp_filepath'].split("/")[-1]
            display_inputs['gene_expression'] = filename_map.get(secured_name, secured_name)
            if ge_input.get('metadata_filepath'):
                meta_name = ge_input['metadata_filepath'].split("/")[-1]
                display_inputs['metadata'] = filename_map.get(meta_name, meta_name)

        # Enrich Layout input
        layout_input = inputs.get('layout', {})
        if layout_input.get('source') == 'FILE' and layout_input.get('layout_filepath'):
            layout_name = layout_input['layout_filepath'].split("/")[-1]
            display_inputs['layout'] = filename_map.get(layout_name, layout_name)
        elif layout_input.get('source') == 'UMAP_GENERATED':
            display_inputs['layout'] = "Generated from data"

        # Add the enriched data to the analysis dictionary
        analysis['display_inputs'] = display_inputs

    user_analyses_sorted = sorted(user_analyses, key=lambda x: x.get("created_at", ""), reverse=True)

    return render_template(
        "index.html",
        user_files=user_files,
        current_user_analyses=user_analyses_sorted,
    )


@main_bp.route("/download_data/<filename>")
@login_required
def download_data(filename):
    # Security check: Ensure the requested file belongs to the current user
    all_users_data = get_all_users_data()
    user_files = all_users_data.get(current_user.id, {}).get("files", [])

    if not any(f['filename'] == filename for f in user_files):
        flash("File not found or access denied.", "error")
        return redirect(url_for('main_routes.index'))

    # Get the user's upload directory path
    user_dir = get_user_specific_data_path(current_user.id)
    if not user_dir:
        flash("Could not locate user directory.", "error")
        return redirect(url_for('main_routes.index'))

    try:
        # Use Flask's secure send_from_directory
        return send_from_directory(user_dir, filename, as_attachment=True)
    except FileNotFoundError:
        flash("File not found on disk.", "error")
        return redirect(url_for('main_routes.index'))


@main_bp.route("/upload_data", methods=["POST"])
@login_required
def upload_data():
    try:
        # Check system resources first
        if not check_system_resources():
            flash("System resources insufficient for file upload. Please try again later.", "error")
            return redirect(url_for("main_routes.index"))

        if "data_file" not in request.files:
            flash("No file part in the request.", "error")
            return redirect(url_for("main_routes.index"))

        file = request.files["data_file"]
        description = request.form.get("description", "").strip()
        file_type = request.form.get("file_type", "Other")  # Default to 'Other'

        if file.filename == "":
            flash("No file selected for upload.", "error")
            return redirect(url_for("main_routes.index"))

        if file and allowed_file(file.filename):
            original_filename = file.filename
            filename = secure_filename(original_filename)

            # --- Decide if QC should be run ---
            should_run_qc = False
            filename_lower = filename.lower()
            if filename_lower.endswith('.h5ad'):
                should_run_qc = True
            elif (filename_lower.endswith('.csv') or filename_lower.endswith('.tsv')) and file_type == "Gene Expression":
                should_run_qc = True

            user_data_storage_path = get_user_specific_data_path(current_user.id)
            if not user_data_storage_path:
                flash("Could not access user storage directory.", "error")
                return redirect(url_for("main_routes.index"))

            destination_file_path = os.path.join(user_data_storage_path, filename)

            # if a file already exists
            all_users_data = get_all_users_data()
            if current_user.id not in all_users_data:
                all_users_data[current_user.id] = {"password": "HASH_PLACEHOLDER", "files": []}

            user_files_list = all_users_data[current_user.id].get("files", [])

            if any(f["filename"] == filename for f in user_files_list):
                flash(
                    f'File "{filename}" already exists. Please rename or delete the existing file.',
                    "warning",
                )
                return redirect(url_for("main_routes.index"))

            try:
                file.seek(0, os.SEEK_END)
                file_size = file.tell()
                file.seek(0)  # Rewind to the beginning so save() works correctly

                with open(destination_file_path, 'wb') as f:
                    chunk_size = 8192
                    while True:
                        chunk = file.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)

                # Find file_type based on extension
                file_extension = os.path.splitext(filename)[1].lower()
                if file_extension == ".h5ad":
                    file_type = "h5ad File"

                # Store file info
                user_files_list.append({
                    "filename": filename,  # The secured filename
                    "original_filename": original_filename,  # Keep original for display if needed
                    "description": description,
                    "path": destination_file_path,  # Store the full path for easier deletion
                    "file_type": file_type,
                    "file_size": file_size,
                    "upload_date": datetime.now().isoformat(),
                    **({"qc_status": "processing"} if should_run_qc else {})
                })

                all_users_data[current_user.id]["files"] = user_files_list
                save_all_users_data(all_users_data)

                if should_run_qc:
                    # --- TRIGGER BACKGROUND QC CALCULATION ---
                    run_in_background(calculate_and_save_qc_metrics, current_user.id, filename, destination_file_path)
                    flash(
                        f'File "{original_filename}" uploaded successfully! QC analysis is running in the background.',
                        "success"
                    )
                else:
                    flash(
                        f'File "{original_filename}" uploaded successfully!',
                        "success",
                    )
            except Exception as e:
                current_app.logger.error(f"Error saving file {filename} for user {current_user.id}: {e}")
                # Cleanup partial file
                if os.path.exists(destination_file_path):
                    try:
                        os.remove(destination_file_path)
                    except OSError:
                        pass
                flash(f"An error occurred during file upload: {str(e)}", "error")
        else:
            flash(f"File type not allowed. Allowed: {', '.join(current_app.config['ALLOWED_EXTENSIONS'])}", "error")

    except RequestEntityTooLarge:
        flash("File too large for upload.", "error")

    except Exception as e:
        current_app.logger.error(f"Unexpected error in upload_data: {e}")
        flash("An unexpected error occurred. Please try again.", "error")

    return redirect(url_for("main_routes.index"))


@main_bp.route("/delete_data/<filename>", methods=["POST"])
@login_required
def delete_data(filename):
    # The filename from the URL is the one stored (already secured)
    secured_filename_to_delete = secure_filename(filename)  # Re-secure to be safe

    all_users_data = get_all_users_data()
    current_user_data = all_users_data.get(current_user.id)

    if not current_user_data:
        flash("User data not found.", "error")
        return redirect(url_for("main_routes.index"))

    user_files_list = current_user_data.get("files", [])
    file_info_to_delete = None

    for f_info in user_files_list:
        if f_info.get("filename") == secured_filename_to_delete:
            file_info_to_delete = f_info
            break

    if file_info_to_delete:
        # Use the stored 'path' for deletion
        file_path_on_disk = file_info_to_delete.get("path")

        if not file_path_on_disk:
            # Fallback: construct path if 'path' wasn't stored (older data format?)
            user_data_storage_path = get_user_specific_data_path(current_user.id)
            if user_data_storage_path:
                file_path_on_disk = os.path.join(
                    user_data_storage_path, secured_filename_to_delete
                )
            else:
                flash(
                    f'Could not determine path for "{secured_filename_to_delete}". Deletion failed.',
                    "error",
                )
                return redirect(url_for("main_routes.index"))

        try:
            if os.path.exists(file_path_on_disk):
                os.remove(file_path_on_disk)
                current_app.logger.info(
                    f"Deleted file {file_path_on_disk} for user {current_user.id}"
                )
            else:
                current_app.logger.warning(
                    f"File {file_path_on_disk} not found on disk for deletion for user {current_user.id}."
                )
                # Still remove from records if not on disk, or flash a specific error

            # Remove from users.json records
            user_files_list.remove(file_info_to_delete)
            all_users_data[current_user.id]["files"] = user_files_list
            save_all_users_data(all_users_data)
            flash(
                f'File "{file_info_to_delete.get("original_filename", secured_filename_to_delete)}" deleted successfully.',
                "success",
            )
        except OSError as e:
            current_app.logger.error(
                f'Error deleting file "{secured_filename_to_delete}" from disk: {e}'
            )
            flash(f"Error deleting file from disk: {e}", "error")
    else:
        flash(f'File record for "{secured_filename_to_delete}" not found.', "error")

    return redirect(url_for("main_routes.index"))


# --- Analysis Routes ---
@main_bp.route("/analysis/create", methods=["GET"])
@login_required
def create_analysis_page():
    all_users_data = get_all_users_data()
    user_files = all_users_data.get(current_user.id, {}).get("files", [])
    return render_template("create_analysis.html", user_files=user_files)


@main_bp.route("/analysis/create", methods=["POST"])
@login_required
def create_analysis():
    try:
        # Check system resources
        if not check_system_resources():
            flash("System resources insufficient for analysis. Please try again later.", "error")
            return redirect(url_for("main_routes.create_analysis_page"))

        current_app.logger.info(f"Create analysis form data: {request.form}")

        # Validate analysis name
        analysis_name = request.form.get("analysis_name").strip()
        if not analysis_name:
            flash("Analysis name is required.", "error")
            return redirect(url_for("main_routes.create_analysis_page"))

        # Validate analysis name length and characters
        if len(analysis_name) > 100:
            flash("Analysis name too long. Maximum 100 characters.", "error")
            return redirect(url_for("main_routes.create_analysis_page"))

        # Check for existing analysis with same name
        all_users_data = get_all_users_data()
        user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
        if any(a.get("name") == analysis_name for a in user_analyses):
            flash(f'Analysis "{analysis_name}" already exists. Please choose a different name.', "error")
            return redirect(url_for("main_routes.create_analysis_page"))

        # Gene Expression Data
        have_h5ad = request.form.get("have_h5ad") == "on"
        selected_h5ad_file = request.form.get("selected_h5ad_file")
        gene_exp_file = request.form.get("gene_exp_file")
        metadata_file = request.form.get("metadata_file")
        species = request.form.get("species") if request.form.get("species") and request.form.get("species") != "select-species" else "auto"

        # 2D Layout Data
        have_2d_layout = request.form.get("have_2d_layout") == "on"
        layout_file_2d = request.form.get("layout_file_2d")

        # Data Filtering Parameters
        data_filtering = {}
        if request.form.get("filter_cells") == "on":
            data_filtering["filter_cells"] = True
            data_filtering["min_genes"] = request.form.get("min_genes", type=int)
        if request.form.get("filter_genes") == "on":
            data_filtering["filter_genes"] = True
            data_filtering["min_cells"] = request.form.get("min_cells", type=int)
        if request.form.get("qc_filter") == "on":
            data_filtering["qc_filter"] = True
            data_filtering["max_mt_pct"] = request.form.get("max_mt_pct", type=int)
        if request.form.get("data_normalize") == "on":
            data_filtering["data_normalize"] = True
            data_filtering["data_normalize_value"] = request.form.get("data_normalize_value", type=int)
        if request.form.get("log_transform") == "on":
            data_filtering["log_transform"] = True

        # UMAP Parameters (only relevant if not have_2d_layout)
        umap_parameters = None
        if not have_2d_layout:
            umap_parameters = {
                "pca_components": request.form.get("pca_components", default=20, type=int),
                "n_neighbors": request.form.get("n_neighbors", default=15, type=int),
                "min_dist": request.form.get("min_dist", default=0.3, type=float),
                "metric": request.form.get("metric", default="euclidean"),
                "random_state": request.form.get("random_state", type=int, default=0)
            }

        # Gene Expression File Validation
        if have_h5ad:
            if not selected_h5ad_file:
                flash("Please select a .h5ad file.", "error")
                return redirect(url_for("main_routes.create_analysis_page"))
        else:
            if not gene_exp_file:
                flash("Please select a gene expression file.", "error")
                return redirect(url_for("main_routes.create_analysis_page"))

        # Validate gene expression metadata file if provided
        metadata_cols = []
        # 2D Layout File Validation
        if have_2d_layout:
            if not layout_file_2d:
                flash("Please select a 2D layout file.", "error")
                return redirect(url_for("main_routes.create_analysis_page"))

        current_user_node = all_users_data.setdefault(
            current_user.id, {"password": "HASH_PLACEHOLDER", "files": [], "analyses": []}
        )

        if "analyses" not in current_user_node:
            current_user_node["analyses"] = []

        analysis_id = str(uuid.uuid4())
        analysis_data_path = get_user_analysis_path(current_user.id, analysis_id)

        if not analysis_data_path:
            flash("Could not create storage for analysis. Please try again.", "error")
            return redirect(url_for("main_routes.create_analysis_page"))

        new_analysis = {
            "id": analysis_id,
            "name": analysis_name,
            "status": "In Progress",
            "created_at": datetime.now().isoformat(),
            "results_path": analysis_data_path,
            "inputs": {
                "gene_expression": {
                    "source": "h5ad" if have_h5ad else "SEPARATE_FILES",
                    "species": species,
                    **(
                        {
                            "h5ad_filepath": get_file_path(
                                selected_h5ad_file, current_user.id
                            )
                        }
                        if have_h5ad and selected_h5ad_file
                        else {}
                    ),
                    **(
                        {"gene_exp_filepath": get_file_path(gene_exp_file, current_user.id)}
                        if not have_h5ad and gene_exp_file
                        else {}
                    ),
                    **(
                        {"metadata_filepath": get_file_path(metadata_file, current_user.id)}
                        if not have_h5ad and metadata_file and metadata_file != "select-metadata-file"
                        else {}
                    ),
                },
                "data_filtering": data_filtering,
                "layout": {
                    "source": "FILE" if have_2d_layout else "UMAP_GENERATED",
                    **(
                        {"layout_filepath": get_file_path(layout_file_2d, current_user.id)}
                        if have_2d_layout and layout_file_2d
                        else {}
                    ),
                    **(
                        {"umap_settings": umap_parameters}
                        if not have_2d_layout and umap_parameters
                        else {}
                    ),
                }
            },
            **({"metadata_cols": metadata_cols} if not have_h5ad and metadata_file else {}),
        }

        current_user_node["analyses"].append(new_analysis)
        save_all_users_data(all_users_data)

        # 1. Run UMAP Pipeline in the background to generate 2D layout
        run_in_background(
            run_umap_pipeline,
            current_user.id,
            analysis_id,
            new_analysis,
            have_2d_layout,
            update_status_fn=update_analysis_status,
            run_analysis_fn=run_tf_analysis
        )

        flash(f'Analysis "{analysis_name}" created successfully and is pending.', "success")
        return redirect(url_for("main_routes.index"))

    except Exception as e:
        current_app.logger.error(f"Error in create_analysis: {e}")
        flash("An unexpected error occurred while creating analysis.", "error")
        return redirect(url_for("main_routes.create_analysis_page"))


@main_bp.route("/analysis/<analysis_id>")
@login_required
def view_analysis(analysis_id):
    all_users_data = get_all_users_data()
    user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])

    analysis_to_view = find_analysis_by_id(user_analyses, analysis_id)

    if not analysis_to_view:
        flash("Analysis not found.", "error")
        return redirect(url_for("main_routes.index"))

    # Render the view_analysis.html template with the analysis data
    return render_template("view_analysis.html", analysis=analysis_to_view)


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


@main_bp.route("/analysis/plot/<analysis_id>", methods=["POST"])
@login_required
def get_plot(analysis_id):
    all_users_data = get_all_users_data()
    user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
    analysis_to_view = find_analysis_by_id(user_analyses, analysis_id)

    if not analysis_to_view:
        return jsonify({"error": "Analysis not found."}), 404

    if analysis_to_view.get("status") != "Completed":
        return jsonify({"error": "Analysis is not yet complete."}), 400

    return generate_scatter_plot_response(analysis_to_view, plot_type=request.json.get("plot_type"))


def get_layout_and_metadata_dfs(analysis, user_id):
    layout_source = analysis.get("inputs", {}).get("layout", {}).get("source")

    layout_filepath = analysis.get("inputs").get("layout").get("layout_filepath")
    plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))

    if layout_source == "FILE":
        metadata_filepath = analysis.get("inputs").get("gene_expression").get("metadata_filepath", {})
        metadata_df = pd.read_csv(get_file_path(metadata_filepath, user_id), index_col=0, sep=infer_delimiter(metadata_filepath))
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
            return jsonify ({"error": "Layout file not found. Delete this analysis and create new analysis"}), 404

        if not os.path.exists(z_score_filepath):
            current_app.logger.error(f"Z-scores file '{z_score_filepath}' not found.")
            return jsonify ({"error": "Z-scores file not found. Delete this analysis and create new analysis"}), 404

        plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))
        gene_exp_levels_df = pd.read_parquet(z_score_filepath, use_threads=True, columns=[gene_name])

        plot_df = plot_df.merge(gene_exp_levels_df, left_index=True, right_index=True)
        plot_df = plot_df.dropna(subset=[gene_name])

        return plot_df
    except Exception as e:
        current_app.logger.error(f"Error reading layout/z-score file for analysis_id '{analysis['id']}': {e}")
        return  jsonify({"error": "Failed to read layout/z-score file. Invalid format."}), 500


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


@main_bp.route("/analysis/metadata-cluster/<analysis_id>", methods=["POST"])
@login_required
def get_metadata_color_by(analysis_id):
    try:
        all_users_data = get_all_users_data()
        user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
        analysis = find_analysis_by_id(user_analyses, analysis_id)

        if not analysis:
            current_app.logger.error(f"Analysis not found: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis not found."}), 404
        if analysis.get("status") != "Completed":
            current_app.logger.warning(f"Analysis not completed yet: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis is not completed yet."}), 400

        metadata_cluster = request.json.get("selected_metadata_cluster", "").strip()
        plot_type = request.json.get("plot_type", "").strip()

        current_app.logger.debug(
            f"Received metadata_cluster='{metadata_cluster}', plot_type='{plot_type}' for analysis_id={analysis_id}"
        )

        metadata_cols = analysis.get("metadata_cols", [])
        if metadata_cluster not in metadata_cols:
            current_app.logger.warning(f"Requested metadata column '{metadata_cluster}' not found.")
            return jsonify({"error": f"Metadata column '{metadata_cluster}' not found."}), 400

        plot_df, metadata_df = get_layout_and_metadata_dfs(analysis, current_user.id)
        if metadata_cluster not in metadata_df.columns:
            current_app.logger.warning(f"Metadata column '{metadata_cluster}' not found.")
            return (
                jsonify({"error": f"Metadata column '{metadata_cluster}' not found in metadata file."}),
                400,
            )

        # Drop column named "Cluster" if it exists
        if "Cluster" in plot_df.columns:
            plot_df = plot_df.drop(columns=["Cluster"])
        # Merge plot_df with metadata_df on index
        plot_df = plot_df.merge(metadata_df[[metadata_cluster]], left_index=True, right_index=True)
        plot_df = plot_df.rename(columns={metadata_cluster: "Cluster"})
        plot_df["Cluster"] = plot_df["Cluster"].astype(str)

        # Choose columns based on plot_type
        traces, title, x_col, y_col = generate_colored_traces(
            plot_df=plot_df, plot_type=plot_type
        )
        layout = {
            "title": title,
            "xaxis": {"title": x_col},
            "yaxis": {"title": y_col},
        }
        graph_data = {
            "data": traces,
            "layout": layout
        }
        return jsonify(graph_data), 200

    except Exception as e:
        current_app.logger.error(
            f"Error in get_metadata_color_by for analysis_id={analysis_id}.",
            exc_info=True,
        )
        return jsonify({"error": "Internal server error."}), 500


@main_bp.route("/analysis/gene-expression/<analysis_id>", methods=["POST"])
@login_required
def get_gene_expression_color_by(analysis_id):
    try:
        all_users_data = get_all_users_data()
        user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
        analysis = find_analysis_by_id(user_analyses, analysis_id)

        if not analysis:
            current_app.logger.warning(f"Analysis not found: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis not found."}), 404
        if analysis.get("status") != "Completed":
            current_app.logger.info(f"Analysis not completed yet: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis is not completed yet."}), 400

        gene_name = request.json.get("selected_gene", "").strip()
        plot_type = request.json.get("plot_type", "").strip()

        current_app.logger.debug(f"Received gene_name='{gene_name}', plot_type='{plot_type}' for analysis_id={analysis_id}")
        if not gene_name:
            current_app.logger.warning(f"No gene expression provided in request for analysis_id={analysis_id}")
            return jsonify({"error": "Gene expression is required."}), 400

        plot_df = get_layout_and_gene_exp_levels_df(analysis, gene_name)

        title = ""
        x_col = ""
        y_col = ""
        if plot_type.lower() == "umap_plot":
            x_col, y_col = UMAP1_COL, UMAP2_COL
            title = f"UMAP Plot Colored by {gene_name}"
        elif plot_type.lower() == "pca_plot":
            x_col, y_col = PCA1_COL, PCA2_COL
            title = f"PCA Plot Colored by {gene_name}"

        traces = {
            "x": plot_df[x_col].tolist(),
            "y": plot_df[y_col].tolist(),
            "mode": "markers",
            "type": "scattergl",
            "marker": {
                "color": plot_df[gene_name].tolist(),
                "colorscale": "Viridis",
                "showscale": True,
                "colorbar": {"title": gene_name},
            },
        }
        layout = {
            "title": title,
            "xaxis": {"title": x_col},
            "yaxis": {"title": y_col},
        }
        graph_data = {
            "data": [traces],
            "layout": layout,
        }
        return jsonify(graph_data), 200
    except Exception as e:
        current_app.logger.error(
            f"Error in get_gene_expression_color_by for analysis_id={analysis_id}: {e}",
            exc_info=True,
        )
        return jsonify({"error": "Internal server error."}), 500


@main_bp.route("/analysis/change-fdr-tf/<analysis_id>", methods=["POST"])
@login_required
def change_fdr_tf(analysis_id):
    try:
        all_users_data = get_all_users_data()
        user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
        analysis = find_analysis_by_id(user_analyses, analysis_id)

        fdr_level = float(request.json.get("fdr_level", 0.05))
        tf_name = request.json.get("tf_name").strip()
        plot_type = request.json.get("plot_type").strip()

        if not analysis:
            current_app.logger.error(f"Analysis not found: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis not found."}), 404
        if analysis.get("status") != "Completed":
            current_app.logger.warning(f"Analysis not completed yet: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis is not completed yet."}), 400

        pvalues_path = analysis.get("pvalues_path")
        layout_filepath = analysis.get("inputs").get("layout").get("layout_filepath")
        activation_path = analysis.get("activation_path")

        if not os.path.exists(pvalues_path):
            current_app.logger.error(f"P-values file '{pvalues_path}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "P-values file not found. Cannot re-run FDR correction."}), 404
        if not os.path.exists(layout_filepath):
            current_app.logger.error(f"Layout file '{layout_filepath}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "Layout file not found. Cannot re-run FDR correction."}), 404
        if not os.path.exists(activation_path):
            current_app.logger.error(f"Activation file '{activation_path}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "Activation file not found. Cannot re-run FDR correction."}), 404

        plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))
        pvalues_df_tf = pd.read_parquet(pvalues_path, use_threads=True, columns=[tf_name])
        activation_tf = pd.read_parquet(activation_path, use_threads=True, columns=[tf_name])

        corrected_pvalues, reject, p_value_thresholds = bh_fdr_correction(pvalues_df_tf, fdr_level)

        final_output = pd.DataFrame(0, index=reject.index, columns=reject.columns)
        final_output[reject] = activation_tf[reject]
        final_output[pvalues_df_tf.isna()] = np.nan
        final_output = final_output.astype("Int64")

        plot_df = plot_df.merge(final_output.astype(object), left_index=True, right_index=True)
        traces, title, x_col, y_col = generate_colored_traces(
            plot_df=plot_df,
            plot_type=plot_type,
            tf_activity=tf_name,
        )
        layout = {
            "title": title,
            "xaxis": {"title": x_col},
            "yaxis": {"title": y_col}
        }
        graph_data = {
            "data": traces,
            "layout": layout,
            "fdr_level": fdr_level,
            "p_value_threshold": p_value_thresholds.iloc[0]
        }
        return jsonify(graph_data), 200

    except Exception as e:
        current_app.logger.error(
            f"Error re-running FDR correction for analysis_id={analysis_id}: {e}",
            exc_info=True,
        )
        return jsonify({"error": "Internal server error."}), 500


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


@main_bp.route("/analysis/change-pvalue-threshold-tf/<analysis_id>", methods=["POST"])
@login_required
def change_pvalue_threshold_tf(analysis_id):
    try:
        all_users_data = get_all_users_data()
        user_analyses = all_users_data.get(current_user.id, {}).get("analyses", [])
        analysis = find_analysis_by_id(user_analyses, analysis_id)

        fdr_level = float(request.json.get("fdr_level"))
        pvalue_threshold = float(request.json.get("pvalue_threshold"))
        tf_name = request.json.get("tf_name").strip()
        plot_type = request.json.get("plot_type").strip()

        if not analysis:
            current_app.logger.error(f"Analysis not found: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis not found."}), 404
        if analysis.get("status") != "Completed":
            current_app.logger.warning(f"Analysis not completed yet: analysis_id={analysis_id}")
            return jsonify({"error": "Analysis is not completed yet."}), 400

        pvalues_path = analysis.get("pvalues_path")
        layout_filepath = analysis.get("inputs").get("layout").get("layout_filepath")
        activation_path = analysis.get("activation_path")

        if not os.path.exists(pvalues_path):
            current_app.logger.error(f"P-values file '{pvalues_path}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "P-values file not found. Cannot re-run FDR correction."}), 404
        if not os.path.exists(layout_filepath):
            current_app.logger.error(f"Layout file '{layout_filepath}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "Layout file not found. Cannot re-run FDR correction."}), 404
        if not os.path.exists(activation_path):
            current_app.logger.error(f"Activation file '{activation_path}' not found for analysis_id '{analysis_id}'.")
            return jsonify({"error": "Activation file not found. Cannot re-run FDR correction."}), 404

        plot_df = pd.read_csv(layout_filepath, index_col=0, sep=infer_delimiter(layout_filepath))
        pvalues_df = pd.read_parquet(pvalues_path, use_threads=True)
        activation_tf = pd.read_parquet(activation_path, use_threads=True, columns=[tf_name])
        fdr = estimate_fdr_for_gene(pvalues_df, tf_name, pvalue_threshold)

        pvalues_df_tf = pvalues_df[[tf_name]].dropna()
        reject_mask = pvalues_df_tf < pvalue_threshold
        reject = reject_mask.where(pvalues_df_tf.notna())
        reject = reject.astype("boolean")

        final_output = pd.DataFrame(0, index=reject.index, columns=reject.columns)
        final_output[reject] = activation_tf[reject]
        final_output[pvalues_df_tf.isna()] = np.nan
        final_output = final_output.astype("Int64")

        plot_df = plot_df.merge(final_output.astype(object), left_index=True, right_index=True)
        traces, title, x_col, y_col = generate_colored_traces(
            plot_df=plot_df,
            plot_type=plot_type,
            tf_activity=tf_name,
        )
        layout = {
            "title": title,
            "xaxis": {"title": x_col},
            "yaxis": {"title": y_col}
        }
        graph_data = {
            "data": traces,
            "layout": layout,
            "fdr_level": fdr_level,
            "p_value_threshold": pvalue_threshold
        }
        return jsonify(graph_data), 200

    except Exception as e:
        current_app.logger.error(
            f"Error re-running FDR correction for analysis_id={analysis_id}: {e}",
            exc_info=True,
        )
        return jsonify({"error": "Internal server error."}), 500


@main_bp.route("/analysis/delete/<analysis_id>", methods=["POST"])
@login_required
def delete_analysis(analysis_id):
    all_users_data = get_all_users_data()
    current_user_node = all_users_data.get(current_user.id)

    if not current_user_node or "analyses" not in current_user_node:
        flash("No analyses found for this user.", "error")
        return redirect(url_for("main_routes.index"))

    analysis_to_delete = None
    analysis_index = -1
    for i, analysis in enumerate(current_user_node["analyses"]):
        if analysis["id"] == analysis_id:
            analysis_to_delete = analysis
            analysis_index = i
            break

    if analysis_to_delete:
        # Attempt to delete the analysis directory from filesystem
        analysis_data_path = analysis_to_delete.get("results_path")
        if analysis_data_path and os.path.exists(analysis_data_path):
            try:
                shutil.rmtree(
                    analysis_data_path
                )  # Deletes the directory and all its contents
                current_app.logger.info(
                    f"Deleted analysis directory: {analysis_data_path}"
                )
            except OSError as e:
                current_app.logger.error(
                    f"Error deleting analysis directory {analysis_data_path}: {e}"
                )
                flash(
                    f"Error deleting analysis files from disk: {e}. Record removed.",
                    "warning",
                )

        # Remove the analysis record from users.json
        current_user_node["analyses"].pop(analysis_index)
        save_all_users_data(all_users_data)
        flash(
            f'Analysis "{analysis_to_delete["name"]}" deleted successfully.', "success"
        )
    else:
        flash("Analysis not found or already deleted.", "error")

    return redirect(url_for("main_routes.index"))


# --- Auth Routes (login, signup, logout) ---
@auth_bp.route("/signup", methods=["GET", "POST"])
def signup():
    if current_user.is_authenticated:
        return redirect(url_for("main_routes.index"))

    if request.method == "POST":
        username = request.form.get("username", "").strip()
        password = request.form.get("password")
        confirm_password = request.form.get("confirm_password")

        if not username or not password or not confirm_password:
            flash("All fields are required.", "error")
            return render_template("signup.html")  # Re-render form with error

        if " " in username:  # Example validation
            flash("Username cannot contain spaces.", "error")
            return render_template("signup.html")

        if password != confirm_password:
            flash("Passwords do not match.", "error")
            return render_template("signup.html")

        all_users_data = get_all_users_data()
        if username in all_users_data:
            flash("Username already exists. Please choose a different one.", "error")
            return render_template("signup.html")

        hashed_password = generate_password_hash(password)
        all_users_data[username] = {
            "password": hashed_password,
            "files": [],
            "analyses": [],
        }
        save_all_users_data(all_users_data)
        get_user_specific_data_path(username)

        flash("Account created successfully! Please login.", "success")
        return redirect(url_for("auth.login"))

    return render_template("signup.html")


@auth_bp.route("/login", methods=["GET", "POST"])
def login():
    if current_user.is_authenticated:
        return redirect(url_for("main_routes.index"))

    if request.method == "POST":
        username = request.form.get("username", "").strip()
        password = request.form.get("password")

        if not username or not password:
            flash("Both username and password are required.", "error")
            return render_template("login.html")

        all_users_data = get_all_users_data()
        user_account_data = all_users_data.get(username)

        if user_account_data and check_password_hash(
            user_account_data.get("password", ""), password
        ):
            user_obj = User(username)  # Create User instance
            login_user(user_obj)  # Flask-Login handles session
            flash("Logged in successfully!", "success")
            next_page = request.args.get("next")
            return redirect(next_page or url_for("main_routes.index"))
        else:
            flash("Invalid username or password.", "error")
            return render_template("login.html")

    return render_template("login.html")


@auth_bp.route("/logout")
@login_required  # Ensure user is logged in to log out
def logout():
    logout_user()  # Clears user session
    flash("You have been logged out.", "info")
    return redirect(url_for("auth.login"))
