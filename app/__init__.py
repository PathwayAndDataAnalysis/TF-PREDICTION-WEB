import json
import os
import signal
from contextlib import contextmanager

from flask import Flask, current_app
from flask_login import LoginManager, UserMixin


# --- User Model ---
class User(UserMixin):
    def __init__(self, id_):
        self.id = id_


# --- Login Manager Setup ---
login_manager = LoginManager()
login_manager.login_view = "auth.login"  # BlueprintName.routeName
login_manager.login_message_category = "info"
login_manager.login_message = "Please log in to access this page."


# --- User Data Helpers (using current_app for context) ---
def get_users_file_path():
    return current_app.config["USERS_FILE"]


def get_upload_folder_root():
    return current_app.config["UPLOAD_FOLDER"]


def get_all_users_data():
    """Loads all users' data from the JSON file."""
    try:
        with open(get_users_file_path(), "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError:
        current_app.logger.error(f"Error decoding JSON from {get_users_file_path()}")
        # Potentially re-initialize or backup the corrupted file
        return {}  # Return empty dict to prevent app crash, but log issue


def save_all_users_data(users_data):
    """Saves all users' data to the JSON file."""
    try:
        with open(get_users_file_path(), "w") as f:
            json.dump(users_data, f, indent=4)
    except IOError as e:
        current_app.logger.error(f"Error saving users to {get_users_file_path()}: {e}")


def find_analysis_by_id(analyses, analysis_id):
    """Finds an analysis by its ID in a list of analyses."""
    return next((a for a in analyses if a.get("id") == analysis_id), None)


def get_file_path(filename, user_id=None):
    """
    Returns the full path to a file in the uploads folder.
    If user_id is provided, returns the path in the user's folder.
    Otherwise, returns the path in the global uploads folder.
    If the file does not exist, returns None.

    Returns:
        str or None: The full file path if the file exists, otherwise None.
    """
    upload_folder = get_upload_folder_root()
    if user_id:
        path = os.path.join(upload_folder, user_id, filename)
    else:
        path = os.path.join(upload_folder, filename)
    if os.path.exists(path):
        return path
    else:
        raise FileNotFoundError(f"File {path} does not exist")


@login_manager.user_loader
def load_user(user_id):
    """Flask-Login hook to load a user by ID."""
    users_data = get_all_users_data()
    if user_id in users_data:
        return User(user_id)
    return None


def create_app(test_config=None):
    # Create and configure the app
    # app = Flask(__name__, instance_relative_config=True) # App name is 'app'
    app = Flask(
        __name__,
        instance_relative_config=True,
        template_folder="../templates",
        static_folder="../static",
    )

    # --- Configuration ---
    app.config.from_mapping(
        SECRET_KEY="dev_secret_key_CHANGE_THIS_IN_PRODUCTION!",  # IMPORTANT!
        # USERS_FILE will be in the instance folder
        USERS_FILE=os.path.join(app.instance_path, "users.json"),
        # UPLOAD_FOLDER will be at the project root level, sibling to 'app' and 'main.py'
        UPLOAD_FOLDER=os.path.join(os.path.dirname(app.root_path), "user_uploads"),
        ALLOWED_EXTENSIONS={
            "h5ad",
            "txt",
            "csv",
            "tsv",
            "gz",
            "zip",
            "rds",
        },  # Added rds
        MAX_CONTENT_LENGTH=500
        * 1024
        * 1024
        * 1024,  # Example: 500 GB limit for uploads
    )

    if test_config is None:
        # Load the instance config, if it exists, when not testing
        # e.g., create tf-pred-webserver/instance/config.py for production secrets
        app.config.from_pyfile("config.py", silent=True)
    else:
        # Load the test config if passed in
        app.config.from_mapping(test_config)

    # Ensure the instance folder exists (for USERS_FILE, config.py)
    try:
        os.makedirs(app.instance_path, exist_ok=True)
    except OSError as e:
        app.logger.error(f"Could not create instance folder: {app.instance_path} - {e}")

    # Ensure UPLOAD_FOLDER exists
    upload_folder_path = app.config["UPLOAD_FOLDER"]
    if not os.path.exists(upload_folder_path):
        try:
            os.makedirs(upload_folder_path, exist_ok=True)
            app.logger.info(f"Upload folder created at {upload_folder_path}")
        except OSError as e:
            app.logger.error(
                f"Could not create upload folder: {upload_folder_path} - {e}"
            )

    # Initialize Flask extensions
    login_manager.init_app(app)

    # Register Blueprints
    from . import routes as app_routes  # Using an alias to avoid confusion

    app.register_blueprint(app_routes.auth_bp)
    app.register_blueprint(app_routes.main_bp)

    # Initialize users.json if it doesn't exist
    # This should be done within app context to access app.config and logger
    with app.app_context():
        users_fpath = get_users_file_path()
        if not os.path.exists(users_fpath):
            save_all_users_data({})  # Create an empty users file
            current_app.logger.info(f"Initialized empty users file at {users_fpath}")

    return app