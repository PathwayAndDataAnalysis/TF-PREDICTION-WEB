from app import create_app  # Import the app factory from the 'app' package

# Create the Flask app instance using the factory
# This allows configuration to be passed if needed, e.g., for different environments
application = create_app()  # Renamed to 'application' for common WSGI server convention

if __name__ == "__main__":
    # For development, Flask's built-in server is fine.
    # For production, use a WSGI server like Gunicorn or uWSGI.
    # Example: gunicorn -w 4 "main:application"
    # The 'application' in "main:application" refers to the app instance created above.

    # The host '0.0.0.0' makes the server accessible from other devices on the network.
    # debug=True should NOT be used in production.
    application.run(host="0.0.0.0", port=5000, debug=True)
