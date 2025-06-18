from flask import jsonify
from app import create_app
import traceback

application = create_app()

# Global error handlers
@application.errorhandler(413)
def too_large(e):
    return jsonify(error="File too large"), 413

@application.errorhandler(500)
def internal_error(e):
    application.logger.error(f"Internal server error: {e}")
    application.logger.error(traceback.format_exc())
    return jsonify(error="Internal server error"), 500

@application.errorhandler(404)
def not_found(e):
    return jsonify(error="Resource not found"), 404

@application.errorhandler(Exception)
def handle_exception(e):
    application.logger.error(f"Unhandled exception: {e}")
    application.logger.error(traceback.format_exc())
    return jsonify(error="An unexpected error occurred"), 500

if __name__ == "__main__":
    # The host '0.0.0.0' makes the server accessible from other devices on the network.
    # debug=True should NOT be used in production.
    application.run(host="0.0.0.0", port=5000, debug=True)
