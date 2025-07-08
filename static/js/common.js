document.addEventListener("DOMContentLoaded", function () {
	const flashMessages = document.querySelectorAll("#flash-message-container .flash-message-item");
	const FADE_DELAY = 3000; // 3 seconds
	const REMOVE_DELAY_AFTER_FADE = 500;

	flashMessages.forEach(function (message, index) {
		setTimeout(function () {
			message.classList.add("fade-out");
			setTimeout(function () {
				if (message.parentNode) {
					message.remove();
				}
			}, REMOVE_DELAY_AFTER_FADE);
		}, FADE_DELAY + index * 300);
	});

});

function addInputValidation(inputId, min, max=null) {
    const input = document.getElementById(inputId);
    if (!input) return;
    input.addEventListener("input", function () {
        let value = parseInt(this.value, 10);

        if (isNaN(value)) {
            this.value = "";
            return;
        }
        if (value < min) this.value = min;
		if (max !== null && value > max) {
			this.value = max;
		}
    });
}
