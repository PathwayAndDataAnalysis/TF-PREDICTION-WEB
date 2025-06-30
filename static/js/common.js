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
