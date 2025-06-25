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

	const dataFileInput = document.getElementById('data_file');
    const geneExpCheckboxContainer = document.getElementById('gene-exp-checkbox-container');

    if (dataFileInput && geneExpCheckboxContainer) {
        dataFileInput.addEventListener('change', () => {
            const file = dataFileInput.files[0];
            if (!file) {
                geneExpCheckboxContainer.classList.add('hidden');
                return;
            }

            const filename = file.name.toLowerCase();
            // Show the checkbox for tsv or csv files
            if (filename.endsWith('.tsv') || filename.endsWith('.csv')) {
                geneExpCheckboxContainer.classList.remove('hidden');
            } else {
                geneExpCheckboxContainer.classList.add('hidden');
            }
        });
    }
});
