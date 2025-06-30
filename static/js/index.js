document.addEventListener("DOMContentLoaded", function() {

	const dataFileInput = document.getElementById('data_file');
    const fileTypeContainer = document.getElementById('file-type-container');
    const uploadModal = document.getElementById('upload-modal');
    const openUploadModalBtn = document.getElementById('open-upload-modal-btn');
    const closeUploadModalBtn = document.getElementById('close-upload-modal-btn');
    const cancelUploadBtn = document.getElementById('cancel-upload-btn');
     const seeMoreButtons = document.querySelectorAll('.see-more-btn');
     const searchInput = document.getElementById('analysis-search-input');
    const cardsContainer = document.getElementById('analysis-cards-container');

    if (dataFileInput && fileTypeContainer) {
        dataFileInput.addEventListener('change', () => {
            const file = dataFileInput.files[0];
            if (!file) {
                fileTypeContainer.classList.add('hidden');
                return;
            }

            const filename = file.name.toLowerCase();
            // Show the checkbox for tsv or csv files
            if (filename.endsWith('.tsv') || filename.endsWith('.csv')) {
                fileTypeContainer.classList.remove('hidden');
            } else {
                fileTypeContainer.classList.add('hidden');
            }
        });
    }

    // Check if we are on a page that has the modal elements
    if (!uploadModal || !openUploadModalBtn || !closeUploadModalBtn || !cancelUploadBtn) {
        return; // Exit if modal elements aren't found
    }

    const openModal = () => {
        uploadModal.classList.remove('hidden');
        uploadModal.classList.add('flex'); // Use flex to enable centering
    };

    const closeModal = () => {
        uploadModal.classList.add('hidden');
        uploadModal.classList.remove('flex');
    };

    openUploadModalBtn.addEventListener('click', openModal);
    closeUploadModalBtn.addEventListener('click', closeModal);
    cancelUploadBtn.addEventListener('click', closeModal);

    // Close the modal if the user clicks on the dark background
    uploadModal.addEventListener('click', (event) => {
        if (event.target === uploadModal) {
            closeModal();
        }
    });

    // Close the modal if the user presses the 'Escape' key
    document.addEventListener('keydown', (event) => {
        if (event.key === 'Escape' && !uploadModal.classList.contains('hidden')) {
            closeModal();
        }
    });

    seeMoreButtons.forEach(button => {
        button.addEventListener('click', (event) => {
            const parentDiv = event.target.closest('div');
            const shortSpan = parentDiv.querySelector('.description-short');
            const fullSpan = parentDiv.querySelector('.description-full');

            // Toggle visibility
            shortSpan.classList.toggle('hidden');
            fullSpan.classList.toggle('hidden');

            // Change button text
            if (fullSpan.classList.contains('hidden')) {
                event.target.textContent = 'See more';
            } else {
                event.target.textContent = 'See less';
            }
        });
    });

    if (searchInput && cardsContainer) {
        searchInput.addEventListener('input', (event) => {
            const searchTerm = event.target.value.toLowerCase();
            const cards = cardsContainer.querySelectorAll('.analysis-card');

            cards.forEach(card => {
                const cardText = card.textContent.toLowerCase();
                if (cardText.includes(searchTerm)) {
                    card.style.display = 'flex'; // Use 'flex' since our card uses it
                } else {
                    card.style.display = 'none';
                }
            });
        });
    }

})