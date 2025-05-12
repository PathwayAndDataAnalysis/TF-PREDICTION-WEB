// JS for create_analysis.html
// Add create_analysis-specific JS here, e.g., form validation, dynamic UI, etc.

function toggleH5adSection() {
	const h5adCheckbox = document.getElementById("have_h5ad");
	const h5adSelectionSection = document.getElementById("h5ad_selection_section");
	const separateFilesSection = document.getElementById("separate_files_section");
	const rawCountsFile = document.getElementById("raw_counts_file");
	const metadataFile = document.getElementById("metadata_file");

	if (h5adCheckbox.checked) {
		h5adSelectionSection.style.display = "block";
		separateFilesSection.style.display = "none";
		rawCountsFile.required = false;
		metadataFile.required = false;
	} else {
		h5adSelectionSection.style.display = "none";
		separateFilesSection.style.display = "block";
		rawCountsFile.required = true;
		metadataFile.required = false;
		const h5adRadios = document.getElementsByName("selected_h5ad_file");
		for (let radio of h5adRadios) {
			radio.checked = false;
		}
	}
}

function toggleLayoutSection() {
	const layoutCheckbox = document.getElementById("have_2d_layout");
	const layoutSelectionSection = document.getElementById("2d_layout_selection_section");
	const umapParamsSection = document.getElementById("umap_parameters_section");
	const layoutFile = document.getElementById("layout_file");

	if (layoutCheckbox.checked) {
		layoutSelectionSection.style.display = "block";
		umapParamsSection.style.display = "none";
		layoutFile.required = true;
	} else {
		layoutSelectionSection.style.display = "none";
		umapParamsSection.style.display = "block";
		layoutFile.required = false;
		layoutFile.value = "";
	}
}

document.addEventListener("DOMContentLoaded", function () {
	toggleH5adSection();
	toggleLayoutSection();
});
