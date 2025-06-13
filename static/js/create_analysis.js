const layoutCheckbox = document.getElementById("have_2d_layout");
const layoutSelectionSection = document.getElementById("2d_layout_selection_section");
const umapParamsSection = document.getElementById("umap_parameters_section");
const layoutFile2D = document.getElementById("layout_file_2d");
const h5adCheckbox = document.getElementById("have_h5ad");
const h5adSelectionSection = document.getElementById("h5ad_selection_section");
const h5adRadios = document.getElementsByName("selected_h5ad_file");
const separateFilesSection = document.getElementById("separate_files_section");
const geneExpFile = document.getElementById("gene_exp_file");
const metadataFile = document.getElementById("metadata_file");

function toggleH5adSection() {
	if (h5adCheckbox.checked) {
		h5adSelectionSection.style.display = "block";
		separateFilesSection.style.display = "none";
		geneExpFile.required = false;
		metadataFile.required = false;
	} else {
		h5adSelectionSection.style.display = "none";
		separateFilesSection.style.display = "block";
		geneExpFile.required = true;
		metadataFile.required = false;

		for (let radio of h5adRadios) {
			radio.checked = false;
		}
	}
}

function toggleLayoutSection() {
	if (layoutCheckbox.checked) {
		layoutSelectionSection.style.display = "block";
		umapParamsSection.style.display = "none";
		layoutFile2D.required = true;
	} else {
		layoutSelectionSection.style.display = "none";
		umapParamsSection.style.display = "block";
		layoutFile2D.required = false;
		layoutFile2D.value = "";
	}
}

document.addEventListener("DOMContentLoaded", function () {
	toggleH5adSection();
	toggleLayoutSection();
});
