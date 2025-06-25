const layoutCheckbox = document.getElementById("have_2d_layout");
const layoutSelectionSection = document.getElementById("2d_layout_selection_section");
const umapParamsSection = document.getElementById("umap_parameters_section");
const layoutFile2D = document.getElementById("layout_file_2d");
const h5adCheckbox = document.getElementById("have_h5ad");
const h5adSelectionSection = document.getElementById("h5ad_selection_section");
const h5adFileSelection = document.getElementById("selected_h5ad_file");
const separateFilesSection = document.getElementById("separate_files_section");
const geneExpFile = document.getElementById("gene_exp_file");
const speciesSelection = document.getElementById("species");
const metadataFile = document.getElementById("metadata_file");

function toggleH5adSection() {
	if (h5adCheckbox.checked) {
		h5adSelectionSection.style.display = "block";
		separateFilesSection.style.display = "none";
		geneExpFile.required = false;
		metadataFile.required = false;

		geneExpFile.value = "select-gene-exp-file";
		speciesSelection.value = "select-species";
		metadataFile.value = "select-metadata-file";
	} else {
		h5adSelectionSection.style.display = "none";
		separateFilesSection.style.display = "block";
		geneExpFile.required = true;
		metadataFile.required = false;

		h5adFileSelection.value = "select-h5ad-file";
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

function renderHistogram(divElement, data, title, xtitle, ytitle) {
	if (!data || !data.bins || !data.counts) {
		divElement.innerHTML = `<div class="text-center text-gray-500 p-4">Data not available for ${title}.</div>`;
		return;
	}

	// Plotly expects bin centers, not edges, for bar charts representing histograms
	const binCenters = data.bins.slice(0, -1).map((b, i) => (b + data.bins[i + 1]) / 2);

	const plotData = [
		{
			x: binCenters,
			y: data.counts,
			type: "bar",
			marker: { color: "#3b82f6" }, // A nice blue color
		},
	];

	const layout = {
		title: { text: title, font: { size: 14 } },
		xaxis: { title: { text: xtitle, font: { size: 12 } }, automargin: true },
		yaxis: { title: { text: ytitle, font: { size: 12 } }, automargin: true },
		margin: { t: 40, b: 50, l: 50, r: 20 },
		bargap: 0.05,
		dragmode: "pan",
		hovermode: "closest",
	};

	Plotly.newPlot(divElement, plotData, layout, {
		responsive: true,
		displayModeBar: false,
		displaylogo: false,
		scrollZoom: true,
	});
}

document.addEventListener("DOMContentLoaded", function () {
	toggleH5adSection();
	toggleLayoutSection();

	const plotsContainer = document.getElementById("qc-plots-container");
	const plotDivs = {
		n_genes_by_counts: document.getElementById("qc-plot-genes-per-cell"),
		n_cells_by_counts: document.getElementById("qc-plot-cells-per-gene"),
		pct_counts_mt: document.getElementById("qc-plot-mt-percent"),
	};

	function showQcPlotsForFile(selectedFilename) {
		if (!selectedFilename || selectedFilename.startsWith("select-")) {
			plotsContainer.classList.add("hidden");
			return;
		}
		for (let file of user_files) {
			if (selectedFilename === file.filename) {
				// Check if qc_metrics exists and is an object
				if (!file.qc_metrics || typeof file.qc_metrics !== "object") {
					plotsContainer.classList.add("hidden");
					return;
				}

				plotsContainer.classList.remove("hidden");
				renderHistogram(
					plotDivs.n_genes_by_counts,
					file.qc_metrics.n_genes_by_counts,
					"Genes per Cell",
					"Number of Genes",
					"Cell Count"
				);
				renderHistogram(
					plotDivs.n_cells_by_counts,
					file.qc_metrics.n_cells_by_counts,
					"Cells per Gene",
					"Number of Cells",
					"Gene Count"
				);
				renderHistogram(
					plotDivs.pct_counts_mt,
					file.qc_metrics.pct_counts_mt,
					"Mitochondrial Content %",
					"% MT",
					"Cell Count"
				);

				document.getElementById("qc-plot-genes-per-cell-stats").innerHTML = `<p class="text-sm mt-2">
						<span class="font-semibold">Mean: </span>${file.qc_metrics.data_summary.mean_n_genes_by_counts}, 
						<span class="font-semibold">Standard Deviation: </span>${file.qc_metrics.data_summary.sd_n_genes_by_counts}, 
						<span class="font-semibold">Min. Value: </span>${file.qc_metrics.data_summary.min_n_genes_by_counts}, 
						<span class="font-semibold">Max. Value: </span>${file.qc_metrics.data_summary.max_n_genes_by_counts}
					</p>`;
				document.getElementById("qc-plot-cells-per-gene-stats").innerHTML = `<p class="text-sm mt-2">
						<span class="font-semibold">Mean: </span>${file.qc_metrics.data_summary.mean_n_cells_by_counts}, 
						<span class="font-semibold">Standard Deviation: </span>${file.qc_metrics.data_summary.sd_n_cells_by_counts}, 
						<span class="font-semibold">Min. Value: </span>${file.qc_metrics.data_summary.min_n_cells_by_counts}, 
						<span class="font-semibold">Max. Value: </span>${file.qc_metrics.data_summary.max_n_cells_by_counts}
					</p>`;
				document.getElementById("qc-plot-mt-percent-stats").innerHTML = `<p class="text-sm mt-2">
						<span class="font-semibold">Mean: </span>${file.qc_metrics.data_summary.mean_pct_counts_mt}, 
						<span class="font-semibold">Standard Deviation: </span>${file.qc_metrics.data_summary.sd_pct_counts_mt}, 
						<span class="font-semibold">Min. Value: </span>${file.qc_metrics.data_summary.min_pct_counts_mt}, 
						<span class="font-semibold">Max. Value: </span>${file.qc_metrics.data_summary.max_pct_counts_mt}
					</p>`;
				return;
			}
		}
		plotsContainer.classList.add("hidden");
	}

	geneExpFile.addEventListener("change", function () {
		showQcPlotsForFile(geneExpFile.value);
	});

	h5adFileSelection.addEventListener("change", function () {
		showQcPlotsForFile(h5adFileSelection.value);
	});
});
