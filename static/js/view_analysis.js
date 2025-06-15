// Use window.analysis for the analysis object
console.log("Analysis object:", window.analysis);

const defaultPointSize = 4;
const defaultOpacity = 0.5;

// Plot configuration panel logic
const plotTypeSelect = document.getElementById("plot-type");
const colorBySelect = document.getElementById("color-by");
const metadataColSelectionDiv = document.getElementById("metadata-column-selection-div");
const metadataColNameSelect = document.getElementById("metadata-column-name");
const tfSelectionDiv = document.getElementById("tf-selection-div");
const tfManualEntryDiv = document.getElementById("tf-manual-entry-div");
const geneEntryDiv = document.getElementById("gene-entry-div");
const tfNameSelect = document.getElementById("tf-name-select");
const tfManualEntryInput = document.getElementById("tf-manual-entry-input");
const geneEntryInput = document.getElementById("gene-entry-input");
const pointSizeSlider = document.getElementById("point-size");
const pointSizeValue = document.getElementById("point-size-value");
const opacitySlider = document.getElementById("opacity");
const opacityValue = document.getElementById("opacity-value");
const showLegendCheckbox = document.getElementById("show-legend");

plotTypeSelect.addEventListener("change", function () {
	if (this.value === "umap_plot") {  // UMAP Plot
		console.log("UMAP plot selected");
		const apiUrl = `/analysis/umap_plot/${window.analysis.id}`;
		getPlotData(this.value, apiUrl)
		.then(() => {
			console.log("UMAP Plot loaded successfully.");
		})
		.catch((err) => {
			console.error("Failed to load plot:", err);
			alert("Failed to load plot. Please try again later.");
		});
	}
	else {  // PCA Plot
		console.log("PCA plot selected");
		const apiUrl = `/analysis/pca_plot/${window.analysis.id}`;
		getPlotData(this.value, apiUrl)
		.then(() => {
			console.log("PCA Plot loaded successfully.");

			colorBySelect.value = "select_cluster_type";
			tfNameSelect.value = "select_tf";
			metadataColSelectionDiv.value = "select_metadata_column";
			colorBySelect.dispatchEvent(new Event("change"));
			tfSelectionDiv.dispatchEvent(new Event("change"));
			metadataColSelectionDiv.dispatchEvent(new Event("change"));

		})
		.catch((err) => {
			console.error("Failed to load plot:", err);
			alert("Failed to load plot. Please try again later.");
		});
	}
});

if (colorBySelect) {
	colorBySelect.addEventListener("change", function () {
		if (this.value === "tf_activity") {
			tfSelectionDiv.classList.remove("hidden");
			tfManualEntryDiv.classList.remove("hidden");
			metadataColSelectionDiv.classList.add("hidden");
			geneEntryDiv.classList.add("hidden");
		}
		else if (this.value === "metadata_columns") {
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.remove("hidden");
			geneEntryDiv.classList.add("hidden");
		}
		else if (this.value === "gene_expression") {
			geneEntryDiv.classList.remove("hidden");
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.add("hidden");
		}
		else {
			geneEntryDiv.classList.remove("hidden");
			geneEntryDiv.classList.add("hidden");
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.add("hidden");
		}
	});
}

function updatePlot(plot_data){
	if (plot_data.data && plot_data.layout) {
		plot_data.layout.dragmode = "pan"; // Set default to pan
		plot_data.layout.hovermode = "closest"; // Set default hover mode
		plot_data.layout.showlegend = showLegendCheckbox.checked;
		plot_data.layout.legend = {
			x: 1,        			// Position legend slightly to the right of the plot area (0-1 is plot area)
			y: 1,            		// Align legend top with plot top
			xanchor: 'right',   	// Anchor the legend's left edge to the x position
			yanchor: 'top',    		// Anchor the legend's top edge to the y position
			traceorder: 'normal', 	// or 'reversed' or 'grouped'
		};
		plot_data.layout.margin = {
			l: 15, // left
			r: 15, // right
			t: 30, // top
			b: 15  // bottom
		}
		plot_data.layout.autosize = true
		plot_data.layout.xaxis ={
			title: plot_data.layout.xaxis.title,
		}
		plot_data.layout.yaxis = {
			title: plot_data.layout.yaxis.title,
		}
		plot_data.data.forEach(trace => {
			trace.marker = trace.marker || {};
			trace.marker.size = pointSizeSlider ? Number(pointSizeSlider.value) : defaultPointSize;
			trace.marker.opacity = opacitySlider ? Number(opacitySlider.value) : defaultOpacity;
		});

		Plotly.newPlot("scatterPlot", plot_data.data, plot_data.layout, {
			responsive: true,
			displayModeBar: true,
			displaylogo: false,
			scrollZoom: true
		});

	} else {
		throw new Error("Received data is not in the expected format.");
	}
}

metadataColNameSelect.addEventListener("change", function () {
	if (this.value !== "select_metadata_column") {
		const apiUrl = `/analysis/metadata-cluster/${window.analysis.id}`;

		getPlotData(
			this.value, apiUrl,
			"POST",
			{ selected_metadata_cluster: this.value, plot_type: plotTypeSelect.value }
		).then(() => {
			console.log("Metadata cluster plot loaded successfully.");
		})
		.catch((err) => {
			console.error("Failed to load plot:", err);
			alert("Failed to load plot. Please try again later.");
		});
	}
});

tfNameSelect.addEventListener("change", function () {
	if (this.value !== "select_tf"){
		const apiUrl = `/analysis/tf-activity/${window.analysis.id}`;

		getPlotData(
			this.value, apiUrl,
			"POST",
			{ selected_tf: this.value, plot_type: plotTypeSelect.value }
		).then(() => {
			console.log("TF activity plot loaded successfully.");
		})
		.catch((err) => {
			console.error("Failed to load plot:", err);
			alert("Failed to load plot. Please try again later.");
		});
	}
});

tfManualEntryInput.addEventListener("keypress", function (event) {
	if (event.key === "Enter") {
		let tf_name = this.value.trim().toUpperCase();
		if (tf_name) {
			const apiUrl = `/analysis/tf-activity/${window.analysis.id}`;

			getPlotData(
				tf_name, apiUrl,
				"POST",
				{selected_tf: tf_name, plot_type: plotTypeSelect.value}
			).then(() => {
				console.log("TF activity plot loaded successfully.");
			})
				.catch((err) => {
					console.error("Failed to load plot:", err);
					alert("Failed to load plot. Please try again later.");
				});
		} else {
			alert("Please enter a valid TF name.");
		}
	}
});

geneEntryInput.addEventListener("keypress", function (event) {
	if (event.key === "Enter") {
		let gene_name = this.value.trim().toUpperCase();
		if (gene_name) {
			const apiUrl = `/analysis/gene-expression/${window.analysis.id}`;

			getPlotData(
				gene_name, apiUrl,
				"POST",
				{ selected_gene: gene_name, plot_type: plotTypeSelect.value }
			).then(() => {
				console.log("Gene expression plot loaded successfully.");
			})
				.catch((err) => {
					console.error("Failed to load plot:", err);
					alert("Failed to load plot. Please try again later.");
				});
		} else {
			alert("Please enter a valid gene name.");
		}
	}
})

if (pointSizeSlider && pointSizeValue) {
	pointSizeSlider.addEventListener("input", function (e) {
		Plotly.restyle("scatterPlot", { "marker.size": Number(e.target.value) });
		pointSizeValue.textContent = e.target.value;
	});
}

if (opacitySlider && opacityValue) {
	opacitySlider.addEventListener("input", function (e) {
		Plotly.restyle("scatterPlot", { "marker.opacity": Number(e.target.value) });
		opacityValue.textContent = e.target.value;
	});
}

if (showLegendCheckbox) {
	showLegendCheckbox.addEventListener("change", function (e) {
		Plotly.relayout("scatterPlot", { showlegend: e.target.checked });
	});
}

const plotConfigForm = document.getElementById("plot-config-form");
if (plotConfigForm) {
	plotConfigForm.addEventListener("submit", function (event) {
		event.preventDefault();
		// Get form data
		const formData = new FormData(this);
		const config = Object.fromEntries(formData.entries());
		console.log("Plot configuration updated:", config);
		// Call a function to re-render the plot with new config
		// e.g., updateUmapPlot(config);
		// alert("Plot update functionality not yet implemented.");
	});
}

// UMAP plot loading and rendering
async function getPlotData(plot_type = "umap_plot", apiUrl, method = "GET", body = null) {
	document.getElementById("plot-loading-spinner").style.display = "block";

	try {
		const response = await fetch(apiUrl, {
			method: method,
			headers: { "Content-Type": "application/json" },
			body: body ? JSON.stringify(body) : null,
		});

		let data;
		try {
			data = await response.json();
			console.log("Plot data", data);
		} catch (jsonError) {
			if (!response.ok) {
				throw new Error(
					`Server returned status ${response.status}: ${response.statusText}. Response was not valid JSON.`
				);
			}
			throw new Error(
				`Successfully fetched, but response was not valid JSON. ${jsonError.message}`
			);
		}
		if (!response.ok) {
			const errorMessage =
				data?.error ||
				`Failed to fetch plot data. Server responded with status ${response.status}.`;
			throw new Error(errorMessage);
		}

		updatePlot(data)
	} catch (error) {
		console.error("Error in getPlotData:", error);
	} finally {
		document.getElementById("plot-loading-spinner").style.display = "none";
	}
}

document.addEventListener("DOMContentLoaded", function () {
	getPlotData(
		"umap_plot",
		`/analysis/umap_plot/${window.analysis.id}`,
		"GET"
	)
		.then(() => {
			console.log("Plot loaded successfully.");
		})
		.catch((err) => {
			console.error("Failed to load plot:", err);
			alert("Failed to load plot. Please try again later.");
		});

	window.addEventListener("resize", () => {
		Plotly.Plots.resize("scatterPlot");
	});
});
