// Use window.analysis for the analysis object
console.log("Analysis object:", window.analysis);
// Add view_analysis-specific JS here, e.g., Plotly rendering, form listeners, etc.

const defaultPointSize = 4;
const defaultOpacity = 0.5;

// Plot configuration panel logic
const plotTypeSelect = document.getElementById("plot-type");
const colorBySelect = document.getElementById("color-by");
const metadataColSelectionDiv = document.getElementById("metadata-column-selection-div");
const tfSelectionDiv = document.getElementById("tf-selection-div");
const pointSizeSlider = document.getElementById("point-size");
const pointSizeValue = document.getElementById("point-size-value");
const opacitySlider = document.getElementById("opacity");
const opacityValue = document.getElementById("opacity-value");
const showLegendCheckbox = document.getElementById("show-legend");

function resetPlotConfig() {
	if (pointSizeSlider) {
		pointSizeSlider.value = defaultPointSize;
		pointSizeValue.textContent = defaultPointSize;
	}
	if (opacitySlider) {
		opacitySlider.value = defaultOpacity;
		opacityValue.textContent = defaultOpacity;
	}
	if (showLegendCheckbox) {
		showLegendCheckbox.checked = false;
	}
}

plotTypeSelect.addEventListener("change", function () {
	if (this.value === "umap_plot") {  // UMAP Plot
		console.log("UMAP plot selected");
		getPlotData(this.value)
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
		getPlotData(this.value)
		.then(() => {
			console.log("PCA Plot loaded successfully.");
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
			metadataColSelectionDiv.classList.add("hidden");
		} else if (this.value === "metadata_columns") {
			tfSelectionDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.remove("hidden");
		}
		else {
			tfSelectionDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.add("hidden");
		}
	});
}

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
		alert("Plot update functionality not yet implemented.");
	});
}

// UMAP plot loading and rendering
async function getPlotData(plot_type = "umap_plot") {
	document.getElementById("plot-loading-spinner").style.display = "block";

	const apiUrl = plot_type === "umap_plot"
	    ? `/analysis/umap_plot/${window.analysis.id}`
	    : `/analysis/pca_plot/${window.analysis.id}`;

	try {
		const response = await fetch(apiUrl, {
			method: "GET",
			headers: { "Content-Type": "application/json" },
		});

		let data;
		try {
			data = await response.json();
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

		if (data.data && data.layout) {
			data.layout.dragmode = "pan"; // Set default to pan
			data.layout.legend = {
				x: -5,
				xanchor: "left",
				y: 10,
				yanchor: "top",
			};
			Plotly.newPlot("scatterPlot", data.data, data.layout, {
				responsive: true,
				displayModeBar: true,
				displaylogo: false,
				scrollZoom: true,
			});

			// Reset plot configuration
			resetPlotConfig();
		} else {
			throw new Error(
				"Received data is not in the expected format (missing 'data' or 'layout' properties)."
			);
		}
	} catch (error) {
		console.error("Error in getPlotData:", error);
	} finally {
		document.getElementById("plot-loading-spinner").style.display = "none";
	}
}

document.addEventListener("DOMContentLoaded", function () {
	getPlotData()
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
