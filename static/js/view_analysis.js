console.log("Analysis object:", window.analysis);

const defaultPointSize = 4;
const defaultOpacity = 0.5;

// Plot configuration panel logic
const plotConfigForm = document.getElementById("plot-config-form");
const plotTitle = document.getElementById("plot-title");
const plotTypeSelect = document.getElementById("plot-type");
const colorBySelect = document.getElementById("color-by");
const metadataColSelectionDiv = document.getElementById("metadata-column-selection-div");
const metadataColNameSelect = document.getElementById("metadata-column-name");
const tfSelectionDiv = document.getElementById("tf-selection-div");
const tfManualEntryDiv = document.getElementById("tf-manual-entry-div");
const geneEntryDiv = document.getElementById("gene-entry-div");
const colorScaleSelectionDiv = document.getElementById("color-scale-selection-div");
const tfNameSelect = document.getElementById("tf-name-select");
const tfManualEntryInput = document.getElementById("tf-manual-entry-input");
const geneEntryInput = document.getElementById("gene-entry-input");
const pointSizeSlider = document.getElementById("point-size");
const pointSizeValue = document.getElementById("point-size-value");
const opacitySlider = document.getElementById("opacity");
const opacityValue = document.getElementById("opacity-value");
const showLegendCheckbox = document.getElementById("show-legend");
const fdrLevelDiv = document.getElementById("fdr-level-div");
const fdrLevel = document.getElementById("fdr-level");
const pValueThresholdDiv = document.getElementById("p-val-threshold-div");
const pValueThreshold = document.getElementById("p-val-threshold");
const reRunFDRCorrectionButton = document.getElementById("re-run-fdr-correction");
const changePValueThreshold = document.getElementById("change-p-value-threshold");
const moreInfoBtn = document.getElementById("more-info-btn");
const modal = document.getElementById("more-info-modal");
const closeModalBtn = document.getElementById("close-modal-btn");
const totalCells = document.getElementById("total-cells");
const plotLoadingSpinner = document.getElementById("plot-loading-spinner")
const changeThresholdTypeDiv = document.getElementById("change-threshold-type-div")
const fdrCorrectionRadio = document.getElementById("fdr-correction-radio")
const pValueThresholdRadio = document.getElementById("p-value-threshold-radio")


function clearAllSelections() {
	if (colorBySelect)
		colorBySelect.value = "select_cluster_type";
	if (metadataColNameSelect){
		metadataColNameSelect.value = "select_metadata_column";
		metadataColSelectionDiv.classList.add("hidden");
	}
	if (tfNameSelect) {
		tfNameSelect.value = "select_tf";
		tfSelectionDiv.classList.add("hidden");
		changeThresholdTypeDiv.classList.add("hidden");
	}
	if (tfManualEntryInput) {
		tfManualEntryInput.value = "";
		tfManualEntryDiv.classList.add("hidden");
	}
	if (geneEntryInput) {
		geneEntryInput.value = "";
		geneEntryDiv.classList.add("hidden");
		colorScaleSelectionDiv.classList.add("hidden");
	}
	if (pointSizeSlider)
		pointSizeSlider.value = defaultPointSize;
	if (opacitySlider)
		opacitySlider.value = defaultOpacity;
	if (fdrLevel)
		fdrLevel.value = "";
	if (pValueThreshold)
		pValueThreshold.value = "";
}

plotTypeSelect.addEventListener("change", function () {

	clearAllSelections();

	const apiUrl = `/analysis/plot/${window.analysis.id}`;
	getPlotData(apiUrl, "POST", {"plot_type": this.value})
	.then(() => {
		console.log("Plot loaded successfully.");
	})
	.catch((err) => {
		alert("Failed to load plot. Please try again later.");
	});
});

if (colorBySelect) {
	colorBySelect.addEventListener("change", function () {
		if (this.value === "tf_activity") {
			tfSelectionDiv.classList.remove("hidden");
			tfManualEntryDiv.classList.remove("hidden");
			metadataColSelectionDiv.classList.add("hidden");
			geneEntryDiv.classList.add("hidden");
			colorScaleSelectionDiv.classList.add("hidden");
			fdrLevelDiv.classList.remove("hidden");
			changeThresholdTypeDiv.classList.remove("hidden");

			// Clear selections
			metadataColNameSelect.value = "select_metadata_column";
			geneEntryInput.value = "";
		}
		else if (this.value === "metadata_columns") {
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.remove("hidden");
			geneEntryDiv.classList.add("hidden");
			colorScaleSelectionDiv.classList.add("hidden");
			fdrLevelDiv.classList.add("hidden");
			changeThresholdTypeDiv.classList.add("hidden");

			// Clear selections
			tfNameSelect.value = "select_tf";
			tfManualEntryInput.value = "";
			fdrLevel.value = "";
			pValueThreshold.value = "";
			geneEntryInput.value = "";

		}
		else if (this.value === "gene_expression") {
			geneEntryDiv.classList.remove("hidden");
			colorScaleSelectionDiv.classList.remove("hidden");
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.add("hidden");
			fdrLevelDiv.classList.add("hidden");
			changeThresholdTypeDiv.classList.add("hidden");

			// Clear selections
			tfNameSelect.value = "select_tf";
			tfManualEntryInput.value = "";
			fdrLevel.value = "";
			pValueThreshold.value = "";
			metadataColNameSelect.value = "select_metadata_column";
		}
		else {
			geneEntryDiv.classList.add("hidden");
			colorScaleSelectionDiv.classList.add("hidden");
			tfSelectionDiv.classList.add("hidden");
			tfManualEntryDiv.classList.add("hidden");
			metadataColSelectionDiv.classList.add("hidden");
			fdrLevelDiv.classList.add("hidden");
			changeThresholdTypeDiv.classList.add("hidden");
		}
	});
}

function updatePlot(plot_data){
	if (plot_data.data && plot_data.layout) {
		total_cells = plot_data.data.reduce((acc, trace) => acc + (trace.x ? trace.x.length : 0), 0);
		totalCells.textContent = `Total Cells: ${total_cells}`;
		console.log(plot_data);

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
			if (trace.cluster && trace.x && trace.x.length !== undefined) {
				const clusterCount = trace.x.length;
				trace.cluster = `${trace.cluster} (${clusterCount})`;
				trace.name = `${trace.name} (${clusterCount})`;
			}
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

		if (plot_data.fdr_level)
			fdrLevel.value = plot_data.fdr_level;

		if (plot_data.p_value_threshold)
			pValueThreshold.value = plot_data.p_value_threshold;

		if (plot_data.layout.title)
			plotTitle.textContent = plot_data.layout.title;

		if(plot_data.p_value_threshold)
			pValueThreshold.value = plot_data.p_value_threshold;

	}
	else
		throw new Error("Received data is not in the expected format.");
}

metadataColNameSelect.addEventListener("change", function () {
	if (this.value !== "select_metadata_column") {
		const apiUrl = `/analysis/metadata-cluster/${window.analysis.id}`;

		updatePlotData(
			apiUrl,
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

		updatePlotData(
			apiUrl,
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

			updatePlotData(
				apiUrl,
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

			updatePlotData(
				apiUrl,
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

colorScaleSelectionDiv.addEventListener("change", function (e) {
	console.log("ColorScaleSelectionDiv changed successfully.");
	Plotly.restyle('scatterPlot', {
		'marker.colorscale': [e.target.value]
    }, [0]);
});

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

reRunFDRCorrectionButton.addEventListener("click", async function (e) {
	try {
		plotLoadingSpinner.style.display = "block";
		const response = await fetch(
			`/analysis/re-run-fdr-correction/${window.analysis.id}`,
			{
				method: "POST",
				headers: {"Content-Type": "application/json"},
				body: JSON.stringify({fdr_level: fdrLevel.value})
			}
		);
		if (!response.ok) {
			alert(
				`Failed to re-run FDR correction. Server responded with status ${response.status}.`
			)
		}else{
			alert("FDR correction re-run successfully. Reload the page to see the updated plot.");
			window.location.reload();
		}
		plotLoadingSpinner.style.display = "none";
	}
	catch (error) {
		console.error("Error re-running FDR correction:", error);
		alert("Failed to re-run FDR correction. Please try again later.");
	}
});

changePValueThreshold.addEventListener("click", async function (e) {
	try {
		plotLoadingSpinner.style.display = "block";
		const response = await fetch(
			`/analysis/change-p-value-threshold/${window.analysis.id}`,
			{
				method: "POST",
				headers: {"Content-Type": "application/json"},
				body: JSON.stringify({p_value_threshold: pValueThreshold.value})
			}
		);
		if (!response.ok) {
			alert(
				`Failed to change p-value threshold. Server responded with status ${response.status}.`
			)
		}else{
			alert("P-value threshold changed successfully. Reload the page to see the updated plot.");
			window.location.reload();
		}
		plotLoadingSpinner.style.display = "none";
	}
	catch (error) {
		console.error("Error changing p-value threshold:", error);
		alert("Failed to change p-value threshold. Please try again later.");
	}
});


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

async function updatePlotData(apiUrl, method, body) {
	plotLoadingSpinner.style.display = "block";

	try {
		// 1. Make the API request
		const response = await fetch(apiUrl, {
			method: method,
			headers: { "Content-Type": "application/json" },
			body: body ? JSON.stringify(body) : null,
		});

		// 2. Check for HTTP errors first. If the response is not ok, process the error and throw.
		if (!response.ok) {
			let errorMessage = `HTTP error ${response.status}: ${response.statusText}`;
			try {
				const errorData = await response.json();
				errorMessage = errorData.error || errorMessage;
			} catch (e) {
				alert("Failed to parse error response as JSON.");
				console.warn("Failed to parse error response as JSON:", e);
			}
			throw new Error(errorMessage);
		}

		// Update coordinates and layout
		const data = await response.json();
		updatePlot(data);
	} catch (error) {
		console.error("Error in updatePlotData:", error);
		alert(`Error loading plot: ${error.message}`);
	} finally {
		plotLoadingSpinner.style.display = "none";
	}
}

async function getPlotData(apiUrl, method, body) {
	plotLoadingSpinner.style.display = "block";

	try {
		// 1. Make the API request
		const response = await fetch(apiUrl, {
			method: method,
			headers: { "Content-Type": "application/json" },
			body: body ? JSON.stringify(body) : undefined,
		});

		// 2. Check for HTTP errors first. If the response is not ok, process the error and throw.
		if (!response.ok) {
			let errorMessage = `HTTP error ${response.status}: ${response.statusText}`;
			try {
				const errorData = await response.json();
				errorMessage = errorData.error || errorMessage;
			} catch (e) {
				alert("Failed to parse error response as JSON.");
				console.warn("Failed to parse error response as JSON:", e);
			}
			throw new Error(errorMessage);
		}

		const data = await response.json();

		updatePlot(data);
	} catch (error) {
		console.error("Error loading plot data:", error);
		alert(`Error loading plot: ${error.message}`);
	} finally {
		plotLoadingSpinner.style.display = "none";
	}
}

document.addEventListener("DOMContentLoaded", function () {

	addInputValidation("fdr-level",0.0, 1.0);
	addInputValidation("p-val-threshold",0.0, 1.0);

	getPlotData(
		`/analysis/plot/${window.analysis.id}`,
		"POST",
		{"plot_type": "umap_plot"}
	)
	.then(() => {
		console.log("Plot loaded successfully.");
	})
	.catch((err) => {
		alert("Failed to load plot. Please try again later.");
	});

	window.addEventListener("resize", () => {
		Plotly.Plots.resize("scatterPlot");
	});

    if (moreInfoBtn && modal && closeModalBtn) {
        moreInfoBtn.addEventListener("click", () => {
            modal.classList.remove("hidden");
        });
        closeModalBtn.addEventListener("click", () => {
            modal.classList.add("hidden");
        });
        modal.addEventListener("click", (e) => {
            if (e.target === modal) {
                modal.classList.add("hidden");
            }
        });
        document.addEventListener("keydown", (e) => {
            if (e.key === "Escape") {
                modal.classList.add("hidden");
            }
        });
    }

	if (fdrCorrectionRadio && pValueThresholdRadio) {
		fdrCorrectionRadio.addEventListener("change", function () {
			if (this.checked) {
				fdrLevelDiv.classList.remove("hidden");
				pValueThresholdDiv.classList.add("hidden");
			}
		});

		pValueThresholdRadio.addEventListener("change", function () {
			if (this.checked) {
				fdrLevelDiv.classList.add("hidden");
				pValueThresholdDiv.classList.remove("hidden");
			}
		});
	}

});
