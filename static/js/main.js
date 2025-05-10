document.addEventListener("DOMContentLoaded", function () {
	const flashMessages = document.querySelectorAll("#flash-message-container .flash-message-item");
	const FADE_DELAY = 3000; // 5 seconds
	const REMOVE_DELAY_AFTER_FADE = 500;

	flashMessages.forEach(function (message, index) {
		// Stagger the appearance slightly if desired, or just use a single timeout
		setTimeout(function () {
			// Add fade-out class for CSS transition
			message.classList.add("fade-out");

			// Remove the element from DOM after the fade-out transition completes
			setTimeout(function () {
				if (message.parentNode) {
					// Check if it hasn't been removed by other means
					message.remove();
				}
			}, REMOVE_DELAY_AFTER_FADE);
		}, FADE_DELAY + index * 300); // Stagger multiple messages slightly
	});

	const colorBySelect = document.getElementById("color-by");
	const geneSelectionDiv = document.getElementById("gene-selection-div");
	const pointSizeSlider = document.getElementById("point-size");
	const pointSizeValue = document.getElementById("point-size-value");
	const opacitySlider = document.getElementById("opacity");
	const opacityValue = document.getElementById("opacity-value");

	// Check if these elements exist before adding event listeners
	// This makes the script safe to include on pages where these elements are not present
	if (colorBySelect) {
		colorBySelect.addEventListener("change", function () {
			if (geneSelectionDiv) {
				if (this.value === "gene_expression") {
					geneSelectionDiv.classList.remove("hidden");
				} else {
					geneSelectionDiv.classList.add("hidden");
				}
			}
		});
	}

	if (pointSizeSlider && pointSizeValue) {
		pointSizeSlider.addEventListener("input", function () {
			pointSizeValue.textContent = this.value;
		});
	}

	if (opacitySlider && opacityValue) {
		opacitySlider.addEventListener("input", function () {
			opacityValue.textContent = this.value;
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

	getPlotData()

	// Update plot when window is resized
    window.addEventListener('resize', function () {
        Plotly.relayout('scatterPlot', {
            width: window.innerWidth * 0.44,
            height: window.innerHeight * 0.78
        });
    });

});

async function getPlotData() {
    const apiUrl = `/analysis/umap_plot/${analysis.id}`;
	try {
		const response = await fetch(apiUrl, { // Your Python endpoint
			method: 'GET',
			headers: {'Content-Type': 'application/json'},
		});

		let data;
		try {
			data = await response.json();
		} catch (jsonError) {
			if (!response.ok) {
				throw new Error(`Server returned status ${response.status}: ${response.statusText}. Response was not valid JSON.`);
			}
			throw new Error(`Successfully fetched, but response was not valid JSON. ${jsonError.message}`);
		}
		if (!response.ok) {
			const errorMessage = data?.error || `Failed to fetch plot data. Server responded with status ${response.status}.`;
			throw new Error(errorMessage);
		}
		console.log("Received data:", data);
		if (data.error) {
			throw new Error(data.error);
		}
		if (data.data && data.layout) {
			Plotly.newPlot("scatterPlot", data.data, data.layout);
		} else {
			throw new Error("Received data is not in the expected format (missing 'data' or 'layout' properties).");
		}
	} catch (error) {
		console.error("Error in getPlotData:", error);
	}
}
