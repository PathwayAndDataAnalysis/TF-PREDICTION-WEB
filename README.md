# Input File Requirements

This document outlines the required and optional input files for using the tool. Please review the formats carefully to ensure successful analysis.

### Quick Reference

| Input Type | Format(s) | Required/Optional | Notes |
| :--- | :--- | :--- | :--- |
| **Gene Expression & Metadata** | `.h5ad` (AnnData) | Required (Method 1) | The recommended, all-in-one format. |
| **Gene Expression Data** | `.tsv`, `.csv` | Required (Method 2) | Cells as rows, genes as columns. |
| **Cell Metadata** | `.tsv`, `.csv` | Optional (Method 2) | Cells as rows, metadata as columns. Highly recommended. |
| **2D Layout** | `.tsv`, `.csv` | Optional | Use if you have a pre-computed layout (e.g., UMAP). |
| **Species Setting** | Tool setting | Required if data is mouse | Select "Mouse" to map genes to human orthologs. |

***

## Data Input Methods

You can provide your data using one of the two methods below.

### Method 1: AnnData Object (Recommended)

The primary and most convenient input format is a Scanpy `anndata` object, which has a `.h5ad` file extension.

*   **File:** `data.h5ad`

This single file is expected to contain the gene expression matrix, cell metadata, and gene metadata. The tool will automatically parse this object to extract the necessary information.

<img src="https://falexwolf.de/img/scanpy/anndata.svg" alt="AnnData Structure" width="400" height="400">
*Figure: The structure of an AnnData object.*

### Method 2: Separate TSV/CSV Files

If you don't have an `anndata` object, you can provide your data as separate plain text files.

#### 1. Gene Expression Data (Required)
A TSV or CSV file containing the expression matrix.

*   **Example:** `data.tsv`
*   **Format:** The matrix must be structured with **cells as rows** and **genes as columns**.

#### 2. Cell Metadata (Optional)
A TSV or CSV file containing metadata for each cell. This file is optional but highly recommended for richer visualization and analysis.

*   **Example:** `metadata.tsv`
*   **Format:** The file must have **cells as rows** and metadata attributes as columns. The order of cells should match the gene expression file.

***

## Optional Inputs & Settings

### Custom 2D Layout

If you have already computed a 2D layout (e.g., from UMAP or t-SNE) and want to use it for visualization instead of running a new analysis, you can provide a layout file.

*   **How to use:** Select the "I have a 2D layout file" option in the tool.
*   **File:** A TSV or CSV file (e.g., `layout.tsv`).
*   **Format:** The file must contain cells as rows. The header must follow this exact format:
    ```
    ,X_umap1,X_umap2,X_pca1,X_pca2,Cluster
    ```
    > **Important:** Note the required **leading comma** before `X_umap1`.

### Species Selection

If your gene expression data is from a mouse model, you must select the appropriate species in the tool's interface.

*   **How to use:** Select "Mouse" from the species dropdown menu.
*   **Functionality:** When this option is selected, the tool will automatically map the mouse gene names in your data to their human orthologs using a built-in `mouse2human` mapping file. If your data is human, you can leave the default setting.