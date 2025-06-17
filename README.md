# TF-PRED Webserver Setup Guide

## Prerequisites

- Python 3.12 or higher (can be installed via UV)
- [UV Package Manager](https://astral.sh/uv) (for Python dependencies)
- Node.js v22.11.x (for frontend dependencies)
- Git

## Installation Steps

### 1. Install UV Package Manager

```bash
# For macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# For Windows (PowerShell)
(Invoke-WebRequest -Uri "https://astral.sh/uv/install.ps1" -UseBasicParsing).Content | pwsh -Command -
```

### 2. Install Python 3.12 via UV

```bash
# Install Python 3.12 using UV
uv python install 3.12
```

### 3. Install Node.js v22.11.

```bash
# Install Node.js v22.11.x
curl -fsSL https://deb.nodesource.com/setup_22.x | sudo -E bash -
sudo apt-get install -y nodejs
```

### 4. Clone and Setup Project

```bash
# Clone the repository
git clone <repository-url>
cd tf-pred-webserver

# Create and activate virtual environment using UV
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

### 5. Setup Frontend Dependencies

```bash
# Install Node.js dependencies
npm install

npm run build:css  # Build CSS
```

## Running the Application

### Development Mode

```bash
# Terminal 1: Run frontend development server
npm run watch:css

# Terminal 2: Run Flask backend
uv run main.py
```



## Input File Requirements

This document outlines the required and optional input files for using the tool. Please review the formats carefully to ensure successful analysis.

### Quick Reference

| Input Type                     | Format(s)         | Required/Optional         | Notes                                                   |
| :----------------------------- | :---------------- | :------------------------ | :------------------------------------------------------ |
| **Gene Expression & Metadata** | `.h5ad` (AnnData) | Required (Method 1)       | The recommended, all-in-one format.                     |
| **Gene Expression Data**       | `.tsv`, `.csv`    | Required (Method 2)       | Cells as rows, genes as columns.                        |
| **Cell Metadata**              | `.tsv`, `.csv`    | Optional (Method 2)       | Cells as rows, metadata as columns. Highly recommended. |
| **2D Layout**                  | `.tsv`, `.csv`    | Optional                  | Use if you have a pre-computed layout (e.g., UMAP).     |
| **Species Setting**            | Tool setting      | Required if data is mouse | Select "Mouse" to map genes to human orthologs.         |

---

## Data Input Methods

You can provide your data using one of the two methods below.

### Method 1: AnnData Object (Recommended)

The primary and most convenient input format is a Scanpy `anndata` object, which has a `.h5ad` file extension.

- **File:** `data.h5ad`

This single file is expected to contain the gene expression matrix, cell metadata, and gene metadata. The tool will automatically parse this object to extract the necessary information.

<img src="https://falexwolf.de/img/scanpy/anndata.svg" alt="AnnData Structure" width="400" height="400">
*Figure: The structure of an AnnData object.*

### Method 2: Separate TSV/CSV Files

If you don't have an `anndata` object, you can provide your data as separate plain text files.

#### 1. Gene Expression Data (Required)

A TSV or CSV file containing the expression matrix.

- **Example:** `data.tsv`
- **Format:** The matrix must be structured with **cells as rows** and **genes as columns**.

#### 2. Cell Metadata (Optional)

A TSV or CSV file containing metadata for each cell. This file is optional but highly recommended for richer visualization and analysis.

- **Example:** `metadata.tsv`
- **Format:** The file must have **cells as rows** and metadata attributes as columns. The order of cells should match the gene expression file.

---

## Optional Inputs & Settings

### Custom 2D Layout

If you have already computed a 2D layout (e.g., from UMAP or t-SNE) and want to use it for visualization instead of running a new analysis, you can provide a layout file.

- **How to use:** Select the "I have a 2D layout file" option in the tool.
- **File:** A TSV or CSV file (e.g., `layout.tsv`).
- **Format:** The file must contain cells as rows. The header must follow this exact format:
  ```
  ,X_umap1,X_umap2,X_pca1,X_pca2,Cluster
  ```
  > **Important:** Note the required **leading comma** before `X_umap1`.

### Species Selection

If your gene expression data is from a mouse model, you must select the appropriate species in the tool's interface.

- **How to use:** Select "Mouse" from the species dropdown menu.
- **Functionality:** When this option is selected, the tool will automatically map the mouse gene names in your data to their human orthologs using a built-in `mouse2human` mapping file. If your data is human, you can leave the default setting.

## Troubleshooting

1. If you encounter permission issues with UV installation, try:

   ```bash
   sudo curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. For Node.js dependency issues:

   ```bash
   rm -rf node_modules
   npm install
   ```

3. If Flask server fails to start:
   - Ensure you're in the virtual environment
   - Check if port 5000 is available
   - Verify all Python dependencies are installed correctly
