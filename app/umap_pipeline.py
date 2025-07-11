import os
import pandas as pd
import scanpy as sc
import psutil
from flask import current_app


def validate_input_files(analysis_data):
    """Validate that all required input files exist and are readable."""
    gene_expr = analysis_data["inputs"]["gene_expression"]

    if gene_expr["source"] == "h5ad":
        filepath = gene_expr["h5ad_filepath"]
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"[UMAP] H5AD file not found: {filepath}")

        # Check file size
        file_size = os.path.getsize(filepath)
        if file_size == 0:
            raise ValueError(f"[UMAP] H5AD file is empty: {filepath}")

        # Try to read the file
        try:
            adata = sc.read_h5ad(filepath)
            if adata.n_obs == 0 or adata.n_vars == 0:
                raise ValueError(f"[UMAP] H5AD file contains no data: {filepath}")
        except Exception as e:
            raise ValueError(f"[UMAP] Failed to read H5AD file {filepath}: {e}")

    else:
        # Validate separate files
        gene_exp_filepath = gene_expr.get("gene_exp_filepath")
        metadata_filepath = gene_expr.get("metadata_filepath")

        if not gene_exp_filepath or not os.path.exists(gene_exp_filepath):
            raise FileNotFoundError(f"[UMAP] Gene expression file not found: {gene_exp_filepath}")

        if metadata_filepath and not os.path.exists(metadata_filepath):
            raise FileNotFoundError(f"[UMAP] Metadata file not found: {metadata_filepath}")


def run_umap_pipeline(
    user_id, analysis_id, analysis_data, have_2d_layout, update_status_fn, run_analysis_fn
):
    try:
        current_app.logger.info(f"[UMAP] Starting UMAP pipeline for user '{user_id}', analysis '{analysis_id}'.")

        # Validate input files first
        validate_input_files(analysis_data)

        # Check system resources
        memory = psutil.virtual_memory()
        if memory.percent > 85:
            raise MemoryError("[UMAP] Insufficient memory for UMAP processing")

        # 1. Load data
        gene_expr = analysis_data["inputs"]["gene_expression"]
        metadata_cols = []

        if gene_expr["source"] == "h5ad":
            current_app.logger.info(f"[UMAP] Loading .h5ad file: {gene_expr['h5ad_filepath']}")
            adata = sc.read_h5ad(gene_expr["h5ad_filepath"])

            # Save metadata obs_keys in analysis
            metadata_cols = adata.obs_keys() if adata.obs_keys() else []
            current_app.logger.info(f"[UMAP] Metadata columns found in .h5ad file: {metadata_cols}")

        else:
            current_app.logger.info(f"[UMAP] Loading gene_exp file: {gene_expr['gene_exp_filepath']}")
            gene_exp = pd.read_csv(gene_expr['gene_exp_filepath'], index_col=0)

            # A raw count matrix should not have NaNs. If it does, it's a missing value. Fill with 0.
            if gene_exp.isnull().values.any():
                current_app.logger.warning(
                    f"Found {gene_exp.isnull().values.sum()} missing values in gene expression matrix. Filling with 0.")
                gene_exp.fillna(0, inplace=True)

            if gene_expr.get("metadata_filepath"):
                current_app.logger.info(f"[UMAP] Loading metadata file: {gene_expr['metadata_filepath']}")
                meta_data = pd.read_csv(gene_expr["metadata_filepath"], index_col=0)

                # Ensure metadata and expression data have the same cells
                common_cells = gene_exp.index.intersection(meta_data.index)
                current_app.logger.info(f"Found {len(common_cells)} common cells between expression and metadata.")

                if len(common_cells) > 0:
                    gene_exp = gene_exp.loc[common_cells]
                    meta_data = meta_data.loc[common_cells]

                    adata = sc.AnnData(gene_exp, obs=meta_data)
                    metadata_cols = meta_data.columns.tolist()
                else:
                    current_app.logger.warning(
                        "[UMAP] No common cells found between gene expression and metadata files. "
                        "Creating AnnData with gene expression only."
                    )
                    adata = sc.AnnData(gene_exp)

            else:
                current_app.logger.info("[UMAP] No metadata file provided, creating AnnData with gene_exp only.")
                adata = sc.AnnData(gene_exp)


        # ------- Check for Ensembl IDs and map to gene symbols -------
        mapped = False
        if adata.var.index.str.startswith("ENSG").mean() > 0.5:
            current_app.logger.info("[UMAP] Detected Ensembl IDs in var index. Converting to gene symbols.")
            common_names = ['gene_symbols', 'symbol', 'gene_name', 'symbols']

            for g_name in common_names:
                if g_name in adata.var.columns:
                    current_app.logger.info(
                        f"[UMAP] Found gene symbols in column '{g_name}'. Converting var index.")
                    adata.var_names = adata.var[g_name]
                    adata.var_names_make_unique()
                    current_app.logger.info(f"[UMAP] Converted var index to gene symbols using column '{g_name}'.")
                    mapped = True
                    break

            if not mapped and "feature_name" in adata.var.columns:
                current_app.logger.warning("[UMAP] No gene symbols found in var columns. Using 'feature_name' as fallback.")
                try:
                    adata.var["gene_symbols"] = adata.var["feature_name"].str.split("_").str[0]
                    adata.var_names = adata.var["gene_symbols"]
                    adata.var_names_make_unique()
                    mapped = True
                except Exception as e:
                    current_app.logger.error(
                        f"[UMAP] Failed to extract gene symbols from 'feature_name': {e}. Using original var index."
                    )

            if not mapped:
                current_app.logger.info("[UMAP] No luck mapping Ensembl IDs to gene symbols. Using auto-detection.")
                if gene_expr.get("species") == "human":
                    ensembl_map_file = "human_gencode_mapping.csv"
                elif gene_expr.get("species") == "mouse":
                    ensembl_map_file = "mouse_gencode_mapping.csv"
                # TODO: Add support for other species if needed
                else:
                    ensembl_map_file = "human_gencode_mapping.csv"

                if ensembl_map_file:
                    script_dir = os.path.dirname(os.path.abspath(__file__))
                    gene_ensembl_map_df = pd.read_csv(os.path.join(script_dir, "..", "prior_data", ensembl_map_file))

                    # gene_ensembl_map_df is dataframe with columns: "ensembl_id", "gene_symbol". Now we can map the Ensembl IDs to gene symbols
                    # TODO: check if adata.var.index has ensembl IDs with dots (e.g. ENSG00000047056.18)
                    adata.var.index = adata.var.index.str.split(".").str[0]
                    adata.var["gene_symbols"] = adata.var.index.map(
                        gene_ensembl_map_df.set_index("ensembl_id")["gene_symbol"]
                    )
                    adata.var_names = adata.var["gene_symbols"]
                    # adata.var_names_make_unique()
                    adata.var.index = adata.var.index.astype(str)  # Convert to string
                    adata.var_names_make_unique()

                    current_app.logger.info("[UMAP] Successfully mapped Ensembl IDs to gene symbols.")


        # ------- Data Filtering and Preprocessing -------
        current_app.logger.info("[UMAP] Starting data filtering and preprocessing...")
        data_filtering = analysis_data["inputs"]["data_filtering"]
        # 1. Filter Cells
        if data_filtering.get("filter_cells"):
            current_app.logger.info(f"[UMAP] Filtering cells with min_genes={data_filtering.get('min_genes', 0)}")
            sc.pp.filter_cells(adata, min_genes=data_filtering.get("min_genes", 0))
            current_app.logger.info(f"Shape after filtering cells: {adata.n_obs} cells × {adata.n_vars} genes")
            if adata.n_obs == 0:
                update_status_fn(user_id=user_id, analysis_id=analysis_id, status="Error", error="No cells left after filtering. Please check your filter settings.")
                raise ValueError("No cells left after filtering. Please check your filter settings.")

        # 2. Filter genes
        if data_filtering.get("filter_genes"):
            current_app.logger.info(f"[UMAP] Filtering genes with min_cells={data_filtering.get('min_cells', 0)}")
            sc.pp.filter_genes(adata, min_cells=data_filtering.get("min_cells", 0))
            current_app.logger.info(f"Shape after filtering genes: {adata.n_obs} cells × {adata.n_vars} genes")
            if adata.n_vars == 0:
                update_status_fn(user_id=user_id, analysis_id=analysis_id, status="Error", error="No genes left after filtering. Please check your filter settings.")
                raise ValueError("No genes left after filtering. Please check your filter settings.")

        # 3. Filter on mitochondrial gene percentage (with species handling)
        if data_filtering.get("qc_filter"):
            current_app.logger.info(f"[UMAP] Filtering genes with min_counts={data_filtering.get('max_mt_pct', 100)}% mitochondrial counts")
            mito_prefix = None
            species = analysis_data["inputs"]["gene_expression"]["species"].lower()

            if species == 'auto':
                # Auto-detect based on gene name patterns
                has_human_mt = any(adata.var_names.str.startswith('MT-'))
                has_mouse_mt = any(adata.var_names.str.startswith('mt-'))
                if has_human_mt:
                    mito_prefix = 'MT-'
                    current_app.logger.info("   -> Auto-detected human mitochondrial genes (prefix 'MT-').")
                elif has_mouse_mt:
                    mito_prefix = 'mt-'
                    current_app.logger.info("   -> Auto-detected mouse mitochondrial genes (prefix 'mt-').")
            elif species == 'human':
                mito_prefix = 'MT-'
            elif species == 'mouse':
                mito_prefix = 'mt-'
            else:
                raise ValueError(f"Invalid species: '{species}'. Choose from 'human', 'mouse', or 'auto'.")

            if mito_prefix:
                adata.var['mt'] = adata.var_names.str.startswith(mito_prefix)

                # --- CRITICAL CHECK ---
                num_mt_genes_found = adata.var['mt'].sum()
                current_app.logger.info(f"  -> Found {num_mt_genes_found} mitochondrial genes with prefix '{mito_prefix}'.")

                if num_mt_genes_found > 0:
                    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
                    mito_threshold = float(data_filtering.get("max_mt_pct", 100)) # Default to 100% (no filtering) if not provided
                    adata = adata[adata.obs.pct_counts_mt < mito_threshold, :].copy()

                    # catastrophic case where all cells are removed
                    if adata.n_obs == 0:
                        current_app.logger.warning("  -> WARNING: All cells were removed by the mitochondrial filter. "
                                                   "Consider increasing the threshold or examining the data quality.")
                        update_status_fn(user_id=user_id, analysis_id=analysis_id, status="Error", error="No cells left after mitochondrial filtering.")
                        raise ValueError("No cells left after mitochondrial filtering. Please check your filter settings.")
            else:
                current_app.logger.info("   -> WARNING: Could not find mitochondrial genes. Skipping mitochondrial filtering step.")

        # 4. Normalize total expression per cell
        if data_filtering.get("data_normalize"):
            current_app.logger.info(f"[UMAP] Normalizing total with target_sum={data_filtering.get('data_normalize_value')}")
            sc.pp.normalize_total(adata, target_sum=data_filtering.get("data_normalize_value"))

        # 5. Apply log1p transformation
        if data_filtering.get("log_transform"):
            current_app.logger.info("[UMAP] Applying log1p transformation.")
            sc.pp.log1p(adata)
        current_app.logger.info(f"\nFinal processed data shape: {adata.shape}")


        # -------- Dimensionality Reduction and UMAP --------
        if not have_2d_layout:
            params = analysis_data["inputs"]["layout"]["umap_settings"]
            current_app.logger.info(f"[UMAP] Pipeline parameters: {params}")

            current_app.logger.info(f"[UMAP] Running PCA with n_comps={params.get('pca_components')}")
            sc.pp.pca(adata, n_comps=params.get("pca_components"))

            current_app.logger.info(
                f"[UMAP] Computing neighbors with n_neighbors={params.get('n_neighbors')}, metric={params.get('metric', 'euclidean')}"
            )
            sc.pp.neighbors(
                adata,
                n_neighbors=params.get("n_neighbors"),
                use_rep='X_pca',  # Use PCA results as input
                metric=params.get("metric", "euclidean")
            )

            current_app.logger.info(f"[UMAP] Running UMAP with min_dist={params.get('min_dist')}")
            if params.get("random_state") is not None and params.get("random_state") != 0:
                current_app.logger.info(f"[UMAP] Using random_state={params.get('random_state')} for reproducibility.")
                sc.tl.umap(
                    adata,
                    min_dist=params.get("min_dist"),
                    random_state=params.get("random_state"),
                )
            else:
                current_app.logger.info("[UMAP] No random_state provided, using default.")
                sc.tl.umap(adata, min_dist=params.get("min_dist"))

            # --- Save UMAP coordinates to results path
            result_path = analysis_data["results_path"]
            os.makedirs(result_path, exist_ok=True)
            umap_df = pd.DataFrame(
                adata.obsm["X_umap"], index=adata.obs_names, columns=["X_umap1", "X_umap2"]
            )
            # if "X_pca" in adata.obsm:
            #     umap_df["X_pca1"] = adata.obsm["X_pca"][:, 0]
            #     umap_df["X_pca2"] = adata.obsm["X_pca"][:, 1]
            if "X_pca" in adata.obsm:
                umap_df = pd.concat(
                    [umap_df,
                     pd.DataFrame(adata.obsm["X_pca"][:, :2], index=adata.obs_names, columns=["X_pca1", "X_pca2"])],
                    axis=1
                )

            umap_df["Cluster"] = 0  # or assign as needed

            umap_csv_path = os.path.join(result_path, "umap_coordinates.csv")
            current_app.logger.info(f"[UMAP] Saving UMAP coordinates to: {umap_csv_path}")
            umap_df.to_csv(umap_csv_path)

            # Update status to Completed
            update_status_fn(user_id=user_id, analysis_id=analysis_id, status="Completed", umap_csv_path=umap_csv_path, metadata_cols=metadata_cols)
            current_app.logger.info(f"[UMAP] UMAP pipeline completed successfully for analysis '{analysis_id}'.")

        # 5. Run analysis function
        current_app.logger.info(f"[UMAP] Running analysis function for analysis '{analysis_id}'.")

        try:
            run_analysis_fn(
                user_id=user_id,
                analysis_id=analysis_id,
                analysis_data=analysis_data,
                adata=adata,
                update_analysis_status_fn=update_status_fn,
            )
        except Exception as e:
            current_app.logger.error(f"[UMAP] Error in analysis function: {e}")
            update_status_fn(user_id, analysis_id, "Failed", error=str(e))
            raise e

    except FileNotFoundError as e:
        current_app.logger.error(f"[UMAP] File not found: {e}")
        update_status_fn(user_id, analysis_id, "Failed", error=f"Input file not found: {e}")
        raise e

    except MemoryError as e:
        current_app.logger.error(f"[UMAP] Memory error: {e}")
        update_status_fn(user_id, analysis_id, "Failed", error="Insufficient memory for processing")
        raise e

    except TimeoutError as e:
        current_app.logger.error(f"[UMAP] Timeout error: {e}")
        update_status_fn(user_id, analysis_id, "Failed", error="Processing timed out")
        raise e

    except Exception as e:
        current_app.logger.error(
            f"[UMAP] UMAP pipeline failed for analysis '{analysis_id}': {e}",
            exc_info=True,
        )
        update_status_fn(user_id, analysis_id, "Failed", error=str(e))
