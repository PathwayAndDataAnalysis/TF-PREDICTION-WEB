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
            metadata_cols = adata.obs_keys()[1:] if adata.obs_keys() else []
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
                gene_exp = gene_exp.loc[common_cells]
                meta_data = meta_data.loc[common_cells]

                adata = sc.AnnData(gene_exp, obs=meta_data)
                metadata_cols = meta_data.columns.tolist()
            else:
                current_app.logger.info("[UMAP] No metadata file provided, creating AnnData with gene_exp only.")
                adata = sc.AnnData(gene_exp)


        # ------- Data Filtering and Preprocessing -------
        current_app.logger.info("[UMAP] Starting data filtering and preprocessing...")
        data_filtering = analysis_data["inputs"]["data_filtering"]
        # 1. Filter Cells
        if data_filtering.get("filter_cells"):
            current_app.logger.info(f"[UMAP] Filtering cells with min_genes={data_filtering.get('filter_cells_value')}")
            sc.pp.filter_cells(adata, min_genes=data_filtering.get("filter_cells_value"))
            current_app.logger.info(f"Shape after filtering cells: {adata.n_obs} cells × {adata.n_vars} genes")
            if adata.n_obs == 0:
                update_status_fn(user_id, analysis_id, "No cells left after filtering. Please check your filter settings.",)
                raise ValueError("No cells left after filtering. Please check your filter settings.")

        # 2. Filter genes expressed in less than `filter_genes_value` cells
        if data_filtering.get("filter_genes"):
            current_app.logger.info(f"[UMAP] Filtering genes with min_cells={data_filtering.get('filter_genes_value')}")
            sc.pp.filter_genes(adata, min_cells=data_filtering.get("filter_genes_value"))
            current_app.logger.info(f"Shape after filtering genes: {adata.n_obs} cells × {adata.n_vars} genes")
            if adata.n_vars == 0:
                update_status_fn(user_id, analysis_id, "No genes left after filtering. Please check your filter settings.",)
                raise ValueError("No genes left after filtering. Please check your filter settings.")

        # 3. Filter on mitochondrial gene percentage (with species handling)
        if data_filtering.get("qc_filter"):
            current_app.logger.info(f"[UMAP] Filtering genes with min_counts={data_filtering.get('qc_filter_value')}")
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
                    mito_threshold = float(data_filtering.get("qc_filter_value"))
                    adata = adata[adata.obs.pct_counts_mt < mito_threshold, :]

                    # catastrophic case where all cells are removed
                    if adata.n_obs == 0:
                        current_app.logger.warning("  -> WARNING: All cells were removed by the mitochondrial filter. "
                                                   "Consider increasing the threshold or examining the data quality.")
                        update_status_fn(user_id, analysis_id, "No cells left after mitochondrial filtering.",)
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
            sc.pp.neighbors(adata, n_neighbors=params.get("n_neighbors"), metric=params.get("metric", "euclidean"))

            current_app.logger.info(f"[UMAP] Running UMAP with min_dist={params.get('min_dist')}")
            sc.tl.umap(adata, min_dist=params.get("min_dist"))

            # Save results
            result_path = analysis_data["results_path"]
            os.makedirs(result_path, exist_ok=True)
            umap_df = pd.DataFrame(
                adata.obsm["X_umap"], index=adata.obs_names, columns=["X_umap1", "X_umap2"]
            )
            if "X_pca" in adata.obsm:
                umap_df["X_pca1"] = adata.obsm["X_pca"][:, 0]
                umap_df["X_pca2"] = adata.obsm["X_pca"][:, 1]

            if "Cluster" in adata.obs:
                current_app.logger.info("[UMAP] Adding cluster labels from AnnData.obs.")
                umap_df["Cluster"] = adata.obs["Cluster"]
            else:
                current_app.logger.info("[UMAP] No cluster labels found, assigning default value 0.")
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
                fdr_level=analysis_data["inputs"].get("fdr_level"),
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
