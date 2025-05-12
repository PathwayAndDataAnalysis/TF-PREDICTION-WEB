import os

import pandas as pd
import scanpy as sc
from flask import current_app

from app import get_file_path


def run_umap_pipeline(
        user_id, analysis_id, analysis_data, users_data_path, update_status_fn
):
    try:
        current_app.logger.info(
            f"[UMAP] Starting UMAP pipeline for user '{user_id}', analysis '{analysis_id}'."
        )

        # 1. Load data
        gene_expr = analysis_data["inputs"]["gene_expression"]
        if gene_expr["source"] == "h5ad":
            current_app.logger.info(
                f"[UMAP] Loading .h5ad file: {get_file_path(gene_expr["h5ad_filename"], user_id)}"
            )
            adata = sc.read_h5ad(get_file_path(gene_expr["h5ad_filename"], user_id))
        else:
            current_app.logger.info(
                f"[UMAP] Loading counts file: {get_file_path(gene_expr['counts_filename'], user_id)}"
            )
            counts = pd.read_csv(get_file_path(gene_expr["counts_filename"], user_id), index_col=0)
            if gene_expr["metadata_filename"]:
                current_app.logger.info(
                    f"[UMAP] Loading metadata file: {get_file_path(gene_expr['metadata_filename'], user_id)}"
                )
                meta = pd.read_csv(get_file_path(gene_expr["metadata_filename"], user_id), index_col=0)
                adata = sc.AnnData(counts, obs=meta)
            else:
                current_app.logger.info(
                    "[UMAP] No metadata file provided, creating AnnData with counts only."
                )
                adata = sc.AnnData(counts)

        # 2. Run pipeline (use parameters from analysis_data["inputs"]["layout"]["umap_settings"])
        params = analysis_data["inputs"]["layout"]["umap_settings"] or {}
        current_app.logger.info(f"[UMAP] Pipeline parameters: {params}")

        if params.get("filter_cells"):
            current_app.logger.info(
                f"[UMAP] Filtering cells with min_genes={params.get('filter_cells_value', 500)}"
            )
            sc.pp.filter_cells(adata, min_genes=params.get("filter_cells_value", 500))

        if params.get("filter_genes"):
            current_app.logger.info(
                f"[UMAP] Filtering genes with min_cells={params.get('filter_genes_value', 100)}"
            )
            sc.pp.filter_genes(adata, min_cells=params.get("filter_genes_value", 100))

        if params.get("qc_filter"):
            current_app.logger.info(
                f"[UMAP] Filtering genes with min_counts={params.get('qc_filter_value', 10)}"
            )
            sc.pp.filter_genes(adata, min_counts=params.get("qc_filter_value", 10))

        if params.get("data_normalize"):
            current_app.logger.info(
                f"[UMAP] Normalizing total with target_sum={params.get('data_normalize_value', 10000)}"
            )
            sc.pp.normalize_total(
                adata, target_sum=params.get("data_normalize_value", 10000)
            )

        if params.get("log_transform"):
            current_app.logger.info("[UMAP] Applying log1p transformation.")
            sc.pp.log1p(adata)

        current_app.logger.info(
            f"[UMAP] Running PCA with n_comps={params.get('pca_components', 20)}"
        )
        sc.pp.pca(adata, n_comps=params.get("pca_components", 20))

        current_app.logger.info(
            f"[UMAP] Computing neighbors with n_neighbors={params.get('n_neighbors', 15)}, metric={params.get('metric', 'euclidean')}"
        )
        sc.pp.neighbors(
            adata,
            n_neighbors=params.get("n_neighbors", 15),
            metric=params.get("metric", "euclidean"),
        )
        current_app.logger.info(
            f"[UMAP] Running UMAP with min_dist={params.get('min_dist', 0.1)}"
        )
        sc.tl.umap(adata, min_dist=params.get("min_dist", 0.1))

        # 3. Save results
        result_path = analysis_data["results_path"]
        os.makedirs(result_path, exist_ok=True)
        umap_df = pd.DataFrame(
            adata.obsm["X_umap"], index=adata.obs_names, columns=["X_umap1", "X_umap2"]
        )
        if "Cluster" in adata.obs:
            current_app.logger.info("[UMAP] Adding cluster labels from AnnData.obs.")
            umap_df["Cluster"] = adata.obs["Cluster"]
        else:
            current_app.logger.info(
                "[UMAP] No cluster labels found, assigning default value 0."
            )
            umap_df["Cluster"] = 0  # or assign as needed

        umap_csv_path = os.path.join(result_path, "umap_coordinates.csv")
        current_app.logger.info(f"[UMAP] Saving UMAP coordinates to: {umap_csv_path}")
        umap_df.to_csv(umap_csv_path)

        # 4. Update status to Completed
        current_app.logger.info(
            f"[UMAP] UMAP pipeline completed successfully for analysis '{analysis_id}'."
        )
        update_status_fn(user_id, analysis_id, "Completed", umap_csv_path)
    except Exception as e:
        # Update status to Failed
        current_app.logger.error(
            f"[UMAP] UMAP pipeline failed for analysis '{analysis_id}': {e}",
            exc_info=True,
        )
        update_status_fn(user_id, analysis_id, "Failed", error=str(e))
