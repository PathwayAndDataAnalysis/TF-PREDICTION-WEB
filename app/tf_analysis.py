import os
import math
import numpy as np
import pandas as pd
from flask import current_app
from joblib import Parallel, delayed
from scipy.sparse import issparse
from scipy.stats import zscore, norm
from tqdm.auto import tqdm
from app.benjamini_hotchberg import bh_fdr_correction


# CORES_USED = 1 # For debugging, use a single core
CORES_USED = max(1, int(os.cpu_count() * 0.8))  # Use 80% of available cores


def std_dev_mean_norm_rank(n_population: int, k_sample: int) -> float:
    """Return σₙ,ₖ – the analytic SD of the mean of *k* normalised ranks.
    Calculates sigma_m_k: the theoretical standard deviation of the mean of k
    normalized ranks (r_i = (R_i - 0.5)/n) sampled without replacement
    from a population of size n.
    """
    if not (isinstance(n_population, int) and n_population > 0):
        raise ValueError("n_population must be a positive integer.")
    if not (isinstance(k_sample, int) and k_sample > 0):
        raise ValueError("k_sample must be a positive integer.")
    if k_sample > n_population:
        raise ValueError(
            "Sample size k_sample cannot exceed population size n_population."
        )

    if k_sample == n_population or n_population == 1:
        return 0.0  # No variability – we sampled everything.

    var = ((n_population + 1) * (n_population - k_sample)) / (12 * n_population ** 2 * k_sample)
    return math.sqrt(max(var, 0.0))


def _process_single_cell(cell_name: str, expression_series: pd.Series, priors_df: pd.DataFrame):
    expr = expression_series.dropna().sort_values(ascending=False)
    n_genes = expr.size

    # Robustness: If no genes with valid expression, return empty scores for this cell.
    if n_genes == 0:
        return cell_name, {}

    ranks = (np.arange(1, n_genes + 1) - 0.5) / n_genes
    rank_df = pd.DataFrame({"TargetGene": expr.index, "Rank": ranks})

    # Merge ranks into a *copy* of priors_df (avoid in‑place edits)
    p = priors_df.merge(rank_df, on="TargetGene", how="left")
    p = p[p["Rank"].notna()]

    # Invert rank for repressors (RegulatoryEffect == 0)
    max_rank_val = (n_genes - 0.5) / n_genes
    p["AdjustedRank"] = np.where(
        p["RegulatoryEffect"].eq(0),
        max_rank_val - p["Rank"],
        p["Rank"],
    )  # AdjustedRank is max_rank_val - Rank for repressors

    # Summarise per TF
    tf_summary = (
        p.groupby("Regulator")
        .agg(
            AvailableTargets=("AdjustedRank", lambda x: x.notna().sum()),
            RankMean=("AdjustedRank", "mean"),
        )
        .reset_index()
    )
    tf_summary = tf_summary[tf_summary["AvailableTargets"] > 0]

    # If tf_summary is empty after filtering, no TFs to score for this cell
    if tf_summary.empty:
        return cell_name, {}

    # Create new column, If RankMean is <0.5, 1 else -1
    tf_summary["ActivationDir"] = np.where(tf_summary["RankMean"] < 0.5, 1, -1)

    # Compute σₙ,ₖ, Z‑score, and p‑value
    tf_summary["Sigma_n_k"] = tf_summary["AvailableTargets"].apply(
        lambda k: std_dev_mean_norm_rank(n_genes, k)
    )
    tf_summary["Z"] = (tf_summary["RankMean"] - 0.5) / tf_summary["Sigma_n_k"].replace(0, np.nan)
    tf_summary["P_two_tailed"] = np.where(
        tf_summary["RankMean"] < 0.5,
        2 * norm.cdf(tf_summary["Z"]),
        2 * (1 - norm.cdf(tf_summary["Z"])),
    )
    # tf_summary["P_two_tailed"] = 2 * norm.sf(tf_summary["Z"])

    # Create a dictionary where each value is a tuple: (p_value, direction)
    # cell_scores_dict = {
    #     regulator: (p_value, direction)
    #     for regulator, p_value, direction in zip(
    #         tf_summary["Regulator"],
    #         tf_summary["P_two_tailed"],
    #         tf_summary["ActivationDir"],
    #     )
    # }
    # return cell_name, cell_scores_dict
    # Instead of a dictionary, return lists of the components.
    # This makes aggregation much faster later.
    regulators = tf_summary["Regulator"].tolist()
    p_values = tf_summary["P_two_tailed"].tolist()
    directions = tf_summary["ActivationDir"].tolist()

    return cell_name, regulators, p_values, directions


def run_tf_analysis(
        user_id, analysis_id, analysis_data, adata, fdr_level, update_analysis_status_fn
):
    try:
        update_analysis_status_fn(user_id=user_id, analysis_id=analysis_id, status="Running analysis")

        current_app.logger.info(f"[TF_ANALYSIS] Analysis data: {analysis_data}")

        # Validate input data
        if adata is None or adata.n_obs == 0:
            raise ValueError("Input data is empty or invalid")
        if adata.n_obs > 200000:  # Limit to 200k cells
            current_app.logger.warning(f"Large dataset detected: {adata.n_obs} cells. This may take a long time.")


        # ────────── Load gene expression data ────────────────────────────────────────
        X = adata.X.toarray() if issparse(adata.X) else adata.X.copy()
        X[X == 0] = np.nan
        X[X == 0.0] = np.nan
        print("Calculating Z-scores...")
        z_mat = zscore(X, axis=0, nan_policy="omit")  # across all cells for a single gene (axis=0)
        adata.layers["z_scores"] = z_mat  # Store Z-scores in AnnData layers
        z_df = pd.DataFrame(z_mat, index=adata.obs_names, columns=adata.var_names)


        # ────────── Load TF‑target priors ──────────────────────────────────────────────
        print("Loading TF-target priors...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        priors = pd.read_csv(os.path.join(script_dir, "..", "prior_data", "causal_priors.tsv"), sep="\t")
        priors = priors[priors["TargetGene"].notna()]
        # Run Analysis
        current_app.logger.info(f"[TF_ANALYSIS] Starting TF analysis for user '{user_id}', analysis '{analysis_id}'.")
        # Pre-define columns for the p_values_df DataFrame
        all_regulators = priors["Regulator"].unique()
        # p_values_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)
        # activation_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)


        # ────────── Parallel processing of cells ──────────
        print(f"Starting analysis for {len(z_df)} cells using {CORES_USED if CORES_USED > 0 else 'all'} cores.")
        tasks = [
            delayed(_process_single_cell)(
                cell_name,
                z_df.loc[cell_name],  # Pass the expression series for the cell
                priors,  # Pass the (entire) priors DataFrame
            )
            for cell_name in z_df.index
        ]

        # `backend="loky"` is robust and default for joblib in many cases, good for cross-platform
        print(f"Running in parallel with CORES_USED={CORES_USED}...")
        cell_results_list = Parallel(n_jobs=CORES_USED, backend="loky", verbose=2)(
            tqdm(tasks, desc="Processing cells")
        )

        # ────────── Aggregate results into two separate DataFrames ──────────
        # print("\nAggregating results...")
        # for cell_name, scores_dict in tqdm(cell_results_list, desc="Aggregating scores"):
        #     if scores_dict:
        #         for regulator, (p_value, direction) in scores_dict.items():
        #             if regulator in p_values_df.columns:
        #                 p_values_df.at[cell_name, regulator] = p_value
        #                 activation_df.at[cell_name, regulator] = direction
        # --- NEW: High-Performance Aggregation ---
        print("\nAggregating results...")
        p_value_records = []
        activation_records = []
        # 1. Unpack results into flat lists of records
        for cell_name, regulators, p_values, directions in tqdm(cell_results_list, desc="Unpacking results"):
            for i, regulator in enumerate(regulators):
                p_value_records.append({'cell': cell_name, 'regulator': regulator, 'p_value': p_values[i]})
                activation_records.append({'cell': cell_name, 'regulator': regulator, 'direction': directions[i]})
        # 2. Build temporary DataFrames from the lists of records (very fast)
        p_values_temp_df = pd.DataFrame.from_records(p_value_records)
        activation_temp_df = pd.DataFrame.from_records(activation_records)
        # 3. Pivot the temporary DataFrames into the desired final shape (cells x regulators)
        # This is the most efficient way to reshape the data.
        if not p_values_temp_df.empty:
            p_values_df = p_values_temp_df.pivot(index='cell', columns='regulator', values='p_value')
            activation_df = activation_temp_df.pivot(index='cell', columns='regulator', values='direction')
        else:
            # Handle case where no TFs were found for any cell
            all_regulators = priors["Regulator"].unique()
            p_values_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)
            activation_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)
        # --- END OF NEW AGGREGATION ---
        current_app.logger.info("[TF_ANALYSIS] Reindexing result DataFrames to match AnnData object order.")
        p_values_df = p_values_df.reindex(adata.obs_names)
        activation_df = activation_df.reindex(adata.obs_names)


        # ────────── Store results in the AnnData object ──────────
        p_values_df.dropna(axis=1, how="all", inplace=True)
        activation_df = activation_df[p_values_df.columns]
        adata.obsm['p_values'] = p_values_df
        adata.obsm['activation'] = activation_df


        # ────────── Run BH correction and store results in the AnnData object ──────────
        update_analysis_status_fn(user_id=user_id, analysis_id=analysis_id, status="Running BH FDR correction")
        print("Running Benjamini-Hochberg FDR correction...")
        _, reject, p_val_thresholds_series = bh_fdr_correction(p_value_df=p_values_df, alpha=fdr_level)
        final_output = pd.DataFrame(0, index=reject.index, columns=reject.columns)
        final_output[reject] = activation_df[reject]
        final_output[p_values_df.isna()] = np.nan
        final_output = final_output.astype("Int64")  # Use a nullable integer type
        adata.obsm['bh_reject'] = final_output
        adata.uns['p_val_thresholds'] = p_val_thresholds_series.to_dict()
        tf_counts = reject.sum().sort_values(ascending=False)


        # ────────── FINAL STEP: Save the single, enriched AnnData object as the final artifact ──────────
        result_path = analysis_data.get("results_path")
        final_artifact_path = os.path.join(result_path, "analysis_results.h5ad")
        adata.write_h5ad(final_artifact_path, compression="gzip")
        current_app.logger.info(f"[TF_ANALYSIS] All results saved to single artifact: {final_artifact_path}")


        # ────────── Update analysis status ──────────
        update_analysis_status_fn(
            user_id=user_id,
            analysis_id=analysis_id,
            status="Completed",
            artifact_path=final_artifact_path,
            tfs=tf_counts.index.tolist()
        )

        current_app.logger.info(f"[TF_ANALYSIS] TF analysis completed for user '{user_id}', analysis '{analysis_id}'")

    except MemoryError as e:
        current_app.logger.error(f"[TF_ANALYSIS] Memory error: {e}")
        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="Error: Insufficient memory"
        )
        raise e

    except TimeoutError as e:
        current_app.logger.error(f"[TF_ANALYSIS] Timeout error: {e}")
        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="Error: Analysis timed out"
        )
        raise e

    except Exception as e:
        current_app.logger.error(
            f"[TF_ANALYSIS] Error in run_tf_analysis for user '{user_id}', analysis '{analysis_id}': {e}"
        )
        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="Error in TF analysis"
        )
        raise e