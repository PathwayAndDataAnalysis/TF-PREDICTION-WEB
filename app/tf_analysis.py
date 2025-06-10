import os

import math
import numpy as np
import pandas as pd
from flask import current_app
from joblib import Parallel, delayed  # For parallel processing
from scipy.sparse import issparse
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm  # For progress bar


###############################################################################
#         Benjamini-Hochberg FDR correction function                          #
###############################################################################


def bh_fdr_correction(
        p_value_df: pd.DataFrame, alpha: float = 0.05
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Performs Benjamini-Hochberg FDR correction on p-values from a file
    and returns both adjusted p-values and rejection status.

    Correction is applied independently to each column (e.g., each TF).
    """
    p_value_df.dropna(axis=1, how="all", inplace=True)

    df_adjusted_pvals = pd.DataFrame(
        index=p_value_df.index, columns=p_value_df.columns, dtype=float
    )  # For q-values
    df_reject_status = pd.DataFrame(
        index=p_value_df.index, columns=p_value_df.columns, dtype="boolean"
    )  # For True/False/NA

    for tf_name in p_value_df.columns:
        original_pvals_for_tf = p_value_df[tf_name].dropna()

        if original_pvals_for_tf.empty:
            continue

        reject_flags, corrected_pvals, _, _ = multipletests(
            pvals=original_pvals_for_tf.values,
            alpha=alpha,
            method="fdr_bh",
            is_sorted=False,
        )

        df_adjusted_pvals.loc[original_pvals_for_tf.index, tf_name] = corrected_pvals
        df_reject_status.loc[original_pvals_for_tf.index, tf_name] = reject_flags

    return df_adjusted_pvals, df_reject_status


###############################################################################
#                               Helper functions                              #
###############################################################################


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

    var = ((n_population + 1) * (n_population - k_sample)) / (
            12 * n_population ** 2 * k_sample
    )
    return math.sqrt(max(var, 0.0))


###############################################################################
#                      Per-cell processing function (for parallelization)     #
###############################################################################


def _process_single_cell(
        cell_name: str, expression_series: pd.Series, priors_df: pd.DataFrame
) -> tuple[str, dict[str, float]]:
    """
    Processes a single cell's expression data to infer TF activity.
    Returns the cell name and a dictionary of TF scores {regulator: p_value}.
    """
    # Drop genes with all‑nan expression, rank remaining genes
    expr = expression_series.dropna().sort_values(ascending=False)
    n_genes = expr.size

    # Robustness: If no genes with valid expression, return empty scores for this cell.
    if n_genes == 0:
        return cell_name, {}

    ranks = (np.arange(1, n_genes + 1) - 0.5) / n_genes
    rank_df = pd.DataFrame({"TargetGene": expr.index, "Rank": ranks})

    # Merge ranks into a *copy* of priors_df (avoid in‑place edits)
    # .merge() creates a new DataFrame, so priors_df itself is not modified.
    p = priors_df.merge(rank_df, on="TargetGene", how="left")
    p = p[p["Rank"].notna()]

    # Invert rank for repressors (RegulatoryEffect == 0)
    max_rank_val = (n_genes - 0.5) / n_genes
    p["AdjustedRank"] = np.where(
        p["RegulatoryEffect"].eq(0),
        max_rank_val - p["Rank"],
        p["Rank"],
    )  # AdjustedRank is 1 - Rank for repressors

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
    # Create new column, If RankMean is <0.5, 1 else -1
    tf_summary["PValueDir"] = np.where(
        tf_summary["RankMean"] < 0.5, 1, -1
    )

    # If tf_summary is empty after filtering, no TFs to score for this cell
    if tf_summary.empty:
        return cell_name, {}

    # Compute σₙ,ₖ, Z‑score, and p‑value
    tf_summary["Sigma_n_k"] = tf_summary["AvailableTargets"].apply(
        lambda k: std_dev_mean_norm_rank(n_genes, k)
    )
    tf_summary["Z"] = (tf_summary["RankMean"] - 0.5) / tf_summary["Sigma_n_k"].replace(
        0, np.nan
    )
    tf_summary["P_two_tailed"] = np.where(
        tf_summary["Z"].abs() < 0.5,
        2 * norm.cdf(tf_summary["Z"].abs()),
        2 * norm.sf(tf_summary["Z"].abs()),
    )
    # Add PValueDir information to the P_two_tailed column
    tf_summary["P_two_tailed"] *= tf_summary["PValueDir"]

    # Create a dictionary of scores for this cell
    cell_scores_dict = pd.Series(
        tf_summary["P_two_tailed"].values, index=tf_summary["Regulator"]
    ).to_dict()

    return cell_name, cell_scores_dict


def run_tf_analysis(
        user_id, analysis_id, analysis_data, adata, update_analysis_status_fn
):
    try:
        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="Running TF analysis"
        )
        n_jobs = 4
        current_app.logger.info(
            f"[TF_ANALYSIS] Running TF analysis for user '{user_id}', analysis '{analysis_id}'."
        )

        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="TF analysis is running"
        )
        current_app.logger.info(f"[TF_ANALYSIS] Analysis data: {analysis_data}")

        # ───── Load gene expression data ────────────────────────────────────────
        X = adata.X.toarray() if issparse(adata.X) else adata.X.copy()
        # Treat explicit zero counts as missing for z‑scoring consistency
        X[X == 0] = np.nan
        print("Calculating Z-scores...")
        z_mat = zscore(X, axis=1, nan_policy="omit")
        z_df = pd.DataFrame(z_mat, index=adata.obs_names, columns=adata.var_names)

        z_scores_path = os.path.join(
            analysis_data.get("results_path", ""), "z_scores.tsv"
        )
        print(f"Saving Z-scores to {z_scores_path}...")
        z_df.to_csv(z_scores_path, sep="\t", index=True)

        # ───── Load TF‑target priors ──────────────────────────────────────────────
        print("Loading TF-target priors...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        priors = pd.read_csv(
            os.path.join(script_dir, "..", "prior_data", "causal_priors.tsv"), sep="\t"
        )
        expected_cols = {"Regulator", "TargetGene", "RegulatoryEffect"}
        if not expected_cols.issubset(priors.columns):
            raise ValueError(f"Prior file must contain columns {expected_cols}")

        # 2. Run Analysis
        current_app.logger.info(
            f"[TF_ANALYSIS] Starting TF analysis for user '{user_id}', analysis '{analysis_id}'."
        )

        # Pre-define columns for the result DataFrame
        all_regulators = priors["Regulator"].unique()
        result = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)

        # ───── Parallel processing of cells ───────────────────────────────────────
        num_cells = len(z_df)
        if num_cells == 0:
            print("No cells to process.")
        else:
            print(
                f"Starting TF activity inference for {num_cells} cells using {n_jobs if n_jobs > 0 else 'all available'} cores..."
            )

            # Prepare arguments for each parallel task
            # Each task processes one cell: (cell_name, expression_series_for_that_cell, priors_dataframe)
            # We pass a copy of priors to each worker; it's read-only within the worker.
            # z_df.loc[cell_name] gives the expression Series for that cell.
            tasks = [
                delayed(_process_single_cell)(
                    cell_name,
                    z_df.loc[cell_name],  # Pass the expression series for the cell
                    priors,  # Pass the (entire) priors DataFrame
                )
                for cell_name in z_df.index
            ]

            # Run tasks in parallel
            # `tqdm` provides a progress bar
            # `backend="loky"` is robust and default for joblib in many cases, good for cross-platform
            if n_jobs == 1:  # Useful for debugging
                print("Running sequentially (n_jobs=1)...")
                cell_results_list = [
                    _process_single_cell(cell_name, z_df.loc[cell_name], priors)
                    for cell_name in tqdm(z_df.index, desc="Processing cells")
                ]
            else:
                print(f"Running in parallel with n_jobs={n_jobs}...")
                cell_results_list = Parallel(n_jobs=n_jobs, backend="loky", verbose=5)(
                    tqdm(tasks, desc="Processing cells")
                )

            # ───── Aggregate results ───────────────────────────────────────────────
            print("\nAggregating results...")
            for cell_name, scores_dict in tqdm(
                    cell_results_list, desc="Aggregating scores"
            ):
                if scores_dict:  # If the dictionary is not empty
                    for regulator, p_value in scores_dict.items():
                        if regulator in result.columns:
                            result.at[cell_name, regulator] = p_value
                        # else:
                        # This case should ideally not happen if all_regulators is correctly derived
                        # print(f"Warning: Regulator '{regulator}' for cell '{cell_name}' not in pre-defined columns. Skipping.")

        # ───── Save results ───────────────────────────────────────────────────────
        result_path = analysis_data.get("results_path", None)
        p_value_sign = result.copy()
        p_value_sign = p_value_sign.map(
            lambda x: 1 if x > 0 else (-1 if x < 0 else np.nan)
        )
        result = result.map(
            lambda x: abs(x) if pd.notna(x) else np.nan
        )
        p_values_path = os.path.join(result_path, "p_values.tsv")
        print(f"Saving results to {p_values_path}...")
        result.to_csv(p_values_path, sep="\t", index=True)

        # ───── Run Benjamini Hotchberg FDR Correction ────────────────────────────────
        update_analysis_status_fn(
            user_id=user_id,
            analysis_id=analysis_id,
            status="Running BH FDR correction",
            # tfs=result.columns.tolist(),
            pvalues_path=p_values_path,
        )
        print("Running Benjamini-Hochberg FDR correction...")
        _, reject = bh_fdr_correction(result)
        # Replace True with 1 and False with -1 others with 0
        reject = reject.map(
            lambda x: 1 if x is True else (4 if x is False else 5)
        )
        reject = reject.multiply(p_value_sign, axis=0)
        reject = reject.replace({4: 0, 5: np.nan, -4: 0, -5: np.nan}).astype("Int64")
        # 1 -> Activated, -1 -> Inhibited, 0 -> Not significant, and NaN -> Not Enough Data
        bh_reject_path = os.path.join(result_path, "bh_reject.tsv")
        print(f"Saving FDR results to {bh_reject_path}...")
        reject.to_csv(bh_reject_path, sep="\t", index=True)

        # Count how many Trues for each TF and sort tfs by this count
        tf_counts = reject.sum().sort_values(ascending=False)

        update_analysis_status_fn(
            user_id=user_id,
            analysis_id=analysis_id,
            status="Completed",
            tfs=tf_counts.index.tolist(),
            bh_reject_path=bh_reject_path,
        )
        current_app.logger.info(
            f"[TF_ANALYSIS] TF analysis completed for user '{user_id}', analysis '{analysis_id}'."
        )

    except Exception as e:
        current_app.logger.error(
            f"[TF_ANALYSIS] Error in run_tf_analysis for user '{user_id}', analysis '{analysis_id}': {e}"
        )
        update_analysis_status_fn(
            user_id=user_id, analysis_id=analysis_id, status="Error in TF analysis"
        )
        raise e
