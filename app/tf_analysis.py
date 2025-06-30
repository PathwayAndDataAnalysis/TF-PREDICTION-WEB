import os
import threading
import time
import math
import numpy as np
import pandas as pd
from flask import current_app
from joblib import Parallel, delayed
from scipy.sparse import issparse
from scipy.stats import zscore, norm
from statsmodels.stats.multitest import multipletests
from tqdm.auto import tqdm
import psutil
import gc

# CORES_USED = 1 # For debugging, use a single core
CORES_USED = max(1, int(os.cpu_count() * 0.8))  # Use 80% of available cores


def check_memory_usage():
    """Check if we have enough memory for processing."""
    memory = psutil.virtual_memory()
    if memory.percent > 90:  # If more than 90% memory is used
        gc.collect()  # Force garbage collection
        memory = psutil.virtual_memory()
        if memory.percent > 90:
            raise MemoryError("Insufficient memory for processing. Please try with smaller data.")


###############################################################################
#         Benjamini-Hochberg FDR correction function                          #
###############################################################################


def bh_fdr_correction(
        p_value_df: pd.DataFrame, alpha: float
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
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
    # Compute threshold per TF
    thresholds = {}

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
        # Store the threshold for this TF
        thresholds[tf_name] = original_pvals_for_tf[reject_flags].max() if reject_flags.any() else np.nan

        df_adjusted_pvals.loc[original_pvals_for_tf.index, tf_name] = corrected_pvals
        df_reject_status.loc[original_pvals_for_tf.index, tf_name] = reject_flags

    thresholds_seris = pd.Series(thresholds, dtype=float)
    return df_adjusted_pvals, df_reject_status, thresholds_seris


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
):
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
    tf_summary["ActivationDir"] = np.where(
        tf_summary["RankMean"] < 0.5, 1, -1
    )
    # Compute σₙ,ₖ, Z‑score, and p‑value
    tf_summary["Sigma_n_k"] = tf_summary["AvailableTargets"].apply(
        lambda k: std_dev_mean_norm_rank(n_genes, k)
    )
    tf_summary["Z"] = (tf_summary["RankMean"] - 0.5) / tf_summary["Sigma_n_k"].replace(
        0, np.nan
    )
    tf_summary["P_two_tailed"] = np.where(
        tf_summary["RankMean"] < 0.5,
        2 * norm.cdf(tf_summary["Z"]),
        2 * (1 - norm.cdf(tf_summary["Z"])),
    )
    # tf_summary["P_two_tailed"] = 2 * norm.sf(tf_summary["Z"])

    # Create a dictionary where each value is a tuple: (p_value, direction)
    cell_scores_dict = {
        regulator: (p_value, direction)
        for regulator, p_value, direction in zip(
            tf_summary["Regulator"],
            tf_summary["P_two_tailed"],
            tf_summary["ActivationDir"],
        )
    }
    return cell_name, cell_scores_dict


def run_tf_analysis(
        user_id, analysis_id, analysis_data, adata, fdr_level, update_analysis_status_fn
):
    try:
        # Check memory before starting
        check_memory_usage()

        update_analysis_status_fn(user_id=user_id, analysis_id=analysis_id, status="Running analysis")

        current_app.logger.info(f"[TF_ANALYSIS] Running analysis for user '{user_id}', analysis '{analysis_id}'.")
        current_app.logger.info(f"[TF_ANALYSIS] Analysis data: {analysis_data}")

        # Validate input data
        if adata is None or adata.n_obs == 0:
            raise ValueError("Input data is empty or invalid")
        if adata.n_obs > 150000:  # Limit to 150k cells
            current_app.logger.warning(f"Large dataset detected: {adata.n_obs} cells. This may take a long time.")

        # Add memory monitoring during processing
        def monitor_memory():
            while True:
                check_memory_usage()
                time.sleep(30)  # Check every 30 seconds

        # Start memory monitoring in background
        memory_thread = threading.Thread(target=monitor_memory, daemon=True)
        memory_thread.start()

        try:
            # ───── Load gene expression data ────────────────────────────────────────
            X = adata.X.toarray() if issparse(adata.X) else adata.X.copy()
            X[X == 0] = np.nan
            print("Calculating Z-scores...")
            z_mat = zscore(X, axis=0, nan_policy="omit") # across all cells for a single gene (axis=0)
            z_df = pd.DataFrame(z_mat, index=adata.obs_names, columns=adata.var_names)

            z_scores_path = os.path.join(analysis_data.get("results_path", ""), "z_scores.csv")
            print(f"Saving Z-scores to {z_scores_path}...")
            z_df.to_csv(z_scores_path, index=True)

            # ───── Load TF‑target priors ──────────────────────────────────────────────
            print("Loading TF-target priors...")
            script_dir = os.path.dirname(os.path.abspath(__file__))
            priors = pd.read_csv(os.path.join(script_dir, "..", "prior_data", "causal_priors.tsv"), sep="\t")
            priors = priors[priors["TargetGene"].notna()]
            expected_cols = {"Regulator", "TargetGene", "RegulatoryEffect"}
            if not expected_cols.issubset(priors.columns):
                raise ValueError(f"Prior file must contain columns {expected_cols}")

            # Run Analysis
            current_app.logger.info(f"[TF_ANALYSIS] Starting TF analysis for user '{user_id}', analysis '{analysis_id}'.")
            # Pre-define columns for the p_values_df DataFrame
            all_regulators = priors["Regulator"].unique()
            p_values_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)
            activation_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)

            # ───── Parallel processing of cells ───────────────────────────────────────
            num_cells = len(z_df)
            if num_cells == 0:
                print("No cells to process.")
            else:
                print(
                    f"Starting TF activity for {num_cells} cells using {CORES_USED if CORES_USED > 0 else 'all available'} cores..."
                )
                tasks = [
                    delayed(_process_single_cell)(
                        cell_name,
                        z_df.loc[cell_name],  # Pass the expression series for the cell
                        priors,  # Pass the (entire) priors DataFrame
                    )
                    for cell_name in z_df.index
                ]

                # Run tasks in parallel, `tqdm` provides a progress bar
                # `backend="loky"` is robust and default for joblib in many cases, good for cross-platform
                if CORES_USED == 1:  # Useful for debugging
                    print("Running sequentially (CORES_USED=1)...")
                    cell_results_list = [
                        _process_single_cell(cell_name, z_df.loc[cell_name], priors)
                        for cell_name in tqdm(z_df.index, desc="Processing cells")
                    ]
                else:
                    print(f"Running in parallel with CORES_USED={CORES_USED}...")
                    cell_results_list = Parallel(n_jobs=CORES_USED, backend="loky", verbose=2)(
                        tqdm(tasks, desc="Processing cells")
                    )

                # ───── Aggregate results into two separate DataFrames ───────────────────
                print("\nAggregating results...")
                for cell_name, scores_dict in tqdm(
                        cell_results_list, desc="Aggregating scores"
                ):
                    if scores_dict:
                        # Unpack the (p_value, direction) tuple for each regulator
                        for regulator, (p_value, direction) in scores_dict.items():
                            if regulator in p_values_df.columns:
                                p_values_df.at[cell_name, regulator] = p_value
                                activation_df.at[cell_name, regulator] = direction

            # ───── Save p_value and activation results ─────────────────────────────────────
            result_path = analysis_data.get("results_path", None)

            p_values_path = os.path.join(result_path, "p_values.csv")
            print(f"Saving p-values to {p_values_path}...")
            p_values_df.to_csv(p_values_path, index=True)
            activation_path = os.path.join(result_path, "activation.csv")
            print(f"Saving activation results to {activation_path}...")
            activation_df.to_csv(activation_path, index=True)

            # ───── Run Benjamini Hochberg FDR Correction ────────────────────────────────
            update_analysis_status_fn(
                user_id=user_id,
                analysis_id=analysis_id,
                status="Running BH FDR correction",
                pvalues_path=p_values_path,
                activation_path=activation_path,
            )
            print("Running Benjamini-Hochberg FDR correction...")
            _, reject, p_val_thresholds_seris = bh_fdr_correction(p_value_df=p_values_df, alpha=fdr_level)
            tf_counts = reject.sum().sort_values(ascending=False)

            final_output = pd.DataFrame(0, index=reject.index, columns=reject.columns)
            # Where the rejection is True, fill in the direction from activation_df
            final_output[reject] = activation_df[reject]

            # Propagate NaNs for cells/TFs where analysis couldn't be run
            final_output[p_values_df.isna()] = np.nan
            final_output = final_output.astype("Int64")  # Use nullable integer type

            bh_reject_path = os.path.join(result_path, "bh_reject.csv")
            print(f"Saving FDR results to {bh_reject_path}...")
            final_output.to_csv(bh_reject_path, index=True)

            p_val_threshold_path = os.path.join(result_path, "p_val_thresholds.csv")
            print(f"Saving p-value thresholds to {p_val_threshold_path}...")
            p_val_thresholds_seris.to_csv(p_val_threshold_path, index=True)

            update_analysis_status_fn(
                user_id=user_id,
                analysis_id=analysis_id,
                status="Completed",
                tfs=tf_counts.index.tolist() if(tf_counts.index.tolist()) else [],
                bh_reject_path=bh_reject_path,
                fdr_level=fdr_level,
                p_val_threshold_path=p_val_threshold_path,
                z_scores_path=z_scores_path
            )
        finally:
            # Stop memory monitoring
            memory_thread.join(timeout=1)

        # ───── Final logging and status update ─────────────────────────────────────
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