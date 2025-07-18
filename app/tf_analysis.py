import os
import math
import numpy as np
import pandas as pd
from flask import current_app
from joblib import Parallel, delayed
from scipy.sparse import issparse
from scipy.stats import zscore, norm
from tqdm.auto import tqdm

from app.benjamini_hotchberg import run_bh_correction_and_save_tfs

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

    var = ((n_population + 1) * (n_population - k_sample)) / (
            12 * n_population ** 2 * k_sample
    )
    return math.sqrt(max(var, 0.0))


def _process_single_cell(
        cell_name: str, expression_series: pd.Series, priors_df: pd.DataFrame
):
    # Drop genes with all‑nan expression, rank remaining genes
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
        p.groupby("Regulator", observed=True)
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

    regulators = tf_summary["Regulator"].tolist()
    p_values = tf_summary["P_two_tailed"].tolist()
    directions = tf_summary["ActivationDir"].tolist()

    return cell_name, regulators, p_values, directions


def run_tf_analysis(user_id, analysis_id, analysis_data, adata, update_analysis_status_fn):
    try:
        update_analysis_status_fn(user_id=user_id, analysis_id=analysis_id, status="Running analysis")
        current_app.logger.info(f"[TF_ANALYSIS] Running analysis for user '{user_id}', analysis data '{analysis_data}'.")

        # Validate input data
        if adata is None or adata.n_obs == 0:
            raise ValueError("Input data is empty or invalid")
        if adata.n_obs > 300000:  # Limit to 300k cells
            current_app.logger.warning(f"Large dataset detected: {adata.n_obs} cells. This may take a long time.")

        # ───── Load gene expression data ────────────────────────────────────────
        X = adata.X.toarray() if issparse(adata.X) else adata.X.copy()
        X[X == 0] = np.nan
        X[X == 0.0] = np.nan
        print("Calculating Z-scores...")
        z_mat = zscore(X, axis=0, nan_policy="omit") # across all cells for a single gene (axis=0)
        z_df = pd.DataFrame(z_mat, index=adata.obs_names, columns=adata.var_names)

        # Replace the slow .to_csv() with this:
        z_scores_path = os.path.join(analysis_data.get("results_path", ""), "z_scores.parquet")
        print(f"Saving Z-scores to {z_scores_path}...")
        z_df.to_parquet(z_scores_path)


        # ───── Load TF‑target priors ──────────────────────────────────────────────
        print("Loading TF-target priors...")
        script_dir = os.path.dirname(os.path.abspath(__file__))
        causal_priors_path = (
            os.path.join(script_dir, "..", "prior_data", "causal_priors.tsv")
            if analysis_data.get("inputs", {}).get("prior_data", {}).get("prior_data_filepath") == "Default"
            else analysis_data.get("inputs", {}).get("prior_data", {}).get("prior_data_filepath")
        )

        effect_map = {"upregulates-expression": 1, "downregulates-expression": 0}
        col_names = ["Regulator", "RegulatoryEffect", "TargetGene"]
        priors = (
            pd.read_csv(
                causal_priors_path,
                sep="\t",
                header=None,
                usecols=[0, 1, 2],
                names=col_names,
            )
            .assign(RegulatoryEffect=lambda df: df["RegulatoryEffect"].map(effect_map))
            .dropna(subset=["RegulatoryEffect"])
            .astype({
                "RegulatoryEffect": "int8",
                "Regulator": "category",
                "TargetGene": "category"
            })
            .reset_index(drop=True)
        )
        # TODO: Add minium number of targets per regulator code
        min_number_of_targets = analysis_data.get("inputs", {}).get("prior_data", {}).get("min_number_of_targets", 3)
        # Only keep regulators with at least `min_number_of_targets` targets
        priors = priors.groupby("Regulator").filter(lambda x: len(x) >= min_number_of_targets)
        # Only keep the rows of priors that have targets in the z_df (z_df.columns)
        priors = priors[priors["TargetGene"].isin(z_df.columns)]
        priors_group = priors.groupby("Regulator", observed=True).agg({"RegulatoryEffect": list, "TargetGene": list})
        # Now remove the rows the length of the lists in RegulatoryEffect is less than min_number_of_targets
        priors_group = priors_group[priors_group["RegulatoryEffect"].apply(len) >= min_number_of_targets]
        # Keep only regulators in priors that is in priors_group.index
        priors = priors[priors["Regulator"].isin(priors_group.index)]

        # ───── Parallel processing of cells ───────────────────────────────────────
        current_app.logger.info(f"[TF_ANALYSIS] Starting TF analysis for user '{user_id}', analysis '{analysis_id}'.")
        print(f"Starting TF activity using {CORES_USED if CORES_USED > 0 else 'all available'} cores.")
        tasks = [
            delayed(_process_single_cell)(
                cell_name,
                z_df.loc[cell_name],  # Pass the expression series for the cell
                priors,  # Pass the (entire) priors DataFrame
            )
            for cell_name in z_df.index
        ]

        # `backend="loky"` is robust and default for joblib in many cases, good for cross-platform
        if CORES_USED == 1:  # Useful for debugging
            print("Running sequentially (CORES_USED=1)...")
            cell_results_list = [
                _process_single_cell(cell_name, z_df.loc[cell_name], priors)
                for cell_name in tqdm(z_df.index, desc="Processing cells")
            ]
        else:
            print(f"Running in parallel with CORES_USED={CORES_USED}...")
            cell_results_list = Parallel(n_jobs=CORES_USED, backend="loky", verbose=2)(tqdm(tasks, desc="Processing cells"))


        # ───── Aggregate results into two separate DataFrames ───────────────────
        print("\nAggregating results...")
        p_value_records = []
        activation_records = []
        # 1. Unpack results into flat lists of records
        for cell_name, regulators, p_values, directions in tqdm(cell_results_list, desc="Unpacking results"):
            for i, regulator in enumerate(regulators):
                p_value_records.append({'cell': cell_name, 'regulator': regulator, 'p_value': p_values[i]})
                activation_records.append(
                    {'cell': cell_name, 'regulator': regulator, 'direction': directions[i]})
        # 2. Build temporary DataFrames from the lists of records (very fast)
        p_values_temp_df = pd.DataFrame.from_records(p_value_records)
        activation_temp_df = pd.DataFrame.from_records(activation_records)
        # 3. Pivot the temporary DataFrames into the desired final shape (cells x regulators)
        # This is the most efficient way to reshape the data.
        if not p_values_temp_df.empty:
            p_values_df = p_values_temp_df.pivot(index='cell', columns='regulator', values='p_value')
            activation_df = activation_temp_df.pivot(index='cell', columns='regulator', values='direction')
        else:
            # Handle a case where no TFs were found for any cell
            all_regulators = priors["Regulator"].unique()
            p_values_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)
            activation_df = pd.DataFrame(index=z_df.index, columns=all_regulators, dtype=float)

        current_app.logger.info("[TF_ANALYSIS] Reindexing result DataFrames to match AnnData object order.")
        p_values_df = p_values_df.reindex(adata.obs_names)
        activation_df = activation_df.reindex(adata.obs_names)


        # ───── Save p_value and activation results ─────────────────────────────────────
        result_path = analysis_data.get("results_path", None)
        p_values_path = os.path.join(result_path, "p_values.parquet")
        print(f"Saving p-values to {p_values_path}...")
        p_values_df.dropna(axis=1, how="all", inplace=True)
        p_values_df.to_parquet(p_values_path)

        activation_path = os.path.join(result_path, "activation.parquet")
        print(f"Saving activation results to {activation_path}...")
        activation_df = activation_df[p_values_df.columns]
        activation_df.to_parquet(activation_path)


        # ───── Run Benjamini Hochberg FDR Correction ────────────────────────────────
        run_bh_correction_and_save_tfs(
            user_id=user_id,
            analysis_id=analysis_id,
            p_values_df=p_values_df,
            update_analysis_status_fn=update_analysis_status_fn,
            p_values_path=p_values_path,
            activation_path=activation_path,
            z_scores_path=z_scores_path,
        )

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
        update_analysis_status_fn(user_id=user_id, analysis_id=analysis_id, status="Error in TF analysis")
        raise e