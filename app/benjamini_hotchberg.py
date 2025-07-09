import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

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
    # Compute a threshold per TF
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


def run_bh_and_save_files(
        user_id,
        analysis_id,
        p_values_df: pd.DataFrame,
        activation_df: pd.DataFrame,
        result_path: str,
        fdr_level: float,
        update_analysis_status_fn: callable,
        p_values_path: str = None,
        activation_path: str = None,
        z_scores_path: str = None,
):
    update_analysis_status_fn(
        user_id=user_id,
        analysis_id=analysis_id,
        status="Running BH FDR correction",
        pvalues_path=p_values_path,
        activation_path=activation_path,
    )
    print("Running Benjamini-Hochberg FDR correction...")
    _, reject, p_val_thresholds = bh_fdr_correction(p_value_df=p_values_df, alpha=fdr_level)
    p_val_thresholds = p_val_thresholds.dropna()
    p_val_thresholds = pd.DataFrame({ "p_thresholds": p_val_thresholds, "fdr": fdr_level })
    tf_counts = reject.sum().sort_values(ascending=False)
    tf_counts = tf_counts[tf_counts > 0]  # Keep only TFs with at least one rejection

    final_output = pd.DataFrame(0, index=reject.index, columns=reject.columns)
    # Where the rejection is True, fill in the direction from activation_df
    final_output[reject] = activation_df[reject]

    # Propagate NaNs for cells/TFs where analysis couldn't be run
    final_output[p_values_df.isna()] = np.nan
    final_output = final_output.astype("Int64")  # Use a nullable integer type

    # bh_reject_path = os.path.join(result_path, "bh_reject.csv")
    bh_reject_path = os.path.join(result_path, "bh_reject.parquet")
    print(f"Saving FDR results to {bh_reject_path}...")
    # final_output.to_csv(bh_reject_path, index=True)
    final_output.to_parquet(bh_reject_path)

    p_val_threshold_path = os.path.join(result_path, "p_val_thresholds.csv")
    print(f"Saving p-value thresholds to {p_val_threshold_path}")
    p_val_thresholds.to_csv(p_val_threshold_path, index=True, index_label="TF")

    update_analysis_status_fn(
        user_id=user_id,
        analysis_id=analysis_id,
        status="Completed",
        tfs=tf_counts.index.tolist() if (tf_counts.index.tolist()) else [],
        bh_reject_path=bh_reject_path,
        fdr_level=fdr_level,
        p_val_threshold_path=p_val_threshold_path,
        z_scores_path=z_scores_path
    )