import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def bh_fdr_correction(p_value_df: pd.DataFrame, alpha: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    p_value_df.dropna(axis=1, how="all", inplace=True)

    df_adjusted_pvals = pd.DataFrame(index=p_value_df.index, columns=p_value_df.columns, dtype=float)
    df_reject_status = pd.DataFrame(index=p_value_df.index, columns=p_value_df.columns, dtype="boolean")

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


def run_bh_correction_and_save_tfs(
        user_id,
        analysis_id,
        p_values_df: pd.DataFrame,
        update_analysis_status_fn: callable,
        p_values_path: str = None,
        activation_path: str = None,
        z_scores_path: str = None,
):
    fdr_level: float = 0.05
    _, reject, _ = bh_fdr_correction(p_value_df=p_values_df, alpha=fdr_level)
    tf_counts = reject.sum().sort_values(ascending=False)
    tf_counts = tf_counts[tf_counts > 0]

    update_analysis_status_fn(
        user_id=user_id,
        analysis_id=analysis_id,
        status="Completed",
        pvalues_path=p_values_path,
        activation_path=activation_path,
        tfs=tf_counts.index.tolist() if (tf_counts.index.tolist()) else [],
        z_scores_path=z_scores_path
    )