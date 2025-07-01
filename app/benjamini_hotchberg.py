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
