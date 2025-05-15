import os
import pandas as pd
import numpy as np
from scipy.stats import zscore
import requests
from scipy.special import erf
from joblib import Parallel, delayed
import math  # Added for math.sqrt, math.pi

# Path to the "uploads" folder (use absolute path for robustness)
# Assuming __file__ is in a script directory, and "uploads" is in its parent's parent.
# Adjust if your structure is different, e.g., os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "uploads")
# Original: os.path.join(os.path.dirname(os.path.dirname(__file__)), "uploads")
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(os.path.dirname(SCRIPT_DIR),
                          "uploads")  # Assumes script is in a subdir, uploads is one level up from script's parent.
# If uploads is sibling to script's parent dir:
# UPLOAD_DIR = os.path.join(SCRIPT_DIR, "..", "..", "uploads") - This matches original intent more closely if __file__ is in proj/src/script.py and uploads is in proj/uploads
# For robustness, let's ensure the UPLOAD_DIR is created if it doesn't exist
if not os.path.exists(UPLOAD_DIR):
    print(f"Warning: UPLOAD_DIR '{UPLOAD_DIR}' does not exist. Attempting to create.")
    try:
        os.makedirs(UPLOAD_DIR, exist_ok=True)
    except OSError as e:
        print(f"Error creating UPLOAD_DIR '{UPLOAD_DIR}': {e}. File caching might fail.")


# --- Analytical Calculation Functions ---
def calculate_std_dev_mean_norm_rank(n_population: int, k_sample: int) -> float:
    """
    Calculates sigma_m_k: the theoretical standard deviation of the mean of k
    normalized ranks (r_i = (R_i - 0.5)/n) sampled without replacement
    from a population of size n.
    """
    if not isinstance(n_population, int) or n_population <= 0:
        raise ValueError("n_population must be a positive integer.")
    if not isinstance(k_sample, int) or k_sample <= 0:
        raise ValueError("k_sample must be a positive integer.")
    if k_sample > n_population:
        raise ValueError(
            f"Sample size k_sample ({k_sample}) cannot be greater than population size n_population ({n_population}).")

    if n_population == 1:  # Hence k_sample must be 1
        return 0.0  # Mean is always 0.5, std dev is 0

    if k_sample == n_population:  # All items sampled, mean is always 0.5
        return 0.0

    # Var(m_k) = ( (n+1)*(n-k) ) / ( 12 * n^2 * k )
    variance_m_k = ((n_population + 1.0) * (n_population - k_sample)) / \
                   (12.0 * n_population ** 2 * k_sample)

    if variance_m_k < -1e-9:  # Allow for small negative due to float precision near zero
        raise ValueError(
            f"Negative variance ({variance_m_k}) calculated for sigma_m_k (n={n_population}, k={k_sample})")
    if variance_m_k < 0:  # Treat very small negatives as zero
        variance_m_k = 0.0

    return math.sqrt(variance_m_k)


def calculate_analytical_properties_S_k(n_population: int, k_sample: int) -> tuple[float, float]:
    """
    Calculates analytical mean (mu_S_k) and std dev (sigma_S_k) of
    S_k = min(m_k, 1-m_k), assuming m_k ~ N(0.5, sigma_m_k^2).
    m_k is the mean of k normalized ranks sampled from n_population.

    Returns:
        tuple: (mu_S_k, sigma_S_k)
    """
    if k_sample == 0 or k_sample > n_population:  # Should be caught by calculate_std_dev_mean_norm_rank
        # This case indicates an issue upstream or invalid k_sample.
        # Return NaN or raise error, consistent with how it's handled.
        return np.nan, np.nan

    sigma_m_k = calculate_std_dev_mean_norm_rank(n_population, k_sample)

    if sigma_m_k == 0.0:  # Happens when k_sample == n_population (or n_population=1)
        # m_k is exactly 0.5. So S_k = min(0.5, 0.5) = 0.5.
        # E[S_k] = 0.5, StdDev(S_k) = 0.
        return 0.5, 0.0

    # For m_k ~ N(0.5, sigma_m_k^2), let Z = m_k - 0.5. So Z ~ N(0, sigma_m_k^2).
    # S_k = min(m_k, 1-m_k) = 0.5 - |Z|.

    # E[|Z|] for Z ~ N(0, sigma_m_k^2)
    mean_abs_Z = sigma_m_k * math.sqrt(2.0 / math.pi)
    mu_S_k = 0.5 - mean_abs_Z

    # Var(|Z|) = E[|Z|^2] - (E[|Z|])^2 = E[Z^2] - (mean_abs_Z)^2
    # E[Z^2] = Var(Z) + (E[Z])^2 = sigma_m_k^2 + 0 = sigma_m_k^2
    var_abs_Z = sigma_m_k ** 2 - (mean_abs_Z ** 2)  # This is sigma_m_k**2 * (1 - 2/math.pi)

    if var_abs_Z < -1e-9:  # Allow for small negative due to float precision
        raise ValueError(
            f"Negative variance ({var_abs_Z}) for |Z| (n={n_population}, k={k_sample}, sigma_m_k={sigma_m_k})")
    if var_abs_Z < 0:
        var_abs_Z = 0.0

    sigma_S_k = math.sqrt(var_abs_Z)

    return mu_S_k, sigma_S_k


def get_analytical_null_properties(total_genes: int):
    """
    Generates analytical mean (mu_S_k) and standard deviation (sigma_S_k)
    for S_k = min(m_k, 1-m_k), for k from 1 to total_genes.
    Saves and loads these properties to/from a file.
    """
    # Filename reflects that properties are for k up to total_genes.
    props_file_name = f"analytical_null_props_S_k_N{total_genes}.npz"
    props_file_path = os.path.join(UPLOAD_DIR, props_file_name)

    if os.path.isfile(props_file_path):
        print(f"Analytical null properties file exists: {props_file_path}. Reading it.")
        data = np.load(props_file_path)
        return data["means_S_k"], data["std_devs_S_k"]

    print(f"Analytical null properties file does not exist: {props_file_path}. Generating it.")

    # Arrays will store properties for k=1, ..., total_genes.
    # Index i corresponds to k = i+1.
    means_S_k_array = np.zeros(total_genes)
    std_devs_S_k_array = np.zeros(total_genes)

    for k_idx in range(total_genes):  # k_idx from 0 to total_genes-1
        k_sample = k_idx + 1  # k_sample from 1 to total_genes

        try:
            mu_S_k, sigma_S_k = calculate_analytical_properties_S_k(total_genes, k_sample)
        except ValueError as e:
            print(f"Error calculating analytical properties for n={total_genes}, k={k_sample}: {e}")
            # Decide on fallback: NaNs, or re-raise
            mu_S_k, sigma_S_k = np.nan, np.nan  # Fill with NaN if error occurs for a specific k
            # Consider if this scenario should halt execution.

        means_S_k_array[k_idx] = mu_S_k
        std_devs_S_k_array[k_idx] = sigma_S_k

    try:
        np.savez_compressed(file=props_file_path, means_S_k=means_S_k_array, std_devs_S_k=std_devs_S_k_array)
        print(f"Analytical null properties file generated and saved successfully: {props_file_path}")
    except Exception as e:
        print(f"Error saving analytical null properties file '{props_file_path}': {e}")
        # Decide if we should raise this error or just proceed with in-memory values for this run.

    return means_S_k_array, std_devs_S_k_array


# --- Original script functions (modified) ---

# distribution_worker is no longer needed.
# get_sd is replaced by get_analytical_null_properties.

def sample_worker(
        sample: pd.DataFrame,
        prior_network: pd.DataFrame,
        means_S_k_dist: np.array,  # New: array of E[S_k]
        std_devs_S_k_dist: np.array,  # New: array of StdDev(S_k) (replaces `distribution`)
):
    sample.dropna(inplace=True)  # Operates on a copy passed by Parallel
    sample["rank"] = sample.rank(ascending=False)
    sample["rank"] = (sample["rank"] - 0.5) / len(sample)
    sample["rev_rank"] = 1 - sample["rank"]

    # Create copies for result columns to avoid SettingWithCopyWarning if prior_network is a slice
    res_rs = pd.Series(np.nan, index=prior_network.index)
    res_valid_target = pd.Series(np.nan, index=prior_network.index)

    for tf_id, tf_row in prior_network.iterrows():
        targets = tf_row["target"]
        actions = tf_row["action"]

        valid_targets_list = [t for t in targets if t in sample.index]
        valid_target_count = len(valid_targets_list)

        if len(targets) < 3 or valid_target_count < 3:
            # res_rs[tf_id] = np.nan # Already NaN by default
            # res_valid_target[tf_id] = np.nan
            continue

        acti_rs_sum = 0.0  # Sum of ranks for activation
        inhi_rs_sum = 0.0  # Sum of ranks for inhibition (or reverse ranks)

        for i, target_gene in enumerate(targets):
            if target_gene in valid_targets_list:  # Check against pre-filtered list
                action = actions[i]
                target_rank = sample.loc[target_gene, "rank"]  # This could be a Series if gene has multiple rows
                avg_target_rank = np.average(target_rank)  # Ensure scalar if multiple mappings existed for a gene

                if action == 1:  # Upregulates
                    acti_rs_sum += avg_target_rank
                    inhi_rs_sum += (1 - avg_target_rank)  # Using reverse rank for the "opposite" score
                else:  # Downregulates (action == -1)
                    inhi_rs_sum += avg_target_rank
                    acti_rs_sum += (1 - avg_target_rank)

        # Calculate mean rank sums for observed consistency
        mean_acti_rs = acti_rs_sum / valid_target_count
        mean_inhi_rs = inhi_rs_sum / valid_target_count

        # The statistic S_k = min(mean_rank_for_direction1, mean_rank_for_direction2)
        # Here, mean_rank_for_direction1 = mean_acti_rs
        # And mean_rank_for_direction2 (how well it matches inhibition) = mean_inhi_rs
        # This `rs_statistic_S_k` is the observed S_k value.
        rs_statistic_S_k = min(mean_acti_rs, mean_inhi_rs)

        # Determine original sign based on which direction was stronger (smaller mean rank sum)
        # If acti_rs_sum/valid_target_count < inhi_rs_sum/valid_target_count, TF acts as activator.
        # The 'rs' value stored should reflect this S_k, but with a sign.
        # If activation is more consistent (mean_acti_rs is smaller), rs is positive.
        # If inhibition is more consistent (mean_inhi_rs is smaller), rs is negative.
        # This means rs will be rs_statistic_S_k if mean_acti_rs < mean_inhi_rs,
        # and -rs_statistic_S_k if mean_inhi_rs < mean_acti_rs.
        # If they are equal, let's default to positive or handle as a tie.
        if mean_acti_rs <= mean_inhi_rs:  # Includes tie, defaults to positive "activation" trend
            res_rs[tf_id] = rs_statistic_S_k
        else:
            res_rs[tf_id] = -rs_statistic_S_k

        res_valid_target[tf_id] = valid_target_count

    # --- P-value calculation section, vectorized ---
    # Make a temporary DataFrame for calculations to avoid repeated indexing
    # Ensure prior_network is not modified if it's a shared object in Parallel
    # The 'res_rs' and 'res_valid_target' are local to this worker instance's copy of prior_network.
    # For safety, let's create a temporary calculation frame.

    calc_df = pd.DataFrame({
        'rs': res_rs,
        'valid_target': res_valid_target
    }).dropna(subset=['rs'])  # Filter for rows where rs was successfully calculated

    if calc_df.empty:
        # Return NaNs for p-values if no TFs were processed
        return pd.Series(np.nan, index=prior_network.index).values

    observed_S_k_values = np.abs(calc_df['rs'])
    # valid_target indices are 1-based counts, convert to 0-based for array indexing
    indices_for_dist = calc_df['valid_target'].astype(int) - 1

    # Fetch expected mean and std_dev from the precomputed distributions
    expected_means = means_S_k_dist[indices_for_dist]
    expected_std_devs = std_devs_S_k_dist[indices_for_dist]

    # Initialize z_vals Series, aligned with calc_df's index
    z_vals = pd.Series(np.nan, index=calc_df.index)

    # Case 1: std_dev is non-zero
    non_zero_std_mask = (expected_std_devs != 0)
    if np.any(non_zero_std_mask):
        z_vals[non_zero_std_mask] = \
            (observed_S_k_values[non_zero_std_mask] - expected_means[non_zero_std_mask]) / \
            expected_std_devs[non_zero_std_mask]

    # Case 2: std_dev is zero
    # This implies k_sample == n_population. mu_S_k should be 0.5.
    # observed_S_k_values should also be 0.5.
    zero_std_mask = (expected_std_devs == 0)
    if np.any(zero_std_mask):
        # If observed S_k is close to expected mean (0.5), z is 0 (no deviation)
        match_mean_mask = zero_std_mask & (np.isclose(observed_S_k_values, expected_means))
        z_vals[match_mean_mask] = 0.0

        # If observed S_k differs from expected mean (0.5) when std_dev is 0:
        # This is a significant deviation. S_k is always <= 0.5.
        # If observed_S_k < expected_means[zero_std_mask] (i.e. < 0.5), it's more "significant".
        # Assign a large negative Z-score.
        mismatch_mean_mask = zero_std_mask & (~np.isclose(observed_S_k_values, expected_means))
        # Smallest S_k is more significant, leads to more negative Z.
        z_vals[mismatch_mean_mask] = -30.0  # Represents very high significance (p-value ~ 0)

    # Calculate p-values using the Z-scores
    # Original p-value: 1 + erf(z / sqrt(2)) which is 2 * Phi(z) for standard normal CDF Phi.
    # This is suitable if z represents deviation in one direction and smaller S_k is more significant.
    p_vals_calculated = 1 + erf(z_vals / np.sqrt(2))

    # Adjust sign of p-value based on the sign of original 'rs'
    # (to carry directionality: positive p-val for activation, negative for inhibition)
    p_vals_signed = np.where(calc_df['rs'] > 0, p_vals_calculated, -p_vals_calculated)

    # Create the final p-value series to return, aligned with original prior_network index
    final_p_values = pd.Series(np.nan, index=prior_network.index)
    final_p_values[calc_df.index] = p_vals_signed

    return final_p_values.values


def run_analysis(
        tfs,
        gene_exp: pd.DataFrame,
        prior_network: pd.DataFrame,  # This is the grouped one
        means_s_k_for_dist: np.array,  # New
        std_devs_s_k_for_dist: np.array,  # New (replaces `distribution`)
) -> pd.DataFrame:
    gene_exp = gene_exp.T  # Samples as rows, genes as columns

    # Ensure prior_network passed to delayed function is a full copy if it's modified,
    # or ensure modifications inside sample_worker don't affect other parallel calls.
    # Pandas DataFrames are generally copy-on-write in this context, but deepcopy can be safer if modifications are complex.
    # The prior_network itself is not modified for its structure, only for assigning results,
    # but sample_worker now operates on its own series for `rs` and `valid_target` before p-value calc.

    parallel = Parallel(n_jobs=-1, verbose=5, backend="multiprocessing")
    output = parallel(
        delayed(sample_worker)(
            pd.DataFrame(row),  # sample data for one experiment
            prior_network.copy(),  # Pass a copy of prior_network to each worker
            means_s_k_for_dist,
            std_devs_s_k_for_dist
        )
        for idx, row in gene_exp.iterrows()  # Iterate over samples
    )
    output_df = pd.DataFrame(output, columns=tfs, index=gene_exp.index)
    return output_df


def main(prior_network: pd.DataFrame, gene_exp: pd.DataFrame):  # iters is no longer used by core logic
    gene_exp = gene_exp.apply(zscore, axis=1, nan_policy="omit")

    # Grouping prior_network network
    # This creates a new DataFrame, so it's safe.
    grouped_prior_network = prior_network.groupby("tf").agg({"action": list, "target": list})
    # grouped_prior_network["updown"] = grouped_prior_network["target"].apply(lambda x: len(x)) # Not used if analytical uses total_genes

    # `total_genes` is the number of genes available in the expression data for ranking.
    num_total_genes = len(gene_exp)  # gene_exp has genes as rows at this point

    # `max_target` was used to size the old distribution.
    # The new analytical distribution properties will be calculated for k=1 up to num_total_genes.
    # max_tf_targets = grouped_prior_network["updown"].max() # Max targets a TF has in prior
    # We need properties for k up to num_total_genes, as valid_target_count in sample_worker
    # can go up to num_total_genes.

    means_S_k_dist, std_devs_S_k_dist = get_analytical_null_properties(
        total_genes=num_total_genes
    )

    print("Running analysis in multiple cores...")
    return run_analysis(
        tfs=grouped_prior_network.index,
        gene_exp=gene_exp,  # Genes as rows, samples as columns
        prior_network=grouped_prior_network,  # Grouped by TF
        means_s_k_for_dist=means_S_k_dist,
        std_devs_s_k_for_dist=std_devs_S_k_dist,
    )


# --- Data Reading and Main Entry Point (mostly unchanged, check UPLOAD_DIR usage) ---
def read_mouse_to_human_mapping_file():
    mth_file = "mouse_to_human.tsv"
    mth_file_path = os.path.join(UPLOAD_DIR, mth_file)  # Uses global UPLOAD_DIR

    if os.path.isfile(mth_file_path):
        print("Mouse to human mapping file exists. Let's read")
        return pd.read_csv(mth_file_path, sep="\t")

    print("Mouse to human mapping file does not exist. Let's download it.\n")
    file_url = "https://www.cs.umb.edu/~kisan/data/mouse_to_human.tsv"

    try:
        response = requests.get(file_url)
        response.raise_for_status()  # Raise an exception for HTTP errors
        with open(mth_file_path, "wb") as f:
            f.write(response.content)
        print("Mouse to human mapping file downloaded successfully.")

    except requests.exceptions.RequestException as e:
        print(f"Failed to download the mouse to human map file: {e}")
        # Consider if this should halt the program or proceed without mapping if not essential
        raise Exception(f"Failed to download the mouse to human map file: {e}")  # Re-raise

    return pd.read_csv(mth_file_path, sep="\t")


def read_data(p_file: str, g_file: str, organism: str, upload_dir_param: str):
    # Note: upload_dir_param is passed but UPLOAD_DIR global is used by read_mouse_to_human
    # For consistency, file paths here should also use a consistently defined UPLOAD_DIR.
    # If upload_dir_param is intended to override global, adjust accordingly.
    # Current script uses global UPLOAD_DIR for cached files.
    # Let's assume p_file and g_file are relative to upload_dir_param.

    p_file_path = os.path.join(upload_dir_param, p_file)

    prior_network_df = pd.read_csv(
        p_file_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["tf", "action", "target"],
        converters={
            "action": lambda x: (
                1 if x.lower() == "upregulates-expression"  # Make comparison case-insensitive
                else -1 if x.lower() == "downregulates-expression"
                else None  # Or pd.NA for more explicit missing value handling
            )
        },
    ).dropna(subset=['action'])  # Drop rows where action couldn't be converted
    prior_network_df = prior_network_df.reset_index(drop=True)

    g_file_path = os.path.join(upload_dir_param, g_file)
    gene_exp_df = pd.read_csv(g_file_path, sep="\t", index_col=0)

    if organism.lower() == "mouse":
        mouse_to_human = read_mouse_to_human_mapping_file()  # Uses global UPLOAD_DIR

        # Make gene IDs in gene_exp_df uppercase for consistent mapping
        # gene_exp_df.index = gene_exp_df.index.str.upper() # Optional: if case varies
        # mouse_to_human["Mouse"] = mouse_to_human["Mouse"].str.upper() # Optional

        gene_exp_df = gene_exp_df.rename(
            index=dict(zip(mouse_to_human["Mouse"], mouse_to_human["Human"]))
        )
        # Filter for human gene symbols that might be lists like "[GENE1,GENE2]"
        gene_exp_df = gene_exp_df[gene_exp_df.index.astype(str).str.match(r"^\[.*\]$")]
        gene_exp_df.index = gene_exp_df.index.astype(str).str.strip("[]")

        # Explode handles multiple human genes mapped from one mouse probe
        # Need to handle cases where index becomes empty string after strip if original was "[]"
        gene_exp_df = gene_exp_df[gene_exp_df.index != '']
        gene_exp_df = (
            gene_exp_df.assign(index=gene_exp_df.index.str.split(","))
            .explode("index")
        )
        # After explode, 'index' column contains individual gene names. Trim spaces.
        gene_exp_df['index'] = gene_exp_df['index'].str.strip()
        gene_exp_df = gene_exp_df[gene_exp_df['index'] != '']  # Remove empty gene names again

        # If multiple probes/mouse genes map to the SAME human gene, need to aggregate.
        # E.g., mean, median, or error if not handled.
        # Group by new index and take mean, assuming numeric data.
        gene_exp_df = gene_exp_df.groupby("index").mean(numeric_only=True)  # Added numeric_only=True
        # gene_exp_df = gene_exp_df.set_index("index") # Not needed if groupby("index") makes it index

    # Filter genes with too many NaNs (e.g. mostly zeros)
    # Ensure gene_exp_df.columns are samples. len(gene_exp_df.columns) is num_samples.
    # Thresh is min non-NA values.
    min_expression_threshold = int(len(gene_exp_df.columns) * 0.05)
    gene_exp_df = gene_exp_df.replace(0.0, np.nan).dropna(thresh=min_expression_threshold)

    # Further ensure no all-NaN rows after operations (e.g. if a group was all NaN)
    gene_exp_df.dropna(how='all', axis=0, inplace=True)

    return prior_network_df, gene_exp_df


def get_pvalues(prior_file: str, gene_file: str, organism: str, upload_dir: str) -> pd.DataFrame:
    # `iters` is kept in signature for API compatibility but not used by the analytical method.
    try:
        prior_net, gene_e = read_data(prior_file, gene_file, organism, upload_dir)

        if gene_e.empty:
            raise ValueError("Gene expression data is empty after processing. Cannot run analysis.")
        if prior_net.empty:
            raise ValueError("Prior network data is empty. Cannot run analysis.")

        p_values_df = main(prior_net, gene_e)  # iters is passed but not used by new main logic
        p_values_df.dropna(axis=1, how="all", inplace=True)  # Drop TF columns if all p-values are NaN
        return p_values_df

    except Exception as e:
        # Log the full traceback for better debugging
        import traceback
        print(f"Error during get_pvalues execution: {e}\n{traceback.format_exc()}")
        raise Exception(f"Failed to run the analysis: {e}") from e


# Example of how to run (if script is made executable)
if __name__ == '__main__':
    # This is a placeholder for a proper CLI or test setup.
    # You would need to create dummy prior_file.tsv and gene_file.tsv in a test_uploads directory.
    print("Running TF activity prediction script (analytical null distribution).")
    print(f"Using UPLOAD_DIR: {UPLOAD_DIR}")

    # Create dummy files for testing if they don't exist
    test_upload_dir = os.path.join(SCRIPT_DIR, "test_uploads")  # Place test_uploads next to the script
    os.makedirs(test_upload_dir, exist_ok=True)

    dummy_prior_file = "dummy_prior.tsv"
    dummy_gene_file = "dummy_genes.tsv"

    if not os.path.exists(os.path.join(test_upload_dir, dummy_prior_file)):
        with open(os.path.join(test_upload_dir, dummy_prior_file), "w") as f:
            f.write("TF1\tupregulates-expression\tG1\n")
            f.write("TF1\tupregulates-expression\tG2\n")
            f.write("TF1\tupregulates-expression\tG3\n")
            f.write("TF2\tdownregulates-expression\tG2\n")
            f.write("TF2\tdownregulates-expression\tG3\n")
            f.write("TF2\tdownregulates-expression\tG4\n")
            f.write("TF3\tupregulates-expression\tG5\n")  # TF with few targets
            f.write("TF3\tupregulates-expression\tG6\n")

    if not os.path.exists(os.path.join(test_upload_dir, dummy_gene_file)):
        with open(os.path.join(test_upload_dir, dummy_gene_file), "w") as f:
            f.write("Gene\tSample1\tSample2\tSample3\n")
            f.write("G1\t10.0\t12.5\t11.0\n")
            f.write("G2\t8.0\t7.5\t9.0\n")
            f.write("G3\t15.0\t14.0\t16.0\n")
            f.write("G4\t5.0\t6.5\t4.0\n")
            f.write("G5\t20.0\t22.0\t21.0\n")
            f.write("G6\t1.0\t0.5\t1.2\n")  # Low expression gene

    print(f"To test, run with dummy files from: {test_upload_dir}")

    try:
        # Ensure the UPLOAD_DIR for cached analytical properties and mouse map exists.
        # The test files are read from test_upload_dir, but cached files go to UPLOAD_DIR.
        os.makedirs(UPLOAD_DIR, exist_ok=True)

        pvals = get_pvalues(
            prior_file=dummy_prior_file,
            gene_file=dummy_gene_file,
            organism="human",  # or "mouse" if you have a mouse_to_human.tsv in UPLOAD_DIR
            upload_dir=test_upload_dir
        )
        print("\n--- P-values ---")
        print(pvals)

        # Test with organism mouse (will try to download mapping file if not present)
        # print("\n--- Testing with organism: MOUSE ---")
        # pvals_mouse = get_pvalues(
        #     prior_file=dummy_prior_file, # Assuming dummy prior uses human gene names for this test
        #     gene_file=dummy_gene_file,   # And dummy genes are mouse probes mapped in a dummy map
        #     organism="mouse",
        #     iters=100,
        #     upload_dir=test_upload_dir
        # )
        # print("\n--- P-values (Mouse) ---")
        # print(pvals_mouse)

    except Exception as e:
        print(f"An error occurred during the example run: {e}")
        import traceback

        print(traceback.format_exc())