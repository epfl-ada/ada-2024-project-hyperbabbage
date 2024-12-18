import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 




from scipy.stats import norm
import numpy as np

def bootstrap(data, stat_func=np.mean, num_resamples=100, alpha=0.05):
    """
    Parameters:
    - data (array-like): Input data for bootstrapping.
    - stat_func (function): Statistic function to compute (default: np.mean).
    - num_resamples (int): Number of bootstrap resamples (default: 100).
    - alpha (float): Confidence level (default: 0.05).

    Returns:
    - bootstraps (array): The bootstrapped estimate of the statistic.
    """
    # Fix the random seed for reproducibility
    np.random.seed(4)

    data = np.array(data)
    n = len(data)

    # Fit data to a normal distribution
    mu, std = norm.fit(data)

    # Generate probabilities from the PDF
    pdf_values = norm.pdf(data, loc=mu, scale=std)
    probabilities = pdf_values / pdf_values.sum()

    # Weighted sampling based on probabilities
    bootstrap_samples = np.array([
        np.random.choice(
            data, 
            size=n, 
            replace=True, 
            #p=probabilities
            )
        for _ in range(num_resamples)
    ])
    
    # Apply the statistical function to each resample
    bootstraps = np.apply_along_axis(stat_func, axis=1, arr=bootstrap_samples)

    return bootstraps





def bootstrap_chemical_param(df, column):
    """
    Parameters: 
    - filtered_df (pandas dataframe) : cancer-related protein df containing all relevant information
    - column (string) : the colu,mn on which boostrapping must be performed

    Returns: 
    - Bootsrapped array for plotting
    """
    filtered_df = df.copy()
    # Step 2: Filter rows where `ec50` is not NaN
    filtered_df = filtered_df[filtered_df[column].notna()]

    # Initialize an empty DataFrame to store bootstrap results
    bootstrap_df = pd.DataFrame(columns=["target_name", "drugbank_drug_name_present", column])

    # distinguish between drug and non-drug ligands
    for target in filtered_df['target_name'].unique():

        group_with = filtered_df[
            (filtered_df['target_name'] == target) & 
            (filtered_df['drugbank_drug_name_present'] == 'With DrugBank Drug')
        ]
        
        group_without = filtered_df[
            (filtered_df['target_name'] == target) & 
            (filtered_df['drugbank_drug_name_present'] == 'Without DrugBank Drug')
        ]
        
        # Bootstrapping
        bootstrap_group_with = bootstrap(group_with[column].values)
        bootstrap_group_without = bootstrap(group_without[column].values)

        # Convert the bootstrap results to DataFrame
        bootstrap_group_with_df = pd.DataFrame(bootstrap_group_with.T, columns=[column])  # Transpose to get shape (1000, 1)
        bootstrap_group_with_df["target_name"] = target
        bootstrap_group_with_df["drugbank_drug_name_present"] = 'With DrugBank Drug'

        bootstrap_group_without_df = pd.DataFrame(bootstrap_group_without.T, columns=[column])  # Transpose to get shape (1000, 1)
        bootstrap_group_without_df["target_name"] = target
        bootstrap_group_without_df["drugbank_drug_name_present"] = 'Without DrugBank Drug'

        # Concatenate the results into the main DataFrame
        bootstrap_df = pd.concat([bootstrap_df, bootstrap_group_with_df, bootstrap_group_without_df], axis=0)

    # Reset the index of the final DataFrame to avoid index duplication
    bootstrap_df.reset_index(drop=True, inplace=True)

    return bootstrap_df
