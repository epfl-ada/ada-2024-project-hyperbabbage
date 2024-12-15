import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Define a basic bootstrap function
def bootstrap(data, n_iterations=1000, sample_size=None):
    if sample_size is None:
        sample_size = len(data)
    resamples = []
    for _ in range(n_iterations):
        sample = np.random.choice(data, size=sample_size, replace=True)
        resamples.append(sample)
    return np.array(resamples)

# Initialize an empty DataFrame to store bootstrap results
bootstrap_df = pd.DataFrame()

for target in filtered_df['target_name'].unique():
    # Subsets for each group based on 'drugbank_drug_name_present'
    group_with = filtered_df[
        (filtered_df['target_name'] == target) & 
        (filtered_df['drugbank_drug_name_present'] == 'With DrugBank Drug')
    ]
    
    group_without = filtered_df[
        (filtered_df['target_name'] == target) & 
        (filtered_df['drugbank_drug_name_present'] == 'Without DrugBank Drug')
    ]
    
    # Perform bootstrapping
    bootstrap_group_with = bootstrap(group_with[column].values)
    bootstrap_group_without = bootstrap(group_without[column].values)

    # Convert the bootstrap results to DataFrame
    bootstrap_group_with_df = pd.DataFrame(bootstrap_group_with, columns=[column])
    bootstrap_group_with_df["target_name"] = target
    bootstrap_group_with_df["drugbank_drug_name_present"] = 'With DrugBank Drug'

    bootstrap_group_without_df = pd.DataFrame(bootstrap_group_without, columns=[column])
    bootstrap_group_without_df["target_name"] = target
    bootstrap_group_without_df["drugbank_drug_name_present"] = 'Without DrugBank Drug'

    # Concatenate the results into the main DataFrame
    bootstrap_df = pd.concat([bootstrap_df, bootstrap_group_with_df, bootstrap_group_without_df], axis=0)

# Step 4: Plot the boxplot
fig, ax = plt.subplots(figsize=(12, 10))
sns.boxplot(
    data=bootstrap_df, 
    x="target_name",
    y=column,  # Use the column variable for flexibility
    hue="drugbank_drug_name_present"
)

# Customize the plot
plt.yscale("log")
plt.xticks(rotation=90)
plt.title("EC50 values by Protein and DrugBank Drug Presence (Balanced Data)", fontsize=16)
plt.xlabel("DrugBank Protein Name", fontsize=14)
plt.ylabel("EC50 (log scale)", fontsize=14)
plt.legend(title="DrugBank Drug Presence", fontsize=12)
plt.tight_layout()
plt.grid("y")
plt.show()



#####################################################################################################

def bootstraped_chemical_param(filtered_df, column="ec50"):
    # Initialize an empty DataFrame to store bootstrap results
    column = "ec50"
    bootstrap_df = pd.DataFrame(columns=["target_name", "drugbank_drug_name_present", column])

    for target in filtered_df['target_name'].unique():
        # Subsets for each hue
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
        # Make sure the result is a 1D array
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

    # View the result
    return bootstrap_df

###################################################################################################################################

import numpy as np

# Select the top N largest `ec50` values for both groups
column = "ec50"
overall_balanced_data = []
balanced_data = pd.DataFrame()  # Initialize an empty DataFrame
for random_state in range(43):
    for target in filtered_df['target_name'].unique():
        # Subsets for each hue
        group_with = filtered_df[
            (filtered_df['target_name'] == target) & 
            (filtered_df['drugbank_drug_name_present'] == 'With DrugBank Drug')
        ]
        
        group_without = filtered_df[
            (filtered_df['target_name'] == target) & 
            (filtered_df['drugbank_drug_name_present'] == 'Without DrugBank Drug')
        ]
        

        # Find the smaller group size
        min_size = min(len(group_with), len(group_without))
        
        # Randomly sample `min_size` rows from both groups
        sampled_with = group_with.nlargest(n=len(group_with), columns='ec50').sample(n=min_size, random_state=random_state)
        sampled_without = group_without.nlargest(n=len(group_without), columns='ec50').sample(n=min_size, random_state=random_state)


        
        # Concatenate sampled data into the balanced dataset
        balanced_data = pd.concat([balanced_data, sampled_with, sampled_without])
    overall_balanced_data.append(balanced_data)
overall_df = pd.concat(overall_balanced_data, axis=0)

# Step 4: Plot the boxplot
fig, ax = plt.subplots(figsize=(12, 10))
sns.boxplot(
    data=overall_df, 
    x="target_name",
    y="ec50",
    hue="drugbank_drug_name_present"
)

# Customize the plot
plt.yscale("log")
plt.xticks(rotation=90)
plt.title("EC50 values by Protein and DrugBank Drug Presence (Balanced Data)", fontsize=16)
plt.xlabel("DrugBank Protein Name", fontsize=14)
plt.ylabel("EC50 (log scale)", fontsize=14)
plt.legend(title="DrugBank Drug Presence", fontsize=12)
plt.tight_layout()
plt.grid("y")
plt.show()
overall_df

############################### Cancer processes ######################

plt.figure(figsize=(14.5,7))
sns.barplot(x="Biological process", y="ki", hue="Biological process", data=df[df.ki.notna()])
plt.xticks(rotation=90)
plt.yscale("log")
# Adjust legend position (outside the plot)
plt.legend(bbox_to_anchor=(0.0, -0.1), loc='upper left', borderaxespad=0.)
plt.xticks([])
plt.grid(axis="y")
plt.title("mean Ki value for each protein-ligand relationship ")



################################################################################################

# Aggregate the data to count the number of non-null 'ki' values for each biological process
df_counts = (
    df[df.ki.notna()]
    .groupby("Biological process")
    .size()
    .reset_index(name="count")
)

plt.figure(figsize=(14.5, 7))

# Create the bar plot using counts
sns.barplot(
    x="Biological process",
    y="count",
    hue="Biological process",
    data=df_counts
)

plt.xticks(rotation=90)
# Adjust legend position (outside the plot)
plt.legend(bbox_to_anchor=(0.0, -0.1), loc="upper left", borderaxespad=0.)
plt.xticks([])
plt.grid(axis="y")
plt.title("Count of Ki values for each protein-ligand relationship")
plt.ylabel("Count")
plt.xlabel("")
plt.show()