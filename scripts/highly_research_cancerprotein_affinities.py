from matplotlib.patches import Circle
from src.data_paths import MERGED, COSMIC_PROTEINS
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
from scripts.utils.Chemical_analysis import bootstrap, bootstrap_chemical_param
import matplotlib.patches as mpatches

#################
#   Class Description:
#################
class Highly_Research_Cancerprotein_Affinities:
    def __init__(self):
        self.merged_df, self.protein_classes = self.load_dataset()
        self.filtered_df = self.initialize_filtered_df()
        self.chemical_param_list = ["ki", "kd", "ec50", "ic50"]
        self.BOXPLOT_MINIMUM_COUNT = 5

    def load_dataset(self):
        merged_df = pd.read_pickle(MERGED)
        protein_classes = pd.read_csv(COSMIC_PROTEINS, sep='\t')

        return merged_df, protein_classes
    
    def initialize_filtered_df(self):
        highly_studied_proteins = [
            'Cytochrome P450 3A4', 
            'Epidermal growth factor receptor', 
            'Proto-oncogene tyrosine-protein kinase Src', 
            'Vascular endothelial growth factor receptor 2', 
            'Adenosine receptor A2a', 'Cytochrome P450 2C9', 
            'Cytochrome P450 1A2', 'Cytochrome P450 2C19', 
            'Cytochrome P450 2D6', 
            'Prostaglandin G/H synthase 1', 
            'Prostaglandin G/H synthase 2',
        ]
        ### extract columns that contain cancer_keywords 
        filtered_df = self.merged_df[
            self.merged_df['target_name'].isin(highly_studied_proteins)
        ]
        filtered_df.reset_index(inplace=True)

        ### Create a new column to distinguish rows with and without "drugbank_drug_name"
        filtered_df.loc[:,'drugbank_drug_name_present'] = filtered_df['drugbank_drug_name'].notna().map({True: 'With DrugBank Drug', False: 'Without DrugBank Drug'})
        return filtered_df
    
    def plot_binding_affinity_corr_mtx(self):
        ### Get only the cancer related proteins
        chemical_param_for_corr = self.filtered_df[self.chemical_param_list]
        chemical_param_for_corr

        ### Get all proteins and calculate the pearsons correlation
        merged_df_chemical_param = self.merged_df[self.chemical_param_list]
        correlation_matrix = merged_df_chemical_param.corr(method='pearson')

        ### Plotting the corr-mtx
        plt.figure(figsize=(15, 15))
        sns.heatmap(correlation_matrix, annot=False, cmap='coolwarm', fmt='.2f', cbar=True, square=True)
        plt.title('Correlation Matrix Heatmap of affinity parameters', size=20)
        plt.show()
    
    def get_colleration_mtx_values(self):
        merged_df_chemical_param = self.merged_df[self.chemical_param_list]
        correlation_matrix = merged_df_chemical_param.corr(method='pearson')
        return correlation_matrix
    
    def plot_binding_measurements_to_drug_presence(self, background_color="white"):
        custom_palette = {
            "Without DrugBank Drug": "royalblue",  ### Blue means Present
            "With DrugBank Drug": "orange",  ### Orange means Absent
        }

        _, axes = plt.subplots(2, 2, figsize=(20, 15))
        axes = axes.flatten()
        for i, col in enumerate(self.chemical_param_list):
            data = self.filtered_df[(self.filtered_df[col].notna())].drop_duplicates()

            counts = data.groupby(["target_name", "drugbank_drug_name_present"]).size().unstack(fill_value=0)
            valid_targets = counts[
                (counts["Without DrugBank Drug"] >= self.BOXPLOT_MINIMUM_COUNT) & (counts["With DrugBank Drug"] >= self.BOXPLOT_MINIMUM_COUNT)
            ].index
            filtered_data = data[data["target_name"].isin(valid_targets)]


            if not filtered_data.empty:

                sns.boxplot(
                    data=filtered_data,
                    x="target_name",
                    y=col,
                    hue="drugbank_drug_name_present",
                    ax=axes[i], 
                    palette=custom_palette,
                )
                
                axes[i].set_yscale("log")
                axes[i].tick_params(axis='x', rotation=30)
                for tick in axes[i].get_xticklabels():
                    tick.set_ha('right')

                axes[i].set_title(f"{col} by Protein and DrugBank Drug Presence", fontsize=16)
                axes[i].set_ylabel(f"{col} (log scale)", fontsize=10)
                axes[i].set_xlabel("")
                
                axes[i].legend(title="", fontsize=10, loc='upper right')
                axes[i].grid(axis="y")
                axes[i].set_facecolor(background_color)

        plt.tight_layout()
        
        plt.show()
    
    def plot_measurement_in_different_conditions(self, col = "ic50"):

        if col not in self.filtered_df.columns:
            raise ValueError(f"{col} is not a valid column name for showing this plot.")
        custom_palette = {
            "Without DrugBank Drug": "royalblue",
            "With DrugBank Drug": "orange",
        }

        temperatures = list(range(19, 30))
        pH = [7.4, 7.5]

        data = self.filtered_df[
            (self.filtered_df[col].notna()) &
            (self.filtered_df["ph"].isin(pH)) &
            (self.filtered_df["temp"].isin(temperatures))
        ].drop_duplicates()

        counts = data.groupby(["target_name", "drugbank_drug_name_present"]).size().unstack(fill_value=0)
        valid_targets = counts[
            (counts["Without DrugBank Drug"] >= self.BOXPLOT_MINIMUM_COUNT) & (counts["With DrugBank Drug"] >= self.BOXPLOT_MINIMUM_COUNT)
        ].index
        filtered_data = data[data["target_name"].isin(valid_targets)]

        if not filtered_data.empty:
            plt.figure(figsize=(8, 8))
            sns.boxplot(
                data=filtered_data,
                x="target_name",
                y=col,
                hue="drugbank_drug_name_present",
                palette=custom_palette
            )
            
            plt.yscale("log")
            
            plt.xticks(rotation=45)
            
            plt.title(f"{col} at {temperatures[0]}-{temperatures[-1]}°C and {pH[0]}-{pH[-1]} pH", fontsize=14)
            plt.ylabel(f"{col} (log scale)", fontsize=12)
            plt.xlabel("", rotation=30)
            plt.legend(title="")
            plt.grid(axis="y")
            
            plt.show()
        else:
            raise ValueError("No valid data available for the plot.")

    
    def plot_measurement_by_protein_drug_presence(self, col = "ic50"):

            if col not in self.filtered_df.columns:
                raise ValueError(f"{col} is not a valid column name for showing this plot.")
            chem_prop_df = bootstrap_chemical_param(self.filtered_df, column=col)

            fig, ax = plt.subplots(figsize=(12, 10))
            sns.boxplot(
                data=chem_prop_df, 
                x="target_name",
                y=col,
                hue="drugbank_drug_name_present"
            )
            plt.yscale("log")
            plt.xticks(rotation=90)
            plt.title(f"{col.upper()} values by Protein and DrugBank Drug Presence (Balanced Data)", fontsize=16)
            plt.xlabel("Protein Name", fontsize=14)
            plt.ylabel(f"{col.upper()} (log scale)", fontsize=14)
            plt.legend(title="DrugBank Drug Presence", fontsize=12)
            plt.tight_layout()
            plt.grid("y")
            plt.show()