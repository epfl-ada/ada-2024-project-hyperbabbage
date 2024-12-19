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
class All_Cancerprotein_Affinities:
    def __init__(self):
        self.merged_df, self.protein_classes = self.load_dataset()
        self.drug_biological_process = self.initialize_filtered_df()

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
        # create lists of the relevant columns to identify cancerous proteins
        cancer_proteins = list(self.protein_classes["Gene"].dropna().values)
        genes_uniprot = list(self.protein_classes["Uniprot"].astype(str).dropna().values)

        ### extract columns that contain cancer_keywords 
        pattern_protein_names = '|'.join(rf"\b{re.escape(term)}\b" for term in cancer_proteins)
        pattern_gene_names =  '|'.join(rf"\b{re.escape(term)}\b" for term in genes_uniprot)
        pattern_highly_studied =  '|'.join(rf"\b{re.escape(term)}\b" for term in highly_studied_proteins)

        ### Apply the patterns to the merged_df
        all_cancer_molec_df = self.merged_df[
            self.merged_df['swissprot_protein_id'].str.contains(pattern_gene_names, case=False, na=False) |
            self.merged_df['drugbank_protein_name'].str.contains(pattern_protein_names, case=False, na=False)|
            self.merged_df['target_name'].str.contains(pattern_highly_studied, case=False, na=False)
        ]
        all_cancer_molec_df.reset_index(inplace=True)

        ### Merge the protein classes to the all_cancer_molec_df and group by the biological process
        df = pd.merge(all_cancer_molec_df, self.protein_classes, left_on ="drugbank_protein_name", right_on="Gene", how="left")
        drug_biological_process = df[df["drugbank_drug_name"].notna()].groupby("Biological process")[["drugbank_drug_name"]].count().reset_index()
        return drug_biological_process
    
    def plot_proportion_of_drugs_per_biological_process(self):
        # Extract sizes and categories
        sizes = self.drug_biological_process["drugbank_drug_name"].values
        categories = self.drug_biological_process["Biological process"].values

        # Calculate total and percentages
        total_size = sum(sizes)
        percentages = (sizes / total_size) * 100

        # Create a new DataFrame with aggregated "Others"
        data = pd.DataFrame({"Biological process": categories, "sizes": sizes, "percentages": percentages})
        others = data[data["percentages"] < 1].sum()  # Sum small categories into "Others"
        data = data[data["percentages"] >= 1]  # Keep large categories
        if not others.empty:
            data = pd.concat([data, pd.DataFrame({"Biological process": ["Others"], "sizes": [others["sizes"]], "percentages": [others["percentages"]]})])

        # Generate a color palette for the updated categories
        categories = data["Biological process"].values
        sizes = data["sizes"].values
        palette = sns.color_palette("inferno", n_colors=len(categories))
        colors = [palette[i] for i in range(len(categories))]
        # Create the pie chart
        plt.figure(figsize=(10, 8))
        plt.pie(
            sizes,
            labels=None,  # Remove default labels
            colors=colors,
            autopct='%1.1f%%',  # Percentage values
            startangle=0,
            wedgeprops={"edgecolor": "black", "linewidth": 0.3},
            pctdistance=1.1,  # Position percentages outside the pie
            textprops={'fontsize': 10}
        )

        # Remove the center to create a donut chart with a black edge
        center_circle = Circle((0, 0), 0.60, color='white', linewidth=1.5)
        plt.gca().add_artist(center_circle)

        # Add custom legend below the chart
        legend_handles = [mpatches.Patch(color=colors[i], label=categories[i]) for i in range(len(categories))]
        plt.legend(
            handles=legend_handles,
            title="Biological Process",
            bbox_to_anchor=(0.5, -0.3),  # Position legend below the chart
            loc='center',
            borderaxespad=0.,
            fontsize=10,
            title_fontsize=12,
            ncol=1  # Arrange legend entries into two columns
        )

        # Add a visually distinct title
        plt.title(
            "Proportion of Drugbank Entries by Biological Process",
            fontsize=16,
            fontweight="bold",
            pad=20
        )

        # Optimize layout for neatness
        #plt.tight_layout()

        # Display the chart
        plt.show()
    
    def get_circle_percentages_categories(self):
        # Extract sizes and categories
        sizes = self.drug_biological_process["drugbank_drug_name"].values
        categories = self.drug_biological_process["Biological process"].values

        # Calculate total and percentages
        total_size = sum(sizes)
        percentages = (sizes / total_size) * 100

        # Create a new DataFrame with aggregated "Others"
        data = pd.DataFrame({"Biological process": categories, "sizes": sizes, "percentages": percentages})
        others = data[data["percentages"] < 1].sum()  # Sum small categories into "Others"
        data = data[data["percentages"] >= 1]  # Keep large categories
        if not others.empty:
            data = pd.concat([data, pd.DataFrame({"Biological process": ["Others"], "sizes": [others["sizes"]], "percentages": [others["percentages"]]})])

        # Generate a color palette for the updated categories
        categories = data["Biological process"].values
        sizes = data["sizes"].values
        percentages = [float(round(i/sum(sizes)*100,1)) for i in sizes]

        return percentages, categories
