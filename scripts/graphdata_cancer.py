import pickle
from src.data_paths import *
import pandas as pd
import numpy as np
import os
import igraph as ig


### creates graph_data.pkl file which contains the nodes and edges of the graph
class Graphdata_Cancer_Creator:
    def __init__(self, subfolder: str):
        self.GRAPH_DATA_PATH = 'graph_data/' + subfolder + '/'
        self.file_name = "graph_data"
        self.df = pd.read_pickle(MERGED)
        self.cosmic_protein = pd.read_csv(COSMIC_PROTEINS, sep='\t')
        self.all_prots_related_to_cancer_drugs = pd.DataFrame()
        self.direct_effect_prots = pd.DataFrame()
    
    ## Create the whole procudeure from preprocessing to saving the igraph object
    def create(self):
        print("Data preprocessing")
        self.data_preprocessing()
        print("Creating graph data")
        self.create_graph_data()
        print("Saving igraph object")
        self.save_igraph()

    def data_preprocessing(self):
        bind_cancer_proteins_by_name = self.df[self.df['drugbank_protein_name'].isin(self.cosmic_protein['Gene'])]
        bind_cancer_proteins_by_uniprot = self.df[self.df['swissprot_protein_id'].isin(self.cosmic_protein['Uniprot'])]

        ## top most studies proteins
        most_studied_proteins = ['Cytochrome P450 3A4', 'Epidermal growth factor receptor', 'Proto-oncogene tyrosine-protein kinase Src', 'Vascular endothelial growth factor receptor 2', 'Adenosine receptor A2a', 'Cytochrome P450 2C9', 'Cytochrome P450 1A2', 'Cytochrome P450 2C19', 'Cytochrome P450 2D6', 'Prostaglandin G/H synthase 1', 'Prostaglandin G/H synthase 2']

        # Add most studied proteins
        bind_cancer_proteins_by_name = pd.concat([bind_cancer_proteins_by_name, self.df[self.df['target_name'].isin(most_studied_proteins)]])

        self.cosmic_protein[self.cosmic_protein['Gene'] == 'EGFR'][['Gene', 'Gene synonym', 'Uniprot']]

        cancer_related_proteins = pd.concat([bind_cancer_proteins_by_name, bind_cancer_proteins_by_uniprot]).drop_duplicates()
        cancer_related_proteins_df = self.df[self.df['target_name'].isin(cancer_related_proteins['target_name'])]

        ligands_related_to_cancer_proteins = cancer_related_proteins_df.dropna(subset=['ligand_name'])

        drugs_related_to_cancer_proteins = ligands_related_to_cancer_proteins.dropna(subset='drugbank_drug_class_superclass')
        all_prots_related_to_cancer_drugs = self.df[self.df['drugbank_drug_name'].isin(drugs_related_to_cancer_proteins['drugbank_drug_name'].unique())]

        direct_prots_related_target_names = cancer_related_proteins['target_name'].unique()

        direct_effect_prots = all_prots_related_to_cancer_drugs[all_prots_related_to_cancer_drugs['target_name'].isin(direct_prots_related_target_names)]

        self.direct_effect_prots = direct_effect_prots
        self.all_prots_related_to_cancer_drugs = all_prots_related_to_cancer_drugs

    def _get_target_name_nonmutant(self, x):
        return x.split('[')[0]
    
    def create_graph_data(self):
        
        self.df['target_name_nonmutant'] = self.df['target_name'].apply(self._get_target_name_nonmutant)
        self.df['target_name'].nunique(), self.df['target_name_nonmutant'].nunique()

        unique_prots = self.df['target_name_nonmutant'].unique()

        protein_to_doi = {prot: set() for prot in unique_prots}
        for _, row in self.df.iterrows():
            protein = row['target_name_nonmutant']
            doi = row['doi']
            protein_to_doi[protein].add(doi)
        
        # Create a list of nodes
        nodes = []
        node_id_map = {}  # To map node names to unique ids
        categories = [{'name': 'Possible Cancer Drug'}, {'name': 'Cancer Protein'}]

        # Process ligand nodes
        ligand_nodes = self.all_prots_related_to_cancer_drugs['drugbank_drug_name'].dropna().unique()
        for idx, ligand in enumerate(ligand_nodes):
            ligand = str(ligand)
            node_id_map[ligand] = idx
            node = {
                'id': str(idx),
                'Label': ligand,
                'category': 0,  # Index of 'ligand' in categories
            }
            nodes.append(node)

        # Process protein nodes
        cancer_protein_nodes = self.direct_effect_prots['target_name'].unique()

        # Start from the length of the ligand nodes (to avoid overlapping ids)
        for idx, protein in enumerate(cancer_protein_nodes, start=len(node_id_map)):
            node_id_map[protein] = idx
            node = {
                'id': str(idx),
                'Label': protein,
                'category': 1,  # Index of 'protein' in categories
            }
            nodes.append(node)

        # Create a list of edges
        edges = []
        for _, row in self.direct_effect_prots.iterrows():
            source_id = node_id_map[row['drugbank_drug_name']]
            target_id = node_id_map[row['target_name']]
            strength = row['ec50'] if pd.isna(row['ic50']) else row['ic50']  # Use ec50 if ic50 is NaN
            
            if pd.isna(strength) or strength < 1e-8 or np.log(strength) < 1e-8:
                continue

            strengths = np.log(np.log(strength))  # Log-transform the value to make it more visually appealing

            edge = {
                'source': str(source_id),
                'target': str(target_id),
                'sourceLabel': row['drugbank_drug_name'],
                'targetLabel': row['target_name'],
                'Weight': strengths  
            }
            edges.append(edge)

        # Modify strength to be in the range [0, 1]
        max_strength = max([edge['Weight'] for edge in edges])
        min_strength = min([edge['Weight'] for edge in edges])

        for i, edge in enumerate(edges):
            edges[i]['Weight'] = (edge['Weight'] - min_strength) / (max_strength - min_strength)
            edges[i]['Weight'] = np.clip(edges[i]['Weight'], 0.4, 1)
            edges[i]['Weight'] = (edges[i]['Weight'] - 0.4) / 0.6

        # Create the graph data
        graph_data = {
            'nodes': nodes,
            'links': edges,
            'categories': categories
        }
        
        self.graph_data = graph_data
    
    def save_igraph(self):

        # Extract node and edge data
        nodes = self.graph_data['nodes']
        edges = self.graph_data['links']
        categories = self.graph_data['categories']

        # Create an igraph graph
        g_cancer = ig.Graph(directed=False)  # Set directed=True if your data is directional
        g_cancer.add_vertices(len(nodes))

        # Add vertex attributes
        for i, node in enumerate(nodes):
            g_cancer.vs[i]['name'] = node['Label']
            g_cancer.vs[i]['category'] = categories[node['category']]['name']

        # Add edges
        edge_list = [(int(e['source']), int(e['target'])) for e in edges]
        g_cancer.add_edges(edge_list)

        # Add edge attributes
        if 'Weight' in edges[0]:
            for i, e in enumerate(edges):
                g_cancer.es[i]['affinity'] = e['Weight']
                g_cancer.es[i]['sourceLabel'] = e['sourceLabel']
                g_cancer.es[i]['targetLabel'] = e['targetLabel']

        for v in g_cancer.vs:
            # If this vertex's category is that of a ligand, set type=True
            v['type'] = (v['category'] == 'Possible Cancer Drug')

        ### Save igraph and edges/nodes to pickle to load into gephi
        if not os.path.exists(self.GRAPH_DATA_PATH):
            os.makedirs(self.GRAPH_DATA_PATH)
        with open(self.GRAPH_DATA_PATH +'igraph_object.pkl', 'wb') as f:
            pickle.dump(g_cancer, f)
