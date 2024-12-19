import pickle
from src.data_paths import *
import pandas as pd
import numpy as np
import os
import igraph as ig


### creates graph_data.pkl file which contains the nodes and edges of the graph
class Graphdata_Gephi_Creator:
    def __init__(self, subfolder: str):
        self.GRAPH_DATA_PATH = 'graph_data/' + subfolder + '/'
        self.file_name = "graph_data"
        self.df = pd.read_pickle(MERGED)
        self.cosmic_protein = pd.read_csv(COSMIC_PROTEINS, sep='\t')
        self.all_prots_related_to_cancer_drugs = pd.DataFrame()
        self.direct_effect_prots = pd.DataFrame()
    
    ## Create the whole procudeure from preprocessing to saving the igraph object
    def create(self):
        self.data_preprocessing()
        self.create_graph_data()
        self.save_csv_for_gephi()

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


        count = {prot: len(dois) for prot, dois in protein_to_doi.items()}
        target_name_to_count = {row['target_name']: count[row['target_name_nonmutant']] for _, row in self.df.iterrows()}


        # Create a list of nodes
        nodes = []
        node_id_map = {}  # To map node names to unique ids
        categories = [{'name': 'Ligand'}, {'name': 'Protein'}]

        # Process ligand nodes
        ligand_nodes = self.df['ligand_name'].dropna().unique()
        for idx, ligand in enumerate(ligand_nodes):
            ligand = str(ligand)
            node_id_map[ligand] = idx
            node = {
                'id': str(idx),
                'Label': ligand,
                'category': 0,  # Index of 'ligand' in categories
                'research_count': np.median(list(target_name_to_count.values()))
            }
            nodes.append(node)

        # Process protein nodes
        protein_nodes = self.df['target_name'].dropna().unique()
        for idx, protein in enumerate(protein_nodes, start=len(node_id_map)):
            node_id_map[protein] = idx
            node = {
                'id': str(idx),
                'Label': protein,
                'category': 1,  # Index of 'protein' in categories
                'research_count': target_name_to_count[protein]
            }
            nodes.append(node)

        # Process edges efficiently
        weights_ic = self.df['ic50'].values
        weights_ec = self.df['ec50'].values
        weights = np.where(pd.isna(weights_ic), weights_ec, weights_ic)

        weights = np.log(np.log(weights))
        max_strength = np.quantile(weights[~np.isnan(weights)], 0.95)
        min_strength = np.quantile(weights[~np.isnan(weights)], 0.05)
        weights = (weights - min_strength) / (max_strength - min_strength)
        weights = np.where(np.isnan(weights), 0, weights)
        weights = np.where(weights < 0.1, pd.NA, weights)

        # Edges is df[['ligand_name', 'traget_name', weights]] where weights is not nan
        # Don't use loop
        edges = self.df[['ligand_name', 'target_name']].copy()
        edges['source'] = edges['ligand_name'].map(node_id_map)
        edges['target'] = edges['target_name'].map(node_id_map)
        edges['sourceLabel'] = edges['ligand_name']
        edges['targetLabel'] = edges['target_name']
        edges['Weight'] = weights
        edges = edges.dropna(subset=['source', 'target', 'Weight'])
        edges['Weight'] = edges['Weight'].clip(0.1, 1)

        edges = edges.to_dict(orient='records')

        # Create the graph data
        graph_data = {
            'nodes': nodes,
            'links': edges,
            'categories': categories
        }
        
        self.graph_data = graph_data
    
    def save_csv_for_gephi(self):

        # From the graph data, extract nodes, edges and categories
        nodes = self.graph_data['nodes']
        edges = self.graph_data['links']
        categories = self.graph_data['categories']

        # Create an igraph graph
        g = ig.Graph(directed=False)  # Set directed=True if your data is directional
        g.add_vertices(len(nodes))

        # Add vertex attributes
        for i, node in enumerate(nodes):
            g.vs[i]['name'] = node['Label']
            g.vs[i]['category'] = categories[node['category']]['name']
            g.vs[i]['research_count'] = node['research_count']
            g.vs[i]['type'] = (node['category'] == 'Ligand')

        # Add edges
        edge_list = [(int(e['source']), int(e['target'])) for e in edges]
        g.add_edges(edge_list)

        ## Add edge attributes to the edges
        if 'Weight' in edges[0]:
            for i, e in enumerate(edges):
                g.es[i]['affinity'] = e['Weight']
                g.es[i]['sourceLabel'] = e['sourceLabel']
                g.es[i]['targetLabel'] = e['targetLabel']


        # # g = g.subgraph_edges(g.es.select(affinity_gt=0), delete_vertices=True)

        ## Recompute nodes and edges
        nodes = []
        edges = []

        for i, v in enumerate(g.vs):
            nodes.append({
                'id': i,
                'Label': v['name'],
                'category': v['category'],
                'type': v['type'],
                'research_count': v['research_count']
            })

        for e in g.es:
            edges.append({
                'source': e.source,
                'target': e.target,
                'Weight': e['affinity'],
                'sourceLabel': e['sourceLabel'],
                'targetLabel': e['targetLabel']
            })
        
        
        all_nodes_df = pd.DataFrame(nodes)
        all_edges_df = pd.DataFrame(edges)

        ### Save igraph and edges/nodes to pickle to load into gephi
        if not os.path.exists(self.GRAPH_DATA_PATH):
            os.makedirs(self.GRAPH_DATA_PATH)
        all_nodes_df.to_csv(self.GRAPH_DATA_PATH + 'all_nodes.csv', index=False)
        all_edges_df.to_csv(self.GRAPH_DATA_PATH + 'all_edges.csv', index=False)
        with open(self.GRAPH_DATA_PATH +'igraph_object.pkl', 'wb') as f:
            pickle.dump(g, f)
