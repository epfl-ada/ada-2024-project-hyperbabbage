from pyecharts.charts import Graph
from pyecharts import options as opts
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from src.data_paths import *

# Prepare your data as before (nodes and edges)

df = pd.read_pickle(MERGED)
pair_counts = df.groupby(['smiles']).size().reset_index(name='count')

enough_connections = pair_counts[pair_counts['count'] > 1000]

df_filtered = df[df['smiles'].isin(enough_connections['smiles'])]

# Create a list of nodes
nodes = []
node_id_map = {}  # To map node names to unique ids
categories = [{'name': 'ligand'}, {'name': 'protein'}]

# Process ligand nodes
ligand_nodes = df_filtered['smiles'].unique()
for idx, ligand in enumerate(ligand_nodes):
    node_id_map[ligand] = idx
    node = {
        'id': str(idx),
        'name': ligand,
        'category': 0,  # Index of 'ligand' in categories
        'symbolSize': 10,  # Adjust as needed
        'value': 1  # You can set this to degree or other measure
    }
    nodes.append(node)

# Process protein nodes
protein_nodes = df_filtered['target_name'].unique()
for idx, protein in enumerate(protein_nodes, start=len(node_id_map)):
    node_id_map[protein] = idx
    node = {
        'id': str(idx),
        'name': protein,
        'category': 1,  # Index of 'protein' in categories
        'symbolSize': 10,  # Adjust as needed
        'value': 1  # You can set this to degree or other measure
    }
    nodes.append(node)

# Create a list of edges
edges = []
for _, row in df_filtered.iterrows():
    source_id = node_id_map[row['smiles']]
    target_id = node_id_map[row['target_name']]
    edge = {
        'source': str(source_id),
        'target': str(target_id),
        'value': row['ic50']  # Use ic50 as edge value
    }
    edges.append(edge)

# Create the graph data
graph_data = {
    'nodes': nodes,
    'links': edges,
    'categories': categories
}


graph = Graph(init_opts=opts.InitOpts(width="100%", height="700px"))
graph.add(
    "",
    nodes,
    edges,
    categories=categories,
    layout="force",
    edge_length=[50, 200],
    repulsion=100,
    linestyle_opts=opts.LineStyleOpts(width=0.5, opacity=0.7),
)
graph.set_series_opts(
    label_opts=opts.LabelOpts(
        is_show=False,  # Hide labels by default
        position="right",
        formatter="{b}",  # Use node name as the label
        #show=False,  # Ensure labels are not rendered statically
    )
)

graph.set_global_opts(
    title_opts=opts.TitleOpts(title="Ligand-Protein Interaction Network"),
    # legend_opts=opts.LegendOpts(is_show=False),
    # tooltip_opts=opts.TooltipOpts(
    #     trigger="item",
    #     formatter=opts.JsCode("""
    #         function (params) {
    #             if (params.dataType === 'node') {
    #                 var categoryIndex = params.data.category;
    #                 var categoryName = params.data.categoryName || categories[categoryIndex].name;
    #                 return params.data.name + '<br/>' + 'Category: ' + categoryName;
    #             } else if (params.dataType === 'edge') {
    #                 return 'Edge from ' + params.data.source + ' to ' + params.data.target;
    #             }
    #         }
    #     """)
    # )
)

graph.render("graph.html")
