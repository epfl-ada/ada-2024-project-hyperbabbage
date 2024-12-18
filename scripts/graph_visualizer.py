from pyecharts.charts import Graph
from pyecharts import options as opts
import numpy as np
import pickle

class Graph_Visualizer():
    def __init__(self, subfolder:str):
        self.GRAPH_DATA_PATH = 'graph_data/' + subfolder + '/'
        self.file_name = "graph_data"
        self.igraph_object = None
    
    def get_igraph_object(self):
        with open(self.GRAPH_DATA_PATH + 'igraph_object.pkl', 'rb') as f:
            self.igraph_object = pickle.load(f)
        if self.igraph_object is None:
            raise Exception("No igraph object found")
    
    def render_html(self, project_on: str):
        proj_nodes, proj_edges = self._bipartite_projection_with_affinity(self.igraph_object, project_on=project_on)
        graph = self._generate_echarts_graph(proj_nodes, proj_edges)
        graph.set_global_opts(
            title_opts=opts.TitleOpts(title="Cancer Drugs Projected Interaction Network")
        )
        graph.render(self.GRAPH_DATA_PATH +f"cancer_{project_on}_graph.html")

    def visualize(self):
        self.get_igraph_object()
        self.render_html("ligand")
        self.render_html("protein")
    
    def _bipartite_projection_with_affinity(self, g, project_on="ligand"):
        """
        Compute a bipartite projection of the given graph `g` onto either the ligand or protein side.
        Then, aggregate 'affinity' values for each projected edge.

        Parameters:
        g : igraph.Graph
            The original bipartite graph with
            - `g.es['affinity']` set on the bipartite edges
        project_on : str, optional, default: 'ligand'
            Which set of nodes to project on. 
            'ligand' projects onto ligand nodes (type=True),
            'protein' projects onto protein nodes (type=False).

        Returns:
        proj_nodes : list of dict
            A list of node dictionaries with 'id' and 'name'.
        proj_edges : list of dict
            A list of edge dictionaries with 'source', 'target', and aggregated 'affinity'.
        """

        # Perform bipartite projection
        # g_ligand_proj: projection onto ligand nodes
        # g_protein_proj: projection onto protein nodes
        
        g_ligand_proj, g_protein_proj = g.bipartite_projection(multiplicity=True)

        # Choose which projection to work with
        if project_on == "ligand":
            g_proj = g_ligand_proj
        elif project_on == "protein":
            g_proj = g_protein_proj
        else:
            raise ValueError("project_on must be either 'ligand' or 'protein'")

        # Create lookup from node name to original index for convenience
        name_to_index = {v['name']: v.index for v in g.vs}

        # Compute aggregated affinity on projected edges
        for e in g_proj.es:
            source_name = g_proj.vs[e.source]['name']
            target_name = g_proj.vs[e.target]['name']

            # Map back to original graph indices
            source_index = name_to_index[source_name]
            target_index = name_to_index[target_name]

            # Find common neighbors in original graph
            source_neighbors = set(g.neighbors(source_index))
            target_neighbors = set(g.neighbors(target_index))
            shared_intermediates = source_neighbors & target_neighbors

            # We do the average of the affinities
            affinities = []
            for p in shared_intermediates:
                eid_source_p = g.get_eid(source_index, p)
                eid_target_p = g.get_eid(target_index, p)

                aff_source_p = g.es[eid_source_p]['affinity']
                aff_target_p = g.es[eid_target_p]['affinity']

                avg_aff = (aff_source_p + aff_target_p) / 2.0
                affinities.append(avg_aff)

            e['affinity'] = sum(affinities)/len(affinities) if affinities else None

        # Retrieve nodes with their names
        proj_nodes = []
        for v in g_proj.vs:
            proj_nodes.append({
                'id': v.index,
                'name': v['name']
            })

        # Retrieve edges with their aggregated affinity
        proj_edges = []
        for e in g_proj.es:
            proj_edges.append({
                'source': e.source,
                'target': e.target,
                'affinity': e['affinity']
            })

        return proj_nodes, proj_edges

    

    def _generate_echarts_graph(self, nodes, edges):
        """
        Generate an ECharts graph from the given nodes and edges.

        Parameters:
        nodes : list of dict
            A list of node dictionaries with 'id', 'name', and 'symbolSize'.
        edges : list of dict
            A list of edge dictionaries with 'source', 'target', and 'lineStyle'.

        Returns:
        graph : pyecharts.charts.Graph
        """

        # Compute node degrees
        node_degrees = {node["id"]: 0 for node in nodes}
        for edge in edges:
            node_degrees[edge["source"]] += 1
            node_degrees[edge["target"]] += 1

        # Assign symbol_size to each node based on its degree
        for node in nodes:
            degree = node_degrees[node["id"]]
            # For example, make the size proportional to degree * 5, fallback to 10 if degree = 0
            node["symbolSize"] = np.clip(degree, 0.1, 10) if degree > 0 else 0.1

        # Assign edge width based on 'affinity' (if available)
        for edge in edges:
            weight = edge.get("affinity", 0.5) * 3
            edge["lineStyle"] = {"width": weight}  # scale it up if needed

        # Create the graph
        graph = Graph(init_opts=opts.InitOpts(width="100%", height="700px"))
        graph.add(
            series_name="",
            nodes=nodes,
            links=edges,
            layout="force",
            edge_length=[100, 250],
            repulsion=200,
            linestyle_opts=opts.LineStyleOpts(opacity=0.5)  # Base line style

        )

        graph.set_series_opts(
            label_opts=opts.LabelOpts(
                is_show=False,  # Hide labels by default
                position="right",
                formatter="{b}"  # Use node name as the label
            )
        )

        return graph
