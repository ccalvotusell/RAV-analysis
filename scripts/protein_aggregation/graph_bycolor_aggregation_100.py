import numpy as np
import pandas as pd
import graph_tool.all as gt
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# Seed for reproducibility
gt.seed_rng(47)

#############################
# Parameters
#############################

# Layout parameters
K = 10
gamma = 100
mu = 5.0

# Scaling parameters for nodes and edges
node_scale = 2  # Base size of nodes
edge_scale = 5   # Base thickness of edges
outline_width = 0.2  # Width of the node outline

#############################
# Load Data
#############################

# Load the distance matrix
distance_matrix_file = "combined_distance_matrix.txt"
distance_matrix = np.loadtxt(distance_matrix_file)

# Load updated segment order mapping with diffusion coefficients
segment_mapping_file = "updated_segment_order_mapping.txt"
segment_data = pd.read_csv(segment_mapping_file)

# Cap distances (filter for important interactions)
distance_threshold = 100  # Only show connections below this distance
distance_cap = 200  # Cap distances for visualization scaling
distance_matrix = np.clip(distance_matrix, 0, distance_cap)

# Extract segment IDs and component types
segment_ids = segment_data["SegmentID"].values
component_types = segment_data["ComponentType"].values

# Create a lookup dictionary for component types
segment_dict = dict(zip(segment_ids, component_types))

#############################
# Create the Graph
#############################

# Initialize the graph
g = gt.Graph(directed=False)

# Add vertices for each segment
vprop_type = g.new_vertex_property("string")  # Protein type (spike, mucin, albumin)
vprop_color = g.new_vertex_property("vector<float>")  # Node fill color
vprop_size = g.new_vertex_property("float")  # Node size (scaled by importance)
vprop_size_expanded = g.new_vertex_property("float")  # Node size (scaled by importance)
vprop_outline_color = g.new_vertex_property("vector<float>")  # Outline color
vprop_shape = g.new_vertex_property("string")  # Node shape property
vorder = g.new_vertex_property("int")  # Vertex drawing order

vertices = {}

# Define custom futuristic colors
color_spike = [0.0, 0.4, 1.0, 1.0]  # Electric blue
color_albumin = [0.0, 1.0, 0.5, 1.0]  # Neon green
color_mucin = [1.0, 0.2, 0.2, 1.0]  # Bright crimson

for i, segment_id in enumerate(segment_ids):
    v = g.add_vertex()
    vertices[segment_id] = v
    vprop_type[v] = segment_dict[segment_id]

    # Set node colors based on component type
    if vprop_type[v] == "spike":
        vprop_color[v] = color_spike
        vprop_outline_color[v] = [0.0, 0.0, 0.0, 1.0]  # Black outline
        vprop_shape[v] = "square"  # Shape for spikes
        vorder[v] = 2  # Spikes drawn last (on top)
    elif vprop_type[v] == "mucin":
        vprop_color[v] = color_mucin
        vprop_outline_color[v] = [0.0, 0.0, 0.0, 1.0]  # Black outline
        vprop_shape[v] = "triangle"  # Shape for mucins
        vorder[v] = 1  # Mucins drawn in the middle
    elif vprop_type[v] == "albumin":
        vprop_color[v] = color_albumin
        vprop_outline_color[v] = [0.0, 0.0, 0.0, 1.0]  # Black outline
        vprop_shape[v] = "circle"  # Shape for albumins
        vorder[v] = 0  # Albumins drawn first (at the bottom)

# Add edges based on the distance matrix
eprop_weight = g.new_edge_property("double")  # Edge weight
eprop_weight_edge = g.new_edge_property("double")  # Edge weight
edges = []

for i in range(len(segment_ids)):
    for j in range(i + 1, len(segment_ids)):  # Only upper triangle
        distance = distance_matrix[i, j]
        if distance < distance_threshold:  # Apply filter
            e = g.add_edge(vertices[segment_ids[i]], vertices[segment_ids[j]])
            eprop_weight[e] = distance / 200
            eprop_weight_edge[e] = distance
            edges.append(e)

#############################
# Scale Node Sizes by Importance
#############################

# Calculate degree centrality (number of connections per node)
degree_centrality = g.degree_property_map("total")
max_degree = max(degree_centrality.a)

# Scale node sizes based on degree centrality and node_scale
for v in g.vertices():
    vprop_size[v] = node_scale * (0.8 + degree_centrality[v] / max_degree) * 50
    vprop_size_expanded[v] = node_scale * (0.1 + degree_centrality[v] / max_degree) * 40

#############################
# Visualize the Graph with Sorted Edges
#############################

# Define edge colors and widths based on weight
ecmap = cm.get_cmap("Greys_r")  # Use get_cmap for compatibility with older versions
eprop_color = g.new_edge_property("vector<float>")
eprop_width = g.new_edge_property("float")

# Store edges and their weights for sorting
sorted_edges = []
for e in g.edges():
    weight = eprop_weight_edge[e]
    normalized_weight = weight / distance_cap  # Normalize weight (short distances = 0, long = 1)
    eprop_color[e] = list(ecmap(normalized_weight))[:3] + [1.0]  # Add alpha channel
    eprop_width[e] = edge_scale * (0.5 + 2 * (1 - normalized_weight))  # Scale width (shorter = thicker)
    sorted_edges.append((e, weight))

# Sort edges by weight in descending order (shorter distances last for rendering)
sorted_edges = sorted(sorted_edges, key=lambda x: x[1], reverse=True)

# Create edge order based on sorted edges
eorder = g.new_edge_property("int")
for i, (e, _) in enumerate(sorted_edges):
    eorder[e] = i

# Generate a force-directed layout with specified parameters
pos = gt.sfdp_layout(
    g,
    eweight=eprop_weight,
    K=K,
    gamma=gamma,
    mu=mu
)

# Draw the graph with custom shapes, colors, and ordering
gt.graph_draw(
    g,
    pos=pos,
    vertex_fill_color=vprop_color,
    vertex_size=vprop_size,
    vertex_halo=False,  # Disable node halo
    vertex_outline_color=vprop_outline_color,  # Outline color
    vertex_shape=vprop_shape,  # Use custom shapes for nodes
    edge_color=eprop_color,
    edge_pen_width=eprop_width,
    eorder=eorder,  # Use sorted edge order
    vorder=vorder,  # Ensures spikes are drawn last (on top)
    output_size=(10000, 10000),  # High-resolution output
    output="diffusion_colormap_graph_100.png",
)

# Draw the graph with custom shapes, colors, and ordering
gt.graph_draw(
    g,
    pos=pos,
    vertex_fill_color=vprop_color,
    vertex_size=vprop_size_expanded,
    vertex_halo=False,  # Disable node halo
    vertex_outline_color=vprop_outline_color,  # Outline color
    vertex_shape=vprop_shape,  # Use custom shapes for nodes
    edge_color=eprop_color,
    edge_pen_width=eprop_width,
    eorder=eorder,  # Use sorted edge order
    vorder=vorder,  # Ensures spikes are drawn last (on top)
    output_size=(10000, 10000),  # High-resolution output
    output="diffusion_colormap_graph_100_expanded.png",
)

#############################
# Cluster Analysis
#############################

def analyze_clusters(g, vprop_type):
    comp, _ = gt.label_components(g)
    clusters = {}

    for v in g.vertices():
        cluster_id = comp[v]
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(vprop_type[v])

    cluster_types = {"alone_albumin": 0, "alone_mucin": 0, "alone_spike": 0,
                     "albumin_cluster": 0, "mucin_cluster": 0, "spike_cluster": 0,
                     "albumin-mucin": 0, "albumin-spike": 0, "mucin-spike": 0,
                     "albumin-mucin-spike": 0}

    total_clusters = len(clusters)
    for members in clusters.values():
        unique_types = set(members)
        if len(members) == 1:
            cluster_types[f"alone_{members[0]}"] += 1
        elif unique_types == {"albumin"}:
            cluster_types["albumin_cluster"] += 1
        elif unique_types == {"mucin"}:
            cluster_types["mucin_cluster"] += 1
        elif unique_types == {"spike"}:
            cluster_types["spike_cluster"] += 1
        elif unique_types == {"albumin", "mucin"}:
            cluster_types["albumin-mucin"] += 1
        elif unique_types == {"albumin", "spike"}:
            cluster_types["albumin-spike"] += 1
        elif unique_types == {"mucin", "spike"}:
            cluster_types["mucin-spike"] += 1
        else:
            cluster_types["albumin-mucin-spike"] += 1

    with open("cluster_analysis.txt", "w") as f:
        for key, count in cluster_types.items():
            percentage = (count / total_clusters) * 100
            f.write(f"{key} {count} {percentage:.2f}\n")


analyze_clusters(g, vprop_type)

