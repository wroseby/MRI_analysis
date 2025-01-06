# FUNCTIONS FOR GRAPH THEORETICAL ANALYSIS OF CORRELATION NETWORKS #

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from graph_tool.all import *
import matplotlib.cm as cm
import scipy

# Test function
def test(a,b):
    c = a + b
    return c

# Function to generate a matrix of absolute partial correlations given wide-form data
# E.g. to get partial correlations of surface areas given parcel surface areas
# requires pandas
def partial_cor(data):
    data_cov = data.cov().values # covariance values
    precision_matrix = np.linalg.inv(data_cov) # invert cov matrix to get precision matrix
    data_pcorr = np.zeros_like(precision_matrix) # initialise zero matrix
    for i in range(precision_matrix.shape[0]): # fill zero matrix
        for j in range(precision_matrix.shape[1]):
            if i == j: # diagonal
                data_pcorr[i, j] = 1.0 # 1s on the diagonal
            else: # off-diagonal
                data_pcorr[i, j] = -precision_matrix[i, j] / np.sqrt(precision_matrix[i, i] * precision_matrix[j, j]) # partial correlation

    data_pcorr = pd.DataFrame(data_pcorr, index=data.columns, columns=data.columns) # turn it into a pandas DataFrame

    return (abs(data_pcorr)) # return the absolute partial correlations

# Function to convert given data to a list of two graphs
# Requires a csv file with mixed data for two groups, one 'Control' and one experimental given by groupname
def data_to_graphs(filename, groupname):
    data = pd.read_csv(filename)  # read data

    pcormats = [partial_cor(data[data['Group'] == 'Control'].iloc[:, 7:]),
                # generate partial correlations from data
                partial_cor(data[data['Group'] == groupname].iloc[:, 7:])]  # and put into list

    graphs = []  # initialise list of graphs
    for pcormat in pcormats:  # for each partial correlation matrix
        graphs.append(Graph(scipy.sparse.lil_matrix(pcormat), directed=False))  # generate graphs
        remove_self_loops(graphs[-1])  # remove self-correlations

    return graphs

# Function to plot distributions of strengths
# Requires a list of two graphs and will plot the distributions for each graph
def plot_strengths(graphs):
    strengths = []
    for graph in graphs: # for each graph
        weight = graph.edge_properties["weight"] # get edge weights
        strengths.append(graph.degree_property_map("total", weight = weight).a) # get total weight for each vertex

    fig = plt.figure(figsize=(10, 6)) # initialise figure
    ax = fig.add_subplot(1, 1, 1) # add a single panel
    counts1, bins1 = np.histogram(strengths[0], bins=30, density=True)  # get the histogram of control strengths
    plt.plot(bins1[:-1], counts1, label='Control') # plot histogram
    counts2, bins2 = np.histogram(strengths[1], bins=30, density=True)  # get the histogram of synesthete strengths
    plt.plot(bins2[:-1], counts2, label='Syn') # plot histogram
    fig.set_facecolor('white') # make background white
    ax.set_facecolor('white') # make background white
    plt.xlabel("Node Strength") # x-axis title
    plt.ylabel("Density") # y-axis title
    plt.legend() # show the legend
    plt.title("Node Strength Distributions") # set the title
    plt.show() # show plot

# Function to calculate and agglomerate network metrics
# currently includes: clustering, efficiency, characteristic path length, mean edge and vertex betweeness; also a histogram of shortest distances
# requires pandas, graph-tool, numpy
def measure_net(graph):
    net_metrics = pd.DataFrame() # initialise dataframe to store results

    clustering = global_clustering(graph, weight=graph.ep.weight)[0]  # global weighted clustering

    weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
    for e in graph.edges(): # for all the graph edges
        weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
    graph.ep['weight_inv'] = weight_inv # set inverted weights as graph edge property
    hist = distance_histogram(graph, weight = graph.ep.weight_inv, bins=[0,0.1]) # get distance histogram with inverted weights
    total_efficiency = 0.0  # initialise total efficiency for graph
    count = 0  # initialise a count
    distances = shortest_distance(graph, weights=graph.ep.weight_inv) # shortest distances between vertices using inverted weights
    distances_sum = 0 # initialise the sum of distances
    for j in range(graph.num_vertices()):  # for all vertices
        for k in range(j + 1, graph.num_vertices()):  # to all other vertices
            if distances[j][k] != 0:  # if shortest distance is not 0
                total_efficiency += 1.0 / distances[j][k]  # add inverse of shortest distance to total efficiency
                count += 1  # and add one count
                distances_sum += distances[j][k]
    L_obs = distances_sum / ((graph.num_vertices() - 1) * graph.num_vertices()) # characteristic (observed) path length
    global_efficiency = total_efficiency / count  # global efficiency of graph is mean of efficiency

    vb, eb = betweenness(graph, weight=graph.ep.weight_inv)  # calculate edge and vertex betweenness

    net_metrics['clustering'] = clustering
    net_metrics['efficiency'] = global_efficiency
    net_metrics['L_obs'] = L_obs
    net_metrics['mean_eb'] = np.mean(eb)
    net_metrics['mean_vb'] = np.mean(vb)

    return(net_metrics, hist)

# Function to plot graphs in the following style:
# anatomical layout, given x and y coordinates
# nodes coloured by according to group membership, given group numbers (in x and y data)
# edges brighter, thicker and closer to top of layer with increasing weight
# requires graph-tool, matplotlib.cm
def draw_graph_anatomical(graph, positions, colours):
    def hex_to_rgb(hex_color): # helper function to convert hex colour to RGB tuple
        hex_color = hex_color.lstrip('#')
        return tuple(int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))

    # set vertex positions and colours:
    pos = graph.new_vertex_property("vector<double>") # create new vertx property for position
    color = graph.new_vertex_property("vector<double>") # create new
    for i, vertex in enumerate(graph.vertices()): # for each vertex
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y']) # get x and y positions
        color[vertex] = hex_to_rgb(colours[positions.loc[i, 'group']]) # set colour by group

    # set edge colours:
    edge_color= graph.new_edge_property("vector<double>") # initialise new edge property to store edge colours
    for e in graph.edges(): # for each edge
        weight_val = graph.ep.weight[e] # acquire edge weight
        rgba = list(cm.viridis(weight_val)) # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2) # set alpha values, scaled
        edge_color[e] = rgba # set edge colour
    graph.ep["color"] = edge_color # set edge colours as edge property

    graph_draw( # draw the graph
        graph,
        pos=pos, # set vertex positions
        vertex_fill_color=color, # set vertex fill colours
        edge_color=graph.ep.color, # refer to edge colours
        edge_pen_width=prop_to_size(graph.ep.weight, 0.1, 3, power = 2, log = False), # set edge widths based on weights
        eorder=graph.ep.weight, # set edge order (largest on top)
        vertex_color = 'black', # set vertex stroke colour
        vertex_size = 8, # set vertex size
        output_size = (800,800) # set vertex colour
    )

# Function to merge correlation weights with physical distances and plot scatterplot, given correlation matrix and physical distances
def weight_by_dist(cormat, dists):
    flatmat = cormat.reset_index().melt(id_vars='index', var_name='area_1',value_name='weight')  # convert correlation matrix to long form

    graph_dist = pd.merge(flatmat, dists, left_on=['index', 'area_1'], right_on=['area_2', 'area_1'], how='left')  # merge distance and weight data
    graph_dist['pair'] = graph_dist.apply(lambda row: tuple(sorted([row['index'], row['area_1']])), axis=1) # recreate 'pair' column by reordering areas
    graph_dist = graph_dist[graph_dist['index'] != graph_dist['area_1']]  # drop self-loops
    edge_dist = graph_dist.drop_duplicates(subset='pair').reset_index(drop=True)  # drop repeated pairs

    fig, axs = plt.subplots(1, 1, figsize=(8, 8))
    axs.scatter(edge_dist['distance'], edge_dist['weight'], alpha=0.8, s=1.5, c=edge_dist['weight'])  # plot of weight over distance
    axs.set_xlabel('Distance')
    axs.set_ylabel('pcor')

    return graph_dist

# Function to calculate SWP and delta for a graph, given the correlation matrix and distances, and observed clustering and path length of graph in question
# creates a distance-based lattice graph and a single random rewire for the null models
def SWP(cormat, dists, C_obs, L_obs):
    def calculate_path_length(graph): # helper function to get characteristic path length
        weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
        for e in graph.edges():
            weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
        graph.ep['weight_inv'] = weight_inv
        distances = shortest_distance(graph, weights=graph.ep.weight_inv)
        distances_sum = 0
        for j in range(graph.num_vertices()):  # for all vertices
            for k in range(j + 1, graph.num_vertices()):  # to all other vertices
                if distances[j][k] != 0:  # if shortest distance is not 0
                    distances_sum += distances[j][k]
        L = distances_sum / ((graph.num_vertices() - 1) * graph.num_vertices())

        return L
    # generate lattice graph #

    # sort values so shortest distances have highest weights
    distances_sort = dists.sort_values(by=['distance'], ascending=[True])[['area_1', 'area_2', 'distance']].reset_index(drop=True) # sort distances
    lattice_weights =  dists['weight'].sort_values(ascending = False).reset_index(drop=True) # get sorted weights
    distances_sort['weight'] = lattice_weights # add sorted weights to sorted distances

    # change distances object back to covariance matrix and fill with sorted weights
    areas = pd.unique(dists[['area_1', 'area_2']].values.ravel('K'))
    cormat_latt = pd.DataFrame(index=areas, columns=areas).fillna(1)
    for _, row in distances_sort.iterrows():
        cormat_latt.at[row['area_1'], row['area_2']] = row['weight']
        cormat_latt.at[row['area_2'], row['area_1']] = row['weight']

        cormat_latt.reindex(index=cormat[0].index, columns=cormat[0].columns) # reorder matrix to match original order

    graph_latt = Graph(scipy.sparse.lil_matrix(cormat_latt), directed=False) # create lattice graph
    remove_self_loops(graph_latt) # remove self-loops

    # generate random graph  #

    non_diag = [cormat.iloc[j, k] for j in range(cormat.shape[0]) for k in range(cormat.shape[1]) if j != k]
    np.random.shuffle(non_diag)
    cormat_rand = cormat.copy()
    v = 0
    for j in range(cormat.shape[0]):
        for k in range(cormat.shape[1]):
            if j != k:
                cormat_rand.iloc[j, k] = non_diag[v]
                v += 1
    cormat_rand = cormat_rand.copy()
    graph_rand = Graph(scipy.sparse.lil_matrix(cormat_rand), directed=False)

    # calculate clustering coefficients #

    C_latt = global_clustering(graph_latt, weight = graph_latt.ep.weight)[0] # standard global clustering based on graphs from partial correlations
    C_rand = global_clustering(graph_rand, weight = graph_rand.ep.weight)[0] # standard global clustering based on graphs from partial correlations


    # calculate characteristic path length

    L_latt = calculate_path_length(graph_latt)
    L_rand = calculate_path_length(graph_rand)

    print(C_latt, C_rand, L_latt, L_rand) # sanity check - makes sense

    # calculate divergences from lattice and random models
    delta_C = (C_latt - C_obs) / (C_latt - C_rand)
    delta_L = (L_obs - L_rand) / (L_latt - L_rand)

    SWP = np.sqrt((delta_C**2 + delta_L**2) / 2) # calculate SWP
    alpha = np.arctan(delta_L / delta_C) # get angle of divergence vector
    delta = (((4 * alpha) / np.pi) - 1) # get delta value
    delta =  np.clip(delta, -1, 1) # set range

    return SWP, delta