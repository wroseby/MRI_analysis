# FUNCTIONS FOR GRAPH THEORETICAL ANALYSIS OF CORRELATION NETWORKS #

# Imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from graph_tool.all import *
import matplotlib.cm as cm
import scipy
import itertools

# Test function
def test(a,b):
    c = a + b
    return c

def partial_cor(data, absolute = True):
    """"
    Calculate partial correlation matrix
    :param data: a pandas DataFrame where columns are variables and rows are subjects
    :param absolute: if True then absolute correlation matrix is returned
    """
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

    if absolute:
        data_pcorr = abs(data_pcorr)

    return (data_pcorr) # return the absolute partial correlations

def data_to_graphs(data, groupname = False, log = False, split_sign = False, absolute = True):
    """
    Convert data file into a graph
    :param data: name of file where data is stored
    :param groupname: if True then data is split by the group name and 'Control'
    :param log: if True then data is log-transformed
    :param split_sign: if True then partial correlation matrix is split into positive and negative
    """
    if isinstance(data, str):
        data = pd.read_csv(data)  # if giving filename, read it in

    if log:
        data.iloc[:,7:] = np.log(data.iloc[:,7:])

    if groupname: # if segregating by group
        data = [data[data['Group'] == 'Control'],
                data[data['Group'] == groupname]]

    else: # if not segregating by group
        data = [data] # put into list

    pcormats = [] # initialise partial correlation matrices

    if split_sign: # if splitting positive and negative
        for datum in data:
            pcormat = partial_cor(datum.iloc[:,7:], absolute = False)
            pcormats_split = [pcormat.clip(lower=0), pcormat.clip(upper=0)]
            pcormats.append(pcormats_split)
        graphs = []
        for i in range(len(pcormats)):
            split_graphs = [Graph(scipy.sparse.lil_matrix(pcormats[i][0]), directed = False), # positive split
                            Graph(scipy.sparse.lil_matrix(abs(pcormats[i][1])), directed = False)] # negative split - make positive for edge weights
            remove_self_loops(split_graphs[0])
            remove_self_loops(split_graphs[1])
            graphs.append(split_graphs)

    else: # if not splitting positive and negative
        for datum in data:
            pcormat = partial_cor(datum.iloc[:,7:], absolute = absolute)
            pcormats.append(pcormat)
        graphs = []
        for i in range(len(pcormats)):
            graphs.append(Graph(scipy.sparse.lil_matrix(pcormats[i]), directed = False))
            remove_self_loops(graphs[i])

    if len(graphs) == 1:
        graphs = graphs[0] # so that we don't get a list if only one graph

    return graphs

def plot_strengths(graphs):
    """
    Plot distributions of node strengths of a list of graphs
    :param graphs: a list of graph-tool graph objects
    :return: a matplotlib figure
    """
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

def measure_net(graph):
    """
    Calculate basic network metrics for a graph
    :param graph: a graph-tool graph object
    :return: a tuple of a DataFrame of network metrics and a histogram of shortest distances
    """

    net_metrics = pd.DataFrame() # initialise dataframe to store results

    clustering = global_clustering(graph, weight=graph.ep.weight)[0]  # global weighted clustering

    weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
    for e in graph.edges(): # for all the graph edges
        weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
    graph.ep['weight_inv'] = weight_inv # set inverted weights as graph edge property
    # hist = distance_histogram(graph, weight = graph.ep.weight_inv, bins=[0,0.5]) # get distance histogram with inverted weights
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

    net_metrics['clustering'] = [clustering]
    net_metrics['efficiency'] = [global_efficiency]
    net_metrics['L_obs'] = [L_obs]
    net_metrics['mean_eb'] = [np.mean(eb)]
    net_metrics['mean_vb'] = [np.mean(vb)]

    net_metrics = pd.DataFrame(net_metrics)

    return(net_metrics)

# Function to plot graphs in the following style:
# anatomical layout, given x and y coordinates
# nodes coloured by according to group membership, given group numbers (in x and y data)
# edges brighter, thicker and closer to top of layer with increasing weight
# requires graph-tool, matplotlib.cm
def draw_graph_anatomical(graph, positions, colours, absolute_edges, output_file):
    """
    Draw a graph with an anatomical layout
    Vertices coloured by group membership; stronger edges are brighter, thicker and closer to top of layer
    requires: graph-tool, matplotlib.cm
    :param graph: a graph-tool graph object
    :param positions: a pandas DataFrame with columns 'x' and 'y'
    :param colours: a dictionary with numerical keys corresponding to hex code colours
    :param absolute_edges: if True, scale edge widths but keeping absolute differences; otherwise use min-max scaling
    :param output_file: where to save the figure
    :return: a figure
    """
    def hex_to_rgb(hex_color): # helper function to convert hex colour to RGB tuple
        hex_color = hex_color.lstrip('#')
        return tuple(int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))

    # set vertex positions and colours:
    pos = graph.new_vertex_property("vector<double>") # create new vertx property for position
    color = graph.new_vertex_property("vector<double>") # create new
    for i, vertex in enumerate(graph.vertices()): # for each vertex
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y']) # get x and y positions
        color[vertex] = hex_to_rgb(colours[positions.loc[i, 'Region']]) # set colour by group

    # set edge colours:
    edge_color= graph.new_edge_property("vector<double>") # initialise new edge property to store edge colours
    for e in graph.edges(): # for each edge
        weight_val = graph.ep.weight[e] # acquire edge weight
        if 'pos' in output_file:
           colourmap = cm.Reds
        elif 'neg' in output_file:
            colourmap = cm.Blues
        else:
            colourmap = cm.viridis
        rgba = list(colourmap(weight_val)) # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2) # set alpha values, scaled
        edge_color[e] = rgba # set edge colour
    graph.ep["color"] = edge_color # set edge colours as edge property

    # set edge widths
    if absolute_edges:
        # scaled absolute values
        edge_width = graph.new_edge_property("double")
        edge_width.a = (graph.ep.weight.a**3)*20

    else:
        # min-max scaling
        edge_width = prop_to_size(graph, graph.ep.weight, mi = 0.1, ma = 4, power = 3, log = False)

    graph_draw( # draw the graph
        graph,
        pos=pos, # set vertex positions
        vertex_fill_color=color, # set vertex fill colours
        edge_color=graph.ep.color, # refer to edge colours
        edge_pen_width=edge_width,
        eorder=graph.ep.weight, # set edge order (largest on top)
        vertex_color = 'black', # set vertex stroke colour
        vertex_size = 10, # set vertex size
        output_size = (800,800), # set vertex colour
        output = output_file
    )

def euclidean_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def get_dists(positions, name):
    """
    Get 2D physical distances between all vertices
    :param positions: a pandas DataFrame with columns 'x' and 'y' and a column of vertex names
    :param name: a string giving the vertex type name
    :return: a pandas DataFrame with distances between all vertices
    """
    area_pairs = list(itertools.product(positions[name], repeat=2))  # get pairs of areas
    distances = []  # initialise distances
    for (area1, area2) in area_pairs:  # for each area in each pair
        x1, y1 = positions[positions[name] == area1][['x', 'y']].values[
            0]  # X and Y values for area 1
        x2, y2 = positions[positions[name] == area2][['x', 'y']].values[
            0]  # X and Y values for area 2
        distance = euclidean_distance(x1, y1, x2, y2)  # euclidean distance of pair
        distances.append({
            'pair': f"{area1}-{area2}",  # create data
            'area_1': area1,
            'area_2': area2,
            'distance': distance
        })
    distances = pd.DataFrame(distances)  # convert to dataframe
    return(distances)

def weight_by_dist(cormat, dists, plot = True):
    """
    Merge edge weights with their physical distances
    :param cormat: a pandas DataFrame of vertex correlations
    :param dists: a pandas DataFrame with distances between all vertices
    :return: a pandas DataFrame of weight-distance combinations, and a scatterplot
    """
    flatmat = cormat.reset_index().melt(id_vars='index', var_name='area_1',value_name='weight')  # convert correlation matrix to long form

    dists['area_1'] = dists['area_1'].astype(str)
    dists['area_2'] = dists['area_2'].astype(str)

    graph_dist = pd.merge(flatmat, dists, left_on=['index', 'area_1'], right_on=['area_2', 'area_1'], how='left')  # merge distance and weight data
    graph_dist['pair'] = graph_dist.apply(lambda row: tuple(sorted([row['index'], row['area_1']])), axis=1) # recreate 'pair' column by reordering areas
    graph_dist = graph_dist[graph_dist['index'] != graph_dist['area_1']]  # drop self-loops
    edge_dist = graph_dist.drop_duplicates(subset='pair').reset_index(drop=True)  # drop repeated pairs

    if plot:
        fig, axs = plt.subplots(1, 1, figsize=(8, 8))
        axs.scatter(edge_dist['distance'], edge_dist['weight'], alpha=0.8, s=1.5, c=edge_dist['weight'])  # plot of weight over distance
        axs.set_xlabel('Distance')
        axs.set_ylabel('pcor')

    return graph_dist

# Function to calculate SWP and delta for a graph, given the correlation matrix and distances, and observed clustering and path length of graph in question
# creates a distance-based lattice graph and a single random rewire for the null models
def SWP(cormat, dists, C_obs, L_obs,report = True):
    """
    Calculate small-world propensity and delta for a graph using a distance-based lattice and single random rewire as null models
    :param cormat: a pandas DataFrame of vertex correlations
    :param dists: a pandas DataFrame with distances between all vertices
    :param C_obs: the observed clustering coefficient for the graph corresponding to cormat
    :param L_obs: the observed path length for the graph corresponding to cormat
    :return: tuple: SWP, delta
    """

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
    lattice_weights = dists['weight'].sort_values(ascending = False).reset_index(drop=True) # get sorted weights
    distances_sort['weight'] = lattice_weights # add sorted weights to sorted distances

    # change distances object back to covariance matrix and fill with sorted weights
    areas = pd.unique(dists[['area_1', 'area_2']].values.ravel('K'))
    cormat_latt = pd.DataFrame(index=areas, columns=areas, dtype=float).fillna(1)
    for _, row in distances_sort.iterrows():
        cormat_latt.at[row['area_1'], row['area_2']] = row['weight']
        cormat_latt.at[row['area_2'], row['area_1']] = row['weight']

        cormat_latt.reindex(index=cormat.index, columns=cormat.columns) # reorder matrix to match original order


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

    if report:
        print(C_latt, C_rand, L_latt, L_rand) # sanity check

    # calculate divergences from lattice and random models
    delta_C = np.clip((C_latt - C_obs) / (C_latt - C_rand), 0, 1)
    delta_L = np.clip((L_obs - L_rand) / (L_latt - L_rand), 0, 1)

    # calculate SWP
    SWP = np.sqrt((delta_C**2 + delta_L**2) / 2) # calculate SWP
    alpha = np.arctan(delta_L / delta_C) # get angle of divergence vector
    delta = (((4 * alpha) / np.pi) - 1) # get delta value
    delta = np.clip(delta, -1, 1) # set range

    # put into DataFrame
    SWP_results = pd.DataFrame(np.array([[SWP, delta]]), columns=['SWP', 'delta'])

    return SWP_results

def hedges_g(data, population_mean):
    """
    Calculate Hedges' g for a one-sample test.
    :param: data (array-like): The sample data.
    :param: population_mean (float): The hypothesized population mean.
    :return: float: Hedges' g.
    """
    # sample statistics
    n = len(data)
    sample_mean = np.mean(data)
    sample_std = np.std(data, ddof=1)  # Use ddof=1 for sample standard deviation

    # get Cohen's d
    cohen_d = (sample_mean - population_mean) / sample_std

    # get the correction factor for Hedges' g
    correction_factor = 1 - (3 / (4 * n - 1))

    # get Hedges' g
    hedges_g = cohen_d * correction_factor
    return hedges_g

def draw_graph_betweenclust(graph, positions, absolute_edges, output_file):
    """
    Draw a graph where vertex size and colour scale with betweenness and local clustering
    :param graph: a graph-tool graph object
    :param positions: a pandas DataFrame with columns 'x' and 'y'
    :param absolute_edges: if True, scale edge widths but keeping absolute differences; otherwise use min-max scaling
    :param output_file: where to save the figure
    :return: a graph figure
    """
    # get inverted weights
    weight_inv = graph.new_edge_property("double") # create new edge property for inverted weights
    for e in graph.edges():
        weight_inv[e] = 1.0 / graph.ep.weight[e] # get inverted weights
    graph.ep['weight_inv'] = weight_inv

    # calculate local clustering coefficients
    clusts = local_clustering(graph, weight=graph.ep.weight).a
    clusts = clusts**7 * 200000 # scale for visualisation - try to get one white in syn, one black in control
    colors = cm.inferno(clusts)[:, :3]
    vertex_color = graph.new_vertex_property("vector<double>")
    vertex_color.set_2d_array(colors.T)

    # calculate betweennesss
    vb, _ = betweenness(graph, weight = graph.ep.weight_inv)
    vb_vals = vb.a
    vb_vals = vb_vals * 200 + 12 # scale for visualisation
    vertex_size = graph.new_vertex_property("double")
    vertex_size.a = vb_vals

    # set vertex positions
    pos = graph.new_vertex_property("vector<double>")
    for i, vertex in enumerate(graph.vertices()):
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y'])

    # set edge colours
    edge_color= graph.new_edge_property("vector<double>") # initialise new edge property to store edge colours
    for e in graph.edges():
        weight_val = graph.ep.weight[e] # acquire edge weight
        rgba = list(cm.viridis(weight_val)) # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2) # set alpha values, scaled
        edge_color[e] = rgba
    graph.ep["color"] = edge_color # set edge colours

    # set edge widths
    if absolute_edges:
        # scaled absolute values
        edge_width = graph.new_edge_property("double")
        edge_width.a = (graph.ep.weight.a**3)*20
    else:
        # min-max scaling
        edge_width = prop_to_size(graph, graph.ep.weight, mi = 0.1, ma = 4, power = 3, log = False)

    # draw graphs
    graph_draw(
        graph,
        pos=pos, # set vertex positions
        vertex_fill_color = vertex_color,
        edge_color = graph.ep.color, # refer to edge colours
        edge_pen_width = edge_width,
        eorder = graph.ep.weight, # set edge order (largest on top)
        vertex_color = 'black',
        vertex_size = vertex_size,
        output_size = (800,800),
        output = output_file
    )

def draw_graph_difference(graphs, positions, absolute_edges, output_file):
    """
    Draw a graph from two input graphs such that the edges represent the difference in correlation, vertex colour represents difference in local clustering and vertex size represents the difference in betweenness
    :param graphs: a list of graph-tool graph objects
    :param positions: a pandas DataFrame with columns 'x' and 'y'
    :param absolute_edges: if True, scale edge widths but keeping absolute differences; otherwise use min-max scaling
    :param output_file: where to save the figure
    :return: a graph figure
    """
    for graph in graphs:
        # get inverted weights
        weight_inv = graph.new_edge_property("double") # create new edge property for inverted weights
        for e in graph.edges():
            weight_inv[e] = 1.0 / graph.ep.weight[e] # get inverted weights
        graph.ep['weight_inv'] = weight_inv

        # calculate local clustering coefficients
        graph.vp['local_clustering'] = local_clustering(graph, weight=graph.ep.weight)
        #calculate betweenness
        vb, eb = betweenness(graph, weight = graph.ep.weight_inv)
        graph.vp['v_between'] = vb

    # get difference values: absolute difference of synesthete - control
    edge_diffs = abs(graphs[1].ep.weight.a - graphs[0].ep.weight.a)
    clust_diffs = abs(graphs[1].vp.local_clustering.a - graphs[0].vp.local_clustering.a)
    between_diffs = abs(graphs[1].vp.v_between.a - graphs[0].vp.v_between.a)

    # create new graph and set properties
    graph = graphs[0].copy()

    # set weights
    graph.ep.weight.a = edge_diffs

    # set vertex colours
    colors = cm.inferno(np.log1p(clust_diffs)*20)[:, :3] # log-transform differences
    vertex_color = graph.new_vertex_property("vector<double>")
    vertex_color.set_2d_array(colors.T)

    # set vertex sizes
    vertex_size = graph.new_vertex_property("double")
    vertex_size.a = np.log1p(between_diffs)*210 + 8

    # set vertex positions
    pos = graph.new_vertex_property("vector<double>")
    for i, vertex in enumerate(graph.vertices()):
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y'])

    # set edge colours
    edge_color= graph.new_edge_property("vector<double>") # initialise new edge property to store edge colours
    for e in graph.edges():
        weight_val = graph.ep.weight[e] # acquire edge weight
        rgba = list(cm.viridis(weight_val)) # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2) # set alpha values
        edge_color[e] = rgba
    graph.ep["color"] = edge_color # set edge colours

    # set edge widths
    if absolute_edges:
        # scaled absolute values
        edge_width = graph.new_edge_property("double")
        edge_width.a = (graph.ep.weight.a**3)*20
    else:
        # min-max scaling
        edge_width = prop_to_size(graph, graph.ep.weight, mi = 0.1, ma = 4, power = 3, log = False)

    # draw graph
    graph_draw(
        graph,
        pos=pos, # set vertex positions
        vertex_fill_color = vertex_color,
        edge_color = graph.ep.color, # refer to edge colours
        edge_pen_width = prop_to_size(graph.ep.weight, 0.01, 4, power = 4, log = False), # set edge widths based on weights
        eorder = graph.ep.weight, # set edge order (largest on top)
        vertex_color = 'black',
        vertex_size = vertex_size,
        output_size = (800,800),
        output=output_file
    )
    
def permutation_test_metrics(data1, data2, observed_metrics, sign=None, n_permutations=1000, test_directions=None):
    """
    Perform a permutation test by randomizing group labels and computing metric differences.

    Parameters:
    - data1, data2: Pandas DataFrames containing observations for each group. Must have a 'Group' column.
    - n_permutations: Number of permutations to perform.

    Returns:
    - A DataFrame with metric differences for each iteration.
    """

    # Compute observed metrics
    observed_diffs = observed_metrics.iloc[0, :] - observed_metrics.iloc[1, :]  # Difference: Group1 - Group2
    print('Observed differences:')
    print(observed_diffs.to_frame().T) # print the observed difference

    # Store permutation results
    null_diffs = pd.DataFrame(columns=observed_diffs.index, index=range(n_permutations))
    combined_data = pd.concat([data1, data2])
    group_labels = combined_data["Group"].unique()
    group1_label, group2_label = group_labels

    for i in range(n_permutations):
        # Shuffle group labels
        permuted_data = combined_data.copy()
        permuted_data["Group"] = np.random.permutation(permuted_data["Group"])

        # Generate graphs from permuted data
        # For positive or negative graphs, split by sign then take the first or second graph for each group
        if sign == 'Positive':
            graphs = data_to_graphs(permuted_data, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][0], graphs[1][0]]
        elif sign == 'Negative':
            graphs = data_to_graphs(permuted_data, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][1], graphs[1][1]]
        # For absolute weights, don't split
        elif sign is None:
            graphs_perm = data_to_graphs(permuted_data, groupname=group2_label)

        # Compute metrics for permuted groups
        permuted_metrics = []  # initialise metric list
        for graph in graphs_perm:  # for each graph and group
            metrics = measure_net(graph)  # get the network metrics
            permuted_metrics.append(metrics)  # add net metrics to list

        permuted_metrics = pd.concat(permuted_metrics,ignore_index=True)
        permuted_diffs = permuted_metrics.iloc[0, :] - permuted_metrics.iloc[1, :]

        # Store null differences
        null_diffs.iloc[i, :] = permuted_diffs.values

    print('Mean null differences:')
    print(null_diffs.mean())

    # Two-tailed P-values, quantile method
    #lower_percentile = diff_df.quantile(0.025)
    #upper_percentile = diff_df.quantile(0.975)
    #p_values = ((observed_diff < lower_percentile) | (observed_diff > upper_percentile)).astype(int)

    # Two-tailed P-values, proportional method
    # p_values = (np.abs(null_diffs) >= np.abs(observed_diffs)).mean()

    # One-tailed P-values
    p_values = pd.DataFrame(columns=observed_diffs.index, index=range(1)) # Store P-values
    effect_sizes = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # Store Cohen's d
    for metric in observed_diffs.index:
        direction = test_directions.get(metric, "greater")  # Default to greater if not specified
        null_distribution = null_diffs[metric].dropna()  # Drop NaNs if any

        # Compute effect sizes
        null_mean = null_distribution.mean()
        null_std = null_distribution.std()
        cohen_d = (observed_diffs[metric] - null_mean) / null_std if null_std > 0 else np.nan
        effect_sizes[metric] = cohen_d

        if direction == "greater":
            # quantile method:
            #threshold = np.percentile(null_diffs[metric], 95)
            #p_values[metric] = np.mean(observed_diff[metric] > threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] >= observed_diffs[metric]).mean()
        elif direction == "less":
            # quantile method:
            #threshold = np.percentile(null_diffs[metric], 5)
            #p_values[metric] = np.mean(observed_diff[metric] < threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] <= observed_diffs[metric]).mean()
        else:
            raise ValueError(f"Invalid test direction for metric '{metric}': {direction}")

    print('P-values:')
    print(p_values)
    print('Effect sizes:')
    print(effect_sizes)

    return p_values, effect_sizes

def permutation_test_SWP(data1, data2, distances, observed_metrics, sign=None, n_permutations=1000, test_directions=None):
    """
    Perform a permutation test by randomizing group labels and computing SWP differences.

    Parameters:
    - data1, data2: Pandas DataFrames containing observations for each group. Must have a 'Group' column.
    - n_permutations: Number of permutations to perform.

    Returns:
    - A DataFrame with metric differences for each iteration.
    """
    # Compute observed metrics
    observed_diffs = observed_metrics.iloc[0, :] - observed_metrics.iloc[1, :]  # Difference: Group1 - Group2
    print('Observed differences:')
    print(observed_diffs.to_frame().T) # print the observed difference

    # Store permutation results
    null_diffs = pd.DataFrame(columns=observed_diffs.index, index=range(n_permutations))
    combined_data = pd.concat([data1, data2])
    group_labels = combined_data["Group"].unique()
    group1_label, group2_label = group_labels

    for i in range(n_permutations):
        # Shuffle group labels
        permuted_data = combined_data.copy()
        permuted_data["Group"] = np.random.permutation(permuted_data["Group"])

        # Split data into groups
        permuted_data1 = permuted_data[permuted_data['Group'] == group1_label]
        permuted_data2 = permuted_data[permuted_data['Group'] == group2_label]

        # Generate graphs from permuted data
        # For positive or negative graphs, split by sign then take the first or second graph for each group
        if sign == 'Positive':
            # Generate positive pcormats
            pcormat1 = partial_cor(permuted_data1.iloc[:, 7:], absolute=False).clip(lower = 0)
            pcormat2 = partial_cor(permuted_data2.iloc[:, 7:], absolute=False).clip(lower = 0)
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs = data_to_graphs(permuted_data, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][0], graphs[1][0]] # select positive graphs
        elif sign == 'Negative':
            # Generate negative pcormats
            pcormat1 = abs(partial_cor(permuted_data1.iloc[:, 7:], absolute=False).clip(upper = 0))
            pcormat2 = abs(partial_cor(permuted_data2.iloc[:, 7:], absolute=False).clip(upper = 0))
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs = data_to_graphs(permuted_data, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][1], graphs[1][1]]
        # For absolute weights, don't split
        elif sign is None:
            # Generate absolute pcormats
            pcormat1 = partial_cor(permuted_data1.iloc[:, 7:])
            pcormat2 = partial_cor(permuted_data2.iloc[:, 7:])
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs_perm = data_to_graphs(permuted_data, groupname=group2_label)

        # Get C_obs and L_obs for each permuted graph
        permuted_metrics = []  # initialise metric list
        for graph in graphs_perm:  # for each graph and group
            metrics = measure_net(graph)  # get the network metrics
            permuted_metrics.append(metrics)  # add net metrics to list

        permuted_metrics = pd.concat(permuted_metrics,ignore_index=True)
        C_obs = permuted_metrics['clustering']
        L_obs = permuted_metrics['L_obs']

        # Calculate SWP and delta for each permuted dataset
        SWP1 = SWP(pcormat1, dists1, C_obs[0], L_obs[0], report=False)
        SWP2 = SWP(pcormat2, dists2, C_obs[1], L_obs[1], report=False)

        permuted_diffs = SWP1 - SWP2

        # Store null differences
        null_diffs.iloc[i, :] = permuted_diffs.values

    print('Mean null differences:')
    print(null_diffs.mean())

    # One-tailed P-values
    p_values = pd.DataFrame(columns=observed_diffs.index, index=range(1)) # Store P-values
    effect_sizes = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # Store Cohen's d
    for metric in observed_diffs.index:
        direction = test_directions.get(metric, "greater")  # Default to greater if not specified
        null_distribution = null_diffs[metric].dropna()  # Drop NaNs if any

        # Compute effect sizes
        null_mean = null_distribution.mean()
        null_std = null_distribution.std()

        cohen_d = (observed_diffs[metric] - null_mean) / null_std if null_std > 0 else np.nan
        effect_sizes[metric] = cohen_d

        if direction == "greater":
            # quantile method:
            # threshold = np.percentile(null_diffs[metric], 95)
            # p_values[metric] = np.mean(observed_diff[metric] > threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] >= observed_diffs[metric]).mean()
        elif direction == "less":
            # quantile method:
            # threshold = np.percentile(null_diffs[metric], 5)
            # p_values[metric] = np.mean(observed_diff[metric] < threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] <= observed_diffs[metric]).mean()
        else:
            raise ValueError(f"Invalid test direction for metric '{metric}': {direction}")

    print('P-values:')
    print(p_values)
    print('Effect sizes:')
    print(effect_sizes)

    return p_values, effect_sizes
