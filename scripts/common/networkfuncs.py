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
    """
    Calculates a partial correlation matrix from data.

    Parameters:
    - data (Pandas DataFrame): Data where columns are variables and rows are subjects.
    - absolute (bool): If true then returns absolute partial correlations.

    Returns:
    - DataFrame: Partial correlation matrix.
    """

    data_cov = data.cov().values  # get covariance values
    precision_matrix = np.linalg.inv(data_cov)  # invert cov matrix to get precision matrix
    data_pcorr = np.zeros_like(precision_matrix)  # initialise zero matrix to fill
    for i in range(precision_matrix.shape[0]):  # for every row
        for j in range(precision_matrix.shape[1]):  # and every column
            if i == j:  # for the diagonal entries
                data_pcorr[i, j] = 1.0  # add 1s
            else:  # for the off-diagonal entries
                data_pcorr[i, j] = (-precision_matrix[i, j] /
                                    np.sqrt(precision_matrix[i, i] * precision_matrix[j, j]))  # add partial correlation

    data_pcorr = pd.DataFrame(data_pcorr, index=data.columns, columns=data.columns)  # convert to pandas DataFrame

    if absolute:  # if requiring absolute partial correlations
        data_pcorr = abs(data_pcorr)  # get absolute values

    return data_pcorr


def data_to_graphs(data, data_cols, groupname=False, log=False, split_sign=False, absolute=True):
    """
    Converts data or data file into a graph-tool graph.

    Parameters:
    - data (DataFrame, str): Data or a filename containing data.
    - data_cols (slice): Indexes of columns that contain brain data.
    - groupname (str): Name of the experimental group for splitting the data.
    - log (bool): If true then log-transforms the data before graph generation.
    - split_sign (bool): If true then creates two separate graphs from positive and negative correlations.

    Returns:
    - graph-tool graph: A single or list of graphs.
    """

    if isinstance(data, str):  # if data argument is a string
        data = pd.read_csv(data)   # load in from filename

    if log:  # if wanting to log-transform data
        data.iloc[:, data_cols] = np.log(data.iloc[:, data_cols])  # log-transform the data columns

    if groupname:  # if segregating by group
        data = [data[data['Group'] == 'Control'],  # split into list of control and experimental
                data[data['Group'] == groupname]]

    else:  # if not segregating by group
        data = [data]  # put into list anyway

    pcormats = []  # initialise partial correlation matrices

    if split_sign:  # if splitting positive and negative
        for datum in data:  # for each data in list
            pcormat = partial_cor(datum.iloc[:, data_cols], absolute=False)  # get signed partial correlations
            pcormats_split = [pcormat.clip(lower=0), pcormat.clip(upper=0)]  # get positive and negative in list
            pcormats.append(pcormats_split)  # put into list
        graphs = []  # initialise graphs
        for i in range(len(pcormats)):  # for each correlation matrix
            split_graphs = [Graph(scipy.sparse.lil_matrix(pcormats[i][0]),
                                  directed=False),  # get graph from positive correlations
                            Graph(scipy.sparse.lil_matrix(abs(pcormats[i][1])),
                                  directed=False)]  # get graph from negative correlations (converted to positive)
            remove_self_loops(split_graphs[0])  # remove self-connections (1s in matrix)
            remove_self_loops(split_graphs[1])  # remove self-connections (1s in matrix)
            graphs.append(split_graphs)  # add to list

    else:  # if not splitting positive and negative
        for datum in data:  # for each data in list
            pcormat = partial_cor(datum.iloc[:, data_cols], absolute=absolute)  # get absolute partial correlations
            pcormats.append(pcormat)  # add to list
        graphs = []  # initialise graphs
        for i in range(len(pcormats)):  # for each partial correlation matrix
            graphs.append(Graph(scipy.sparse.lil_matrix(pcormats[i]), directed=False))  # generate graph
            remove_self_loops(graphs[i])  # remove self-connections (1s in matrix)

    if len(graphs) == 1:  # if there is one graph in list
        graphs = graphs[0]  # remove from list to return single graph object

    return graphs

def plot_strengths(graphs):
    """
    Plots distributions of node strengths for a list of graphs.

    Parameters:
    - graphs (list): A list of graph-tool graphs.

    Returns:
    - plt: A matplotlib figure.
    """

    strengths = []  # initialise list of strength values
    for graph in graphs:  # for each graph
        weight = graph.edge_properties["weight"]  # get edge weights
        strengths.append(graph.degree_property_map("total", weight=weight).a)  # get total weight for each vertex

    fig = plt.figure(figsize=(10, 6))  # initialise figure
    ax = fig.add_subplot(1, 1, 1)  # add a single panel
    counts1, bins1 = np.histogram(strengths[0], bins=30, density=True)  # get the histogram of control strengths
    plt.plot(bins1[:-1], counts1, label='Control')  # plot histogram
    counts2, bins2 = np.histogram(strengths[1], bins=30, density=True)  # get the histogram of synesthete strengths
    plt.plot(bins2[:-1], counts2, label='Syn')  # plot histogram
    fig.set_facecolor('white')  # make background white
    ax.set_facecolor('white')  # make background white
    plt.xlabel("Node Strength")  # x-axis title
    plt.ylabel("Density")  # y-axis title
    plt.legend()  # show the legend
    plt.title("Node Strength Distributions")  # set the title
    plt.show()  # show plot

def measure_net(graph, normalise=False, randomise=False, seed=91939):
    """
    Calculates basic network metrics for a graph.

    Parameters:
    - graph (graph-tool graph): A graph.

    Returns:
    - DataFrame: Network metrics for a graph.
    """

    net_metrics = pd.DataFrame()  # initialise dataframe to store results

    if normalise:  # if normalising weights
        norm_indices = np.where((graph.ep.weight.a != 1) & (graph.ep.weight.a != 0))[0]  # get off-diagonal indexes
        non_diag = graph.ep.weight.a[norm_indices]  # get the off-diagonal values
        weight_norm = ((non_diag - min(non_diag)) / (max(non_diag) - min(non_diag))) + 0.001  # normalise
        graph.ep.weight.a[norm_indices] = weight_norm  # replace weights with normalised weights

    if randomise:  # if randomising weights
        for i in range(10):
            shuffle_indices = np.where(graph.ep.weight.a != 1)[0]  # get off-diagonal indexes
            non_diag = graph.ep.weight.a[shuffle_indices]  # get the off-diagonal values
            np.random.seed(seed)  # set seed for consistency
            np.random.shuffle(non_diag)  # randomly shuffle
            graph.ep.weight.a[shuffle_indices] = non_diag  # replace

            # initialise lists for randomised graph metrics
            clusterings = []
            L_obss = []
            efficiencies = []
            vbs = []
            ebs = []

            clusterings.append(global_clustering(graph, weight=graph.ep.weight)[0])  # calculate global weighted clustering
            weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
            for e in graph.edges():  # for all the graph edges
                weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
            graph.ep['weight_inv'] = weight_inv  # set inverted weights as graph edge property
            total_efficiency = 0.0  # initialise total efficiency for graph
            count = 0  # initialise a count
            distances = shortest_distance(graph, weights=graph.ep.weight_inv)  # shortest distances between vertices
            distances_sum = 0  # initialise the sum of distances
            for j in range(graph.num_vertices()):  # for all vertices
                for k in range(j + 1, graph.num_vertices()):  # to all other vertices
                    if distances[j][k] != 0:  # if shortest distance is not 0
                        total_efficiency += 1.0 / distances[j][
                            k]  # add inverse of shortest distance to total efficiency
                        count += 1  # and add one count
                        distances_sum += distances[j][k]  # add shortest distance to sum
            L_obss.append(distances_sum / (
                        (graph.num_vertices() - 1) * graph.num_vertices()))  # characteristic (observed) path length
            efficiencies.append(total_efficiency / count)  # global efficiency of graph is mean of efficiency
            vbs.append(betweenness(graph, weight=graph.ep.weight_inv)[0])
            ebs.append(betweenness(graph, weight=graph.ep.weight_inv)[1])

        clustering = np.mean(clusterings)
        L_obs = np.mean(L_obss)
        global_efficiency = np.mean(efficiencies)
        vb = np.mean(vbs)
        eb = np.mean(ebs)

    else:
        clustering = global_clustering(graph, weight=graph.ep.weight)[0]  # calculate global weighted clustering
        weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
        for e in graph.edges():  # for all the graph edges
            weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
        graph.ep['weight_inv'] = weight_inv  # set inverted weights as graph edge property
        total_efficiency = 0.0  # initialise total efficiency for graph
        count = 0  # initialise a count
        distances = shortest_distance(graph, weights=graph.ep.weight_inv)  # shortest distances between vertices
        distances_sum = 0  # initialise the sum of distances
        for j in range(graph.num_vertices()):  # for all vertices
            for k in range(j + 1, graph.num_vertices()):  # to all other vertices
                if distances[j][k] != 0:  # if shortest distance is not 0
                    total_efficiency += 1.0 / distances[j][k]  # add inverse of shortest distance to total efficiency
                    count += 1  # and add one count
                    distances_sum += distances[j][k]  # add shortest distance to sum
        L_obs = distances_sum / ((graph.num_vertices() - 1) * graph.num_vertices())  # characteristic (observed) path length
        global_efficiency = total_efficiency / count  # global efficiency of graph is mean of efficiency
        vb, eb = betweenness(graph, weight=graph.ep.weight_inv)  # calculate edge and vertex betweenness

    # store results in DataFrame
    net_metrics['clustering'] = [clustering]
    net_metrics['efficiency'] = [global_efficiency]
    net_metrics['L_obs'] = [L_obs]
    net_metrics['mean_eb'] = [np.mean(eb)]
    net_metrics['mean_vb'] = [np.mean(vb)]
    net_metrics = pd.DataFrame(net_metrics)

    return(net_metrics)

def plot_network_distances(graph, label, colour):
    """
    Plots the distribution of network distances for a graph.

    Parameters:
    - graph (graph-tool graph): A graph.

    Returns:
    - plt: A matplotlib figure.
    """

    weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
    for e in graph.edges(): # for all the graph edges
        weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
    graph.ep['weight_inv'] = weight_inv # set inverted weights as graph edge property
    hist = distance_histogram(graph, weight = graph.ep.weight_inv, bins=[0,0.5])   # get distance histogram with inverted weights

    plt.plot(hist[1][:-1], hist[0], label=label, color=colour)  # create plot from histogram
    plt.xlabel("Distance")  # set x-axis label
    plt.ylabel("Density")  # set y-axis label
    plt.title("Distribution of shortest network distances")  # set title
    plt.legend()  # show legend
    plt.grid(alpha=0.3)  # make the grid slightly transparent


def draw_graph_anatomical(graph, positions, colours, absolute_edges, output_file, scale=(2, 20, 3)):
    """
    Draws a graph with an anatomical layout, with scaled edge widths and group-coloured vertices.

    Parameters:
    - graph (graph-tool graph): a graph.
    - positions (DataFrame): two-dimensional coordinate data for vertices with columns 'x' and 'y'.
    - colours (dict): numerical keys corresponding to hex code colours.
    - absolute_edges (bool): if True, scale edge widths but keeping absolute differences; otherwise use min-max scaling.
    - output_file (str): file to save the figure.
    - scale (tuple): values for minimum, maximum and power for scaling edge widths.

    Returns:
    - (plt): A matplotlib figure.
    """

    # Helper function to convert hex colour to RGB tuple
    def hex_to_rgb(hex_color):
        hex_color = hex_color.lstrip('#') # get rid of hash symbol
        return tuple(int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4)) # calculate RGB values

    # Set vertex positions and colours
    pos = graph.new_vertex_property("vector<double>")  # create new vertex property for position
    color = graph.new_vertex_property("vector<double>")  # create new vertex property for colour
    for i, vertex in enumerate(graph.vertices()):  # for each vertex
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y'])  # get x and y positions
        color[vertex] = hex_to_rgb(colours[positions.loc[i, 'Region']])  # set colour by group

    # Set edge colours
    edge_color= graph.new_edge_property("vector<double>")  # initialise new edge property to store edge colours
    for e in graph.edges():  # for each edge
        weight_val = graph.ep.weight[e]  # acquire edge weight
        if 'pos' in output_file:  # if creating a positive graph
           colourmap = cm.Reds  # use red colour gradient
        elif 'neg' in output_file:  # if creating a negative graph
            colourmap = cm.Blues  # use blue colour gradient
        else:  # otherwise
            colourmap = cm.viridis  # use the viridis colour gradient for absolute edges
        rgba = list(colourmap(weight_val))  # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2)  # set alpha values, scaled to enhance visibility
        edge_color[e] = rgba  # set edge colour
    graph.ep["color"] = edge_color  # set edge colours as edge property

    # Set edge widths
    if absolute_edges:  # if using absolute partial correlations
        edge_width = graph.new_edge_property("double")  # create new edge property for width
        edge_width.a = (graph.ep.weight.a ** scale[0]) * scale[1] # desired non-linear scaling
        #edge_width = prop_to_size(prop=graph.ep.weight, mi=scale[0], ma=scale[1], power=scale[2], log=False)  # use desired min-max scaling

    else:  # if using signed partial correlations
        edge_width = prop_to_size(prop=graph.ep.weight, mi=scale[0], ma=scale[1], power=scale[2], log=False)  # use desired min-max scaling

    # Draw the graph
    graph_draw(
        graph,  # input graph
        pos=pos,  # set vertex positions
        vertex_fill_color=color,  # set vertex fill colours
        edge_color=graph.ep.color,  # refer to edge colours
        edge_pen_width=edge_width,  # set edge widths
        eorder=graph.ep.weight,  # set edge order (largest on top)
        vertex_color='black',  # set vertex stroke colour
        vertex_size=11,  # set vertex size
        output_size=(800, 800),  # set plot size
        output=output_file,  # set where to save the plot,
        bg_color='white'
    )

def euclidean_distance(x1, y1, x2, y2):
    """
    Calculates the two-dimensional Euclidean distance between two points.

    Parameters:
    - x1 (float): x coordinate of first point.
    - y1 (float): y coordinate of first point.
    - x2 (float): x coordinate of second point.
    - y2 (float): y coordinate of second point.

    Returns:
    - float: Euclidean distance between two points.
    """
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def get_dists(positions, name):
    """
    Calculates two-dimensional physical distances between all vertices of a network.

    Parameters:
    - positions (DataFrame): Two-dimensional coordinate data for vertices with columns 'x' and 'y'.
    - name (string): Name of the vertex type.

    Returns:
    - DataFrame: two-dimensional physical distances between all vertices.
    """

    area_pairs = list(itertools.product(positions[name], repeat=2))  # get pairs of areas
    distances = []  # initialise distances
    for (area1, area2) in area_pairs:  # for each area in each pair of areas
        x1, y1 = positions[positions[name] == area1][['x', 'y']].values[0]  # x and y values for area 1
        x2, y2 = positions[positions[name] == area2][['x', 'y']].values[0]  # x and y values for area 2
        distance = euclidean_distance(x1, y1, x2, y2)  # get Euclidean distance of pair
        # add distance to list
        distances.append({
            'pair': f"{area1}-{area2}",
            'area_1': area1,
            'area_2': area2,
            'distance': distance
        })
    distances = pd.DataFrame(distances)  # convert to dataframe
    return distances


def permutation_test_dist(dists1, dists2, groupname, test_direction, n_permutations):
    # Combine the data
    combined_dists = pd.concat([dists1, dists2])
    combined_dists['Group'] = ['Control'] * len(dists1) + [groupname] * len(dists2)

    # Observed difference
    observed_diffs = dists1['weight'] - dists2['weight']
    observed_diff = observed_diffs.mean()
    print('Mean observed difference:')
    print(observed_diff)

    # Run permutations
    null_diffs = np.zeros(n_permutations)
    for i in range(n_permutations):
        permuted = combined_dists.copy()
        permuted['Group'] = np.random.permutation(permuted['Group'])

        mean_control = permuted[permuted['Group'] == 'Control']['weight'].mean()
        mean_treatment = permuted[permuted['Group'] == groupname]['weight'].mean()
        null_diffs[i] = mean_control - mean_treatment

    print('Mean null difference:')
    print(np.mean(null_diffs))

    # Calculate effect size
    null_mean = np.mean(null_diffs)
    null_std = np.std(null_diffs, ddof=1)
    effect_size = (observed_diff - null_mean) / null_std if null_std > 0 else np.nan

    # P-value
    if test_direction == "greater":
        p_value = np.mean(null_diffs >= observed_diff)
    elif test_direction == "less":
        p_value = np.mean(null_diffs <= observed_diff)
    else:
        raise ValueError("Invalid test direction")

    # Display results
    print('P-values:')
    print(p_value)
    print('Effect sizes:')
    print(effect_size)

    return p_value, effect_size


def weight_by_dist(cormat, dists, plot=True):
    """
    Combines graph edge weights with their physical distances.

    Parameters:
    - cormat (DataFrame): A correlation matrix containing the edge weights used to construct a graph.
    - dists (DataFrame): Physical distances between all pairs of vertices.
    - plot (bool): If True, plot the edge weights as a function of distance.

    Returns:
    - DataFrame: Merged weight-distance combinations.
    - plt: Scatterplot of edge weights as a function of distance.
    """

    # Prepare correlation matrix for merging
    flatmat = cormat.reset_index().melt(id_vars='index',
                                                 var_name='area_2',
                                                 value_name='weight')
    flatmat = flatmat.rename(columns={'index': 'area_1'})  # rename index to area_1
    flatmat['area_1'] = flatmat['area_1'].astype(int)  # ensure area names are integers
    flatmat['area_2'] = flatmat['area_2'].astype(int)  # ensure area names are integers
    flatmat['area_pair'] = flatmat.apply(lambda row: tuple(sorted([row['area_1'], row['area_2']])), axis=1)  # create unique pair combinations by ordering
    flatmat = flatmat.drop_duplicates(subset='area_pair').reset_index(drop=True)  # drop duplicate pairs
    flatmat = flatmat[flatmat['area_1'] != flatmat['area_2']]  # drop self-loops
    flatmat.drop(columns=['area_1', 'area_2'], inplace=True)  # drop area columns

    # Prepare distances for merging
    dists['area_1'] = dists['area_1'].astype(int)  # ensure area names are integers
    dists['area_2'] = dists['area_2'].astype(int)  # ensure area names are integers
    dists['area_pair'] = dists.apply(lambda row: tuple(sorted([row['area_1'], row['area_2']])), axis=1)  # create unique pair combinations by ordering
    dists = dists.drop_duplicates(subset='area_pair').reset_index(drop=True)  # drop duplicate pairs
    dists = dists[dists['area_1'] != dists['area_2']]  # drop self-loops
    dists.drop(columns=['area_1', 'area_2'], inplace=True)  # drop area columns

    # Merge weights with distances
    graph_dist = pd.merge(flatmat, dists, on=['area_pair'], how='inner')  # merge by the unique area pairs
    graph_dist['area_1'] = graph_dist['area_pair'].str[0]  # re-acquire area_1
    graph_dist['area_2'] = graph_dist['area_pair'].str[1]  # re-acquire area_2
    graph_dist['area_pair'] = graph_dist['area_pair'].apply(lambda x: f"{int(x[0])}-{int(x[1])}")  # convert back to string with dash

    if plot:  # if a plot of weight-distance combinations is desired
        fig, axs = plt.subplots(1, 1, figsize=(8, 8))  # initialise plot
        axs.scatter(graph_dist['distance'], graph_dist['weight'],  # plot of weight over distance
                    alpha=0.8, s=1.5, c=graph_dist['weight'])  # set transparency and colour
        axs.set_xlabel('Distance')  # set x-axis label
        axs.set_ylabel('pcor')  # set y-axis label

    return graph_dist


def SWP(cormat, dists, C_obs, L_obs, report=True):
    """
    Calculates small-world propensity (SWP) and delta of a graph.
    Uses a distance-based lattice and a single random rewire as null models, preserving edge weight distribution.

    Parameters:
    - cormat (DataFrame): A correlation matrix containing the edge weights used to construct the graph.
    - dists (DataFrame): Physical distances between all pairs of vertices of the graph.
    - C_obs (float): The observed clustering coefficient of the graph.
    - L_obs (float): The observed path length of the graph.
    - report (bool): If True, prints C and L for null models and delta C and delta L.

    Returns:
    - DataFrame: SWP and delta values.
    """

    # Helper function to return characteristic path length of a graph
    def calculate_path_length(graph):
        weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
        for e in graph.edges():  # for each edge in the graph
            weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
        graph.ep['weight_inv'] = weight_inv  # set as edge property
        distances = shortest_distance(graph, weights=graph.ep.weight_inv)  # get shortest distances
        distances_sum = 0  # initialise a sum of shortest distances
        for j in range(graph.num_vertices()):  # for all vertices
            for k in range(j + 1, graph.num_vertices()):  # to all other vertices
                if distances[j][k] != 0:  # if shortest distance is not 0
                    distances_sum += distances[j][k]  # add distance to sum
        L = distances_sum / ((graph.num_vertices() - 1) * graph.num_vertices())  # characteristic path length

        return L

    # Generate lattice graph null model
    distances_sort = dists.sort_values(by=['distance'], ascending=[True])  # sort distances from shortest to largest
    distances_sort = distances_sort[['area_1', 'area_2', 'distance']].reset_index(drop=True)  # get relevant columns
    lattice_weights = dists['weight'].sort_values(ascending=False).reset_index(drop=True)  # get sorted weights
    distances_sort['weight'] = lattice_weights  # combine sorted distances and sorted weights

    areas = pd.unique(dists[['area_1', 'area_2']].values.ravel('K'))  # get unique area names
    cormat_latt = pd.DataFrame(index=areas, columns=areas, dtype=float).fillna(1)  # initialise a lattice matrix
    for _, row in distances_sort.iterrows():  # for each sorted distance
        cormat_latt.at[row['area_1'], row['area_2']] = row['weight']  # add partial correlations to matrix
        cormat_latt.at[row['area_2'], row['area_1']] = row['weight']  # add the symmetric entry
    cormat_latt.reindex(index=cormat.index, columns=cormat.columns)  # reorder matrix to match original order

    graph_latt = Graph(scipy.sparse.lil_matrix(cormat_latt), directed=False)  # create lattice graph
    remove_self_loops(graph_latt)  # remove self-loops

    # Generate random graph null model
    C_rands = []  # initialise random clustering coefficients
    L_rands = []  # initialise random path lengths

    non_diag = [cormat.iloc[j, k]  # take entries from the correlation matrix
                for j in range(cormat.shape[0])  # from rows
                for k in range(cormat.shape[1])  # and from columns
                if j != k]  # but not the diagonal

    for i in range(10):  # multiple randomisations
        np.random.shuffle(non_diag)  # randomly shuffle the weights in place
        cormat_rand = cormat.copy()  # copy the structure of the input matrix
        v = 0  # set a counter for vertices
        for j in range(cormat.shape[0]):  # for all rows
            for k in range(cormat.shape[1]):  # to all columns
                if j != k:  # if not diagonal entry
                    cormat_rand.iloc[j, k] = non_diag[v]  # then populate with shuffled weights
                    v += 1  # add one to counter
        graph_rand = Graph(scipy.sparse.lil_matrix(cormat_rand), directed=False)  # create random graph
        C_rands.append(global_clustering(graph_rand, weight=graph_rand.ep.weight)[0])  # add clustering to list
        L_rands.append(calculate_path_length(graph_rand))  # add path length to list

    # Get network metrics of null models
    C_latt = global_clustering(graph_latt, weight = graph_latt.ep.weight)[0]  # global clustering, lattice graph
    C_rand = np.mean(C_rands)  # global clustering, random graph
    L_latt = calculate_path_length(graph_latt)  # path length, lattice graph
    L_rand = np.mean(L_rands)  # path length, random graph

    if report: # if we want outputs for sanity checking; these should be similar to observed metrics
        print('C_latt, C_rand, L_latt, L_rand')  # print metric names
        print(C_latt, C_rand, L_latt, L_rand)  # print null metrics

    # Calculate divergences of observed graph from lattice and random models
    delta_C = np.clip(((C_latt - C_obs) / (C_latt - C_rand)), a_min=0, a_max=1)  # clustering divergence
    delta_L = np.clip(((L_obs - L_rand) / (L_latt - L_rand)), a_min=0, a_max=1)  # path length divergence

    if report: # if we want outputs for sanity checking; should be between 0 and 1
        print('delta_C, delta_L')  # print divergence names
        print(delta_C, delta_L)  # print divergences

    # Calculate SWP and delta
    SWP = np.sqrt((delta_C**2 + delta_L**2) / 2)  # calculate SWP
    alpha = np.arctan(delta_L / delta_C)  # get angle of divergence vector
    delta = (((4 * alpha) / np.pi) - 1)  # get delta value
    delta = np.clip(delta, -1, 1)  # set range for delta
    SWP_results = pd.DataFrame(np.array([[SWP, delta]]), columns=['SWP', 'delta'])  # place into DataFrame

    return SWP_results

def cohend(data1, data2):
    """
    Calculates Cohen's d effect size for two the difference between two means.

    Parameters:
    - data1 (array-like): Data for one group.
    - data2 (array-like): Data for another group.

    Returns:
    - float: Cohen's d effect size.
    """

    n1 = len(data1)  # sample number for first group
    n2 = len(data2)  # sample number for second group
    mean1 = np.mean(data1)  # mean of first group
    mean2 = np.mean(data2)  # mean of second group
    std1 = np.std(data1, ddof=1)  # standard deviation for first group
    std2 = np.std(data2, ddof=1)  # standard deviation for second group

    stdpool = np.sqrt(((n1 - 1)*std1 + (n2 - 1)*std2) / (n1+n2-2))  # get pooled standard deviation
    cohen_d = (mean1 - mean2) / stdpool  # calculate Cohen's d

    return cohen_d

def hedges_g(data, population_mean):
    """
    Calculates Hedges' g for a one-sample test.

    Parameters:
    - data (array-like): The sample data.
    - population_mean (float): The population mean for comparison.

    Returns:
    - float: Hedges' g for a one-sample test.
    """

    n = len(data)  # sample size for data
    sample_mean = np.mean(data)  # get sample mean
    sample_std = np.std(data, ddof=1)  # get sample standard deviation

    cohen_d = (sample_mean - population_mean) / sample_std  # calculate one-sample Cohen's d
    correction_factor = 1 - (3 / (4 * n - 1))  # calculate correction factor for Hedges' g
    hedges_g = cohen_d * correction_factor  # calculate Hedges' g

    return hedges_g

def draw_graph_betweenclust(graph, positions, absolute_edges, output_file):
    """
    Draws a graph where vertex size and colour scale with betweenness and local clustering.

    Parameters:
    - graph (graph-tool graph): A graph.
    - positions (DataFrame): Two dimensional coordinate data with columns 'x' and 'y'.
    - absolute_edges (bool): If True, scale edge widths but keeping absolute differences; otherwise use min-max scaling.
    - output_file (str): Where to save the graph figure.

    Returns:
    - plt: A graph figure.
    """

    weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
    for e in graph.edges():  # for each edge in graph
        weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
    graph.ep['weight_inv'] = weight_inv  # set as edge property

    # Vertex colour by local clustering
    clusts = local_clustering(graph, weight=graph.ep.weight).a  # calculate weighted local clustering
    clusts = clusts**7 * 200000  # scale for visualisation
    colors = cm.inferno(clusts)[:, :3]  # set colours using inferno gradient
    vertex_color = graph.new_vertex_property("vector<double>")  # create a new vertex property for colour
    vertex_color.set_2d_array(colors.T)  # set vertex colours

    # Vertex size by betweenness
    vb, _ = betweenness(graph, weight = graph.ep.weight_inv)  # get vertex betweenness
    vb_vals = vb.a  # get betweenness values
    vb_vals = vb_vals * 200 + 12  # scale for visualisation
    vertex_size = graph.new_vertex_property("double")  # create new vertex property for size
    vertex_size.a = vb_vals  # set vertex sizes

    # Set vertex positions
    pos = graph.new_vertex_property("vector<double>")  # create new vertex property for positions
    for i, vertex in enumerate(graph.vertices()):  # for each vertex in graph
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y'])  # set two-dimensional position

    # Set edge colours
    edge_color = graph.new_edge_property("vector<double>")  # initialise new edge property to store edge colours
    for e in graph.edges():  # for each edge in graph
        weight_val = graph.ep.weight[e]  # acquire edge weight
        rgba = list(cm.viridis(weight_val))  # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2)  # set alpha values, scaled
        edge_color[e] = rgba  # set edge colour
    graph.ep["color"] = edge_color  # set edge colours to property

    # Set edge widths
    if absolute_edges:  # if using absolute edge weights
        edge_width = graph.new_edge_property("double")  # create new vertex property for edge width
        edge_width.a = (graph.ep.weight.a**3)*20  # non-linear scaling
    else:  # if not using absolute edge weights
        edge_width = prop_to_size(graph, graph.ep.weight, mi = 0.1, ma = 4, power = 3, log = False)  # min-max scaling

    # Draw graph
    graph_draw(
        graph,  # input graph
        pos=pos,  # set vertex positions
        vertex_fill_color=vertex_color,  # set vertex fill colour
        edge_color=graph.ep.color,  # set edge colours
        edge_pen_width=edge_width,  # set edge widths
        eorder=graph.ep.weight,  # set edge order (largest on top)
        vertex_color='black',  # set vertex stroke colour
        vertex_size=vertex_size,  # set vertex size
        output_size=(800,800),  # set size of plot
        output=output_file  # where to save the plot
    )

def draw_graph_difference(graphs, positions, absolute_edges, output_file):
    """
    Draws a graph representing the differences between two input graphs.
    Edges represent the difference in correlation, vertex colour local clustering and vertex size betweenness.

    Parameters:
    - graphs (list): A list of two graph-tool graphs.
    - positions (DataFrame): Two dimensional coordinate data with columns 'x' and 'y'.
    - absolute_edges (bool): If True, scale edge widths according to absolute weights; otherwise use min-max scaling.
    - output_file (str): Where to save the graph figure.

    Returns:
    - plt: A graph figure.
    """

    for graph in graphs:  # for each of two graphs
        weight_inv = graph.new_edge_property("double")  # create new edge property for inverted weights
        for e in graph.edges():  # for each edge in graph
            weight_inv[e] = 1.0 / graph.ep.weight[e]  # get inverted weights
        graph.ep['weight_inv'] = weight_inv  # add as property map

        graph.vp['local_clustering'] = local_clustering(graph, weight=graph.ep.weight)  # local clustering
        vb, eb = betweenness(graph, weight=graph.ep.weight_inv)  # betweenness
        graph.vp['v_between'] = vb  # set betweenness as property map

    # Get value differences betweenness two graphs
    edge_diffs = abs(graphs[1].ep.weight.a - graphs[0].ep.weight.a)  # edge weights
    clust_diffs = abs(graphs[1].vp.local_clustering.a - graphs[0].vp.local_clustering.a)  # local clustering
    between_diffs = abs(graphs[1].vp.v_between.a - graphs[0].vp.v_between.a)  # betweenness

    # Create new graph and set properties
    graph = graphs[0].copy()

    # Set edge weights
    graph.ep.weight.a = edge_diffs  # set weights as differences

    # Set vertex colours
    colors = cm.inferno(np.log1p(clust_diffs)*20)[:, :3]  # log-transform to enhance differences
    vertex_color = graph.new_vertex_property("vector<double>")  # set as vertex property
    vertex_color.set_2d_array(colors.T)  # add colours

    # Set vertex sizes
    vertex_size = graph.new_vertex_property("double")  # create vertex property for size
    vertex_size.a = np.log1p(between_diffs)*210 + 8  # log1p transform and scaling to enhance

    # Set vertex positions
    pos = graph.new_vertex_property("vector<double>")  # create vertex property for position
    for i, vertex in enumerate(graph.vertices()):  # for each vertex
        pos[vertex] = (positions.loc[i, 'x'], positions.loc[i, 'y'])  # set two-dimensional position

    # Set edge colours
    edge_color = graph.new_edge_property("vector<double>")  # create edge property for edge colours
    for e in graph.edges():  # for each edge
        weight_val = graph.ep.weight[e]  # acquire edge weight
        rgba = list(cm.viridis(weight_val))  # get colour from chosen matplotlib colourmap
        rgba[3] = (weight_val * 2)  # set alpha values, scaled to enhance opacity
        edge_color[e] = rgba  # set the colour
    graph.ep["color"] = edge_color  # set edge colours to property map

    # Set edge widths
    if absolute_edges:  # if using absolute edges
        edge_width = graph.new_edge_property("double")  # create edge property for edge width
        edge_width.a = (graph.ep.weight.a**3)*20  # non-linear scaling for width
    else:  # if using signed edges
        edge_width = prop_to_size(graph, graph.ep.weight, mi=0.1, ma=4, power=3, log=False)  # use min-max scaling

    # Draw graph
    graph_draw(
        graph, # input graph
        pos=pos,  # use vertex positions
        vertex_fill_color=vertex_color,  # use vertex colours
        edge_color=graph.ep.color,  # use to edge colours
        edge_pen_width=edge_width,  # use edge widths
        eorder=graph.ep.weight,  # set edge order (largest on top)
        vertex_color='black',  # set vertex stroke colour
        vertex_size=vertex_size,  # set vertex size
        output_size=(800, 800),  # set figure size
        output=output_file  # where to save the figure
    )
    
def permutation_test_metrics(data1, data2, data_cols, observed_metrics, sign=None, n_permutations=1000, test_directions=None, normalise=False, randomise=False):
    """
    Performs a permutation test of network metrics by randomizing group labels and computing metric differences.

    Parameters:
    - data1 (DataFrame): Wide-form raw data for the control group. Must have a 'Group' column.
    - data2 (DataFrame): Wide-form raw data for the experimental group. Must have a 'Group' column.
    - data_cols (slice): Indexes of columns in data that contain brain data.
    - observed_metrics (DataFrame): Wide-form network metrics for both groups, where control group is row 0 and experimental group is row 1.
    - sign (str): Whether to form graphs with positive weights ('Positive'), negative weights ('Negative'), or both (None).
    - n_permutations (int): Number of permutations to perform to construct the null distribution.
    - test_directions (dict): Directions for one-sided tests. Direction refers to control group.
    - normalise (bool): If True, use normalised edge weights for metric calculation.
    - randomise (bool): If True, use randomised edge weights for metric calculation.

    Returns:
    - A DataFrame with metric differences for each iteration.
    """

    # Get the observed group differences
    observed_diffs = observed_metrics.iloc[0, :] - observed_metrics.iloc[1, :]  # group 1 - group 2
    print('Observed differences:')  # print title
    print(observed_diffs.to_frame().T)  # print the observed differences

    # Run permutations
    null_diffs = pd.DataFrame(columns=observed_diffs.index, index=range(n_permutations))  # to store the null diffs
    combined_data = pd.concat([data1, data2])  # combine data into one DataFrame
    group_labels = combined_data["Group"].unique()  # get the group names
    group1_label, group2_label = group_labels  # get the group names separately

    for i in range(n_permutations):  # for each permutation
        permuted_data = combined_data.copy()  # copy the combined data
        permuted_data["Group"] = np.random.permutation(permuted_data["Group"])  # shuffle group labels

        # Generate graphs from permuted data
        if sign == 'Positive':  # if positive-only
            graphs = data_to_graphs(permuted_data, data_cols=data_cols, groupname=group2_label, split_sign=True)  # convert to graphs
            graphs_perm = [graphs[0][0], graphs[1][0]]  # take positive-only graphs
        elif sign == 'Negative':  # if negative-only
            graphs = data_to_graphs(permuted_data, data_cols=data_cols, groupname=group2_label, split_sign=True)  # convert to graphs
            graphs_perm = [graphs[0][1], graphs[1][1]]  # take negative-only graphs
        elif sign is None:  # if using absolute weights
            graphs_perm = data_to_graphs(permuted_data, data_cols=data_cols, groupname=group2_label)  # convert to graph

        # Compute metrics for permuted groups
        permuted_metrics = []  # initialise metric list
        for graph in graphs_perm:  # for each graph and group
            metrics = measure_net(graph, normalise=normalise, randomise=randomise)  # get the network metrics
            permuted_metrics.append(metrics)  # add net metrics to list

        permuted_metrics = pd.concat(permuted_metrics, ignore_index=True)  # put null metrics into DataFrame
        permuted_diffs = permuted_metrics.iloc[0, :] - permuted_metrics.iloc[1, :]  # calculate the difference
        null_diffs.iloc[i, :] = permuted_diffs.values  # store the null difference

    print('Mean null differences:')  # print title
    print(null_diffs.mean().to_frame().T)  # print the average null differences for comparison

    # Two-tailed P-values, quantile method
    #lower_percentile = diff_df.quantile(0.025)
    #upper_percentile = diff_df.quantile(0.975)
    #p_values = ((observed_diff < lower_percentile) | (observed_diff > upper_percentile)).astype(int)

    # Two-tailed P-values, proportional method
    # p_values = (np.abs(null_diffs) >= np.abs(observed_diffs)).mean()

    # One-tailed tests
    p_values = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # to store P-values
    effect_sizes = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # to store Cohen's d
    for metric in observed_diffs.index:  # for each metric
        direction = test_directions.get(metric, "greater")  # default test direction to greater if not specified
        null_distribution = null_diffs[metric].dropna()  # drop NaNs if any

        # Calculate effect sizes
        null_mean = null_distribution.mean()  # get the null distribution mean
        null_std = null_distribution.std()  # get the null distribution standard deviation
        cohen_d = (observed_diffs[metric] - null_mean) / null_std if null_std > 0 else np.nan  # Cohen's d
        effect_sizes[metric] = cohen_d  # store Cohen's d

        # Calculate P-values
        if direction == "greater":  # if direction of control group is greater
            # quantile method:
            #threshold = np.percentile(null_diffs[metric], 95)
            #p_values[metric] = np.mean(observed_diff[metric] > threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] >= observed_diffs[metric]).mean()  # P = null diffs > observed diffs
        elif direction == "less":
            # quantile method:
            #threshold = np.percentile(null_diffs[metric], 5)
            #p_values[metric] = np.mean(observed_diff[metric] < threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] <= observed_diffs[metric]).mean()  # P = null diffs < observed diffs
        else:
            raise ValueError(f"Invalid test direction for metric '{metric}': {direction}")  # if not greater or less

    # Print results
    print('P-values:')
    print(p_values)
    print('Effect sizes:')
    print(effect_sizes)

    return p_values, effect_sizes

def permutation_test_SWP(data1, data2, data_cols, distances, observed_metrics, sign=None, n_permutations=1000, test_directions=None):
    """
    Performs a permutation test of SWP by randomizing group labels and computing SWP differences.

    Parameters:
    - data1 (DataFrame): Wide-form raw data for the control group. Must have a 'Group' column.
    - data2 (DataFrame): Wide-form raw data for the experimental group. Must have a 'Group' column.
    - data_cols (slice): Indexes of columns in data that contain brain data.
    - observed_metrics (DataFrame): Wide-form network metrics for both groups, where control group is row 0 and experimental group is row 1.
    - sign (str): Whether to form graphs with positive weights ('Positive'), negative weights ('Negative'), or both (None).
    - n_permutations (int): Number of permutations to perform to construct the null distribution.
    - test_directions (dict): Directions for one-sided tests. Direction refers to control group.

    Returns:
    - A DataFrame with metric differences for each iteration.
    """

    # Get the observed group differences
    observed_diffs = observed_metrics.iloc[0, :] - observed_metrics.iloc[1, :]  # group 1 - group 2
    print('Observed differences:')  # print title
    print(observed_diffs.to_frame().T)  # print the observed difference

    # Run permutations
    null_diffs = pd.DataFrame(columns=observed_diffs.index, index=range(n_permutations))  # to store null differences
    combined_data = pd.concat([data1, data2])  # combine data into one DataFrame
    group_labels = combined_data["Group"].unique()  # get group names
    group1_label, group2_label = group_labels  # get individual group names

    for i in range(n_permutations):  # for each permutation
        permuted_data = combined_data.copy()  # copy the data
        permuted_data["Group"] = np.random.permutation(permuted_data["Group"])  # shuffle the group labels

        # Split data into groups
        permuted_data1 = permuted_data[permuted_data['Group'] == group1_label]  # group 1
        permuted_data2 = permuted_data[permuted_data['Group'] == group2_label]  # group 2

        # Generate correlation matrices, distances and graphs from permuted data
        if sign == 'Positive':  # if positive-only graphs
            # Generate positive pcormats
            pcormat1 = partial_cor(permuted_data1.iloc[:, data_cols], absolute=False).clip(lower = 0)
            pcormat2 = partial_cor(permuted_data2.iloc[:, data_cols], absolute=False).clip(lower = 0)
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs = data_to_graphs(permuted_data, data_cols=data_cols, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][0], graphs[1][0]] # select positive graphs
        elif sign == 'Negative':  # if negative-only graphs
            # Generate negative pcormats
            pcormat1 = abs(partial_cor(permuted_data1.iloc[:, data_cols], absolute=False).clip(upper = 0))
            pcormat2 = abs(partial_cor(permuted_data2.iloc[:, data_cols], absolute=False).clip(upper = 0))
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs = data_to_graphs(permuted_data, data_cols = data_cols, groupname=group2_label, split_sign=True)
            graphs_perm = [graphs[0][1], graphs[1][1]]
        elif sign is None:  # if using absolute weights
            # Generate absolute pcormats
            pcormat1 = partial_cor(permuted_data1.iloc[:, data_cols])
            pcormat2 = partial_cor(permuted_data2.iloc[:, data_cols])
            # Get weight-distance combinations
            dists1 = weight_by_dist(pcormat1, distances, plot=False)
            dists2 = weight_by_dist(pcormat2, distances, plot=False)
            # Generate graphs
            graphs_perm = data_to_graphs(permuted_data, data_cols=data_cols, groupname=group2_label)

        # Get clustering and path length for each permuted graph
        permuted_metrics = []  # initialise metric list
        for graph in graphs_perm:  # for each graph and group
            metrics = measure_net(graph)  # get the network metrics
            permuted_metrics.append(metrics)  # add net metrics to list
        permuted_metrics = pd.concat(permuted_metrics, ignore_index=True)  # put metrics into DataFrame
        C_obs = permuted_metrics['clustering']  # get clustering
        L_obs = permuted_metrics['L_obs']  # get path length

        # Calculate SWP and delta for each permuted dataset
        SWP1 = SWP(pcormat1, dists1, C_obs[0], L_obs[0], report=False)
        SWP2 = SWP(pcormat2, dists2, C_obs[1], L_obs[1], report=False)

        # Get the differences and store
        permuted_diffs = SWP1 - SWP2
        null_diffs.iloc[i, :] = permuted_diffs.values

    # Print average null differences
    print('Mean null differences:')
    print(null_diffs.mean().to_frame().T)

    # One-tailed P-tests
    p_values = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # to store P-values
    effect_sizes = pd.DataFrame(columns=observed_diffs.index, index=range(1))  # to store Cohen's d
    for metric in observed_diffs.index:  # for each metric
        direction = test_directions.get(metric, "greater")  # set test direction to greater if not specified
        null_distribution = null_diffs[metric].dropna()  # drop NaNs if any

        # Calculate effect sizes
        null_mean = null_distribution.mean()  # get means of null differences
        null_std = null_distribution.std()  # get standard deviations of null differences
        cohen_d = (observed_diffs[metric] - null_mean) / null_std if null_std > 0 else np.nan  # Cohen's d
        effect_sizes[metric] = cohen_d  # add to effect sizes

        # Calculate P-values
        if direction == "greater":  # if direction of control group is greater
            # quantile method:
            # threshold = np.percentile(null_diffs[metric], 95)
            # p_values[metric] = np.mean(observed_diff[metric] > threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] >= observed_diffs[metric]).mean()  # P = null diffs > observed diffs
        elif direction == "less":  # if direction of control group is lesser
            # quantile method:
            # threshold = np.percentile(null_diffs[metric], 5)
            # p_values[metric] = np.mean(observed_diff[metric] < threshold)
            # proportional method:
            p_values[metric] = (null_diffs[metric] <= observed_diffs[metric]).mean()  # P = null diffs < observed diffs
        else:
            raise ValueError(f"Invalid test direction for metric '{metric}': {direction}")

    # Print results
    print('P-values:')
    print(p_values)
    print('Effect sizes:')
    print(effect_sizes)

    return p_values, effect_sizes
