# %%
# load the packages
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import numpy.typing as npt
import pandas as pd
import circlify as circ
import igraph as ig
import matplotlib.pyplot as plt
import os


# %%
# Small functions that return some basic information about the original data
def View_cloneID(cloneID_cellID: dict) -> list:
    '''
    The function called to return cloneIDs of the data.
    '''
    cloneIDs = list(cloneID_cellID.keys())
    return cloneIDs


# %%
# Some small functions that return basic information about the clustering process 

def n_layers(subclones: dict, cloneID: int) -> int:
    '''
    The number of iteration for AGNES clustering. 
    '''
    storage = subclones[cloneID]  # define the storage data for cloneID
    x = len(storage.keys())
    return x


def n_clusters(subclones: dict, cloneID: int) -> list:
    '''
    The number of clsuters for each iteration. 
    '''
    storage = subclones[cloneID]  # define the storage data for cloneID
    a = []
    for i in storage.keys():
        a.append(len(storage[i][0]))
    return a


def n_clusters_name(subclones: dict, cloneID: int) -> list:
    '''
    The name of each cluster for each iteration. 
    '''
    storage = subclones[cloneID]  # define the storage data for cloneID
    b = []
    for i in storage.keys():
        b.append(list(storage[i][0].keys()))
    return b


# %%
# Transformation of clustering results for visualization 

# Choice 1: circles
def circle(storage: dict, times: int, pointer = 0) -> dict:
    '''
    Change the clustering results to nested dictionary data structure,
    so that the user can do the visualization as nested circles. 
    This function is mainly based on recursive algorithm. 
    '''
    if times == 1:  # simplest situation 
        cluster = storage[times][1][pointer]
        circles = []
        count = 0
        for key in cluster:  # only consider id and datum for each cell 
            circles.append({})
            # add cell index with cellID 
            circles[count]['id'] = [key, storage[0][0][key][0]]  
            circles[count]['datum'] = len(storage[times-1][0][key])
            count += 1
        return circles
    
    cluster = storage[times][1][pointer]  # common situations 
    circles = []
    count = 0  # used as index to locate key (clu_name) data to 
    for key in cluster:  # consider id, datum and children clusters 
        circles.append({})
        circles[count]['id'] = key
        circles[count]['datum'] = len(storage[times-1][0][key])  
        circles[count]['children'] = circle(storage, times=times-1, pointer=key)
        count += 1

    return circles


def clone_circle(subclones: dict, cloneID: int) -> dict:
    '''
    The function called to create nested dictionary data for circles visualization. 
    Mainly used for one clone data. 
    '''
    storage = subclones[cloneID]  # extract out the clone data 

    if len(storage) == 1:  # if only one cell contained within this clone 
        id = [0, storage[len(storage)-1][0][0][0]]
        return [{'id': id, 'datum': 1}]
    else:
        # create a wrapper for the whole clone, then add circles data 
        circles = [{'id': 0, 'datum': len(storage[len(storage)-1][0][0]), 'children': []}]
        sub = circle(storage=storage, times=len(storage)-1, pointer = 0)
        circles[0]['children'] = sub
        return circles


def prepare_circle(load: pd.DataFrame, subclones: dict) -> dict:
    '''
    The function called to prepare circle plot visualization data.
    '''
    # prepare to record data
    cellID = load.cellID.unique()
    circle = [{'id': 0, 'datum': len(cellID), 'children': []}]
    # transfer clustering results to the data that can be used for circle plot visualization
    mediate = []
    for cloneID in subclones.keys():
        single_clone_circles = clone_circle(subclones, cloneID)
        mediate.extend(single_clone_circles)

    circle[0]['children'] = mediate

    return circle


# %%
def circle_allclones(clone_circle: list):
    '''
    The circle packing plot shows the clones only
    '''
    # create the object
    circles = circ.circlify(
        clone_circle,
        show_enclosure=False,
        target_enclosure=circ.Circle(x=0, y=0, r=1)
    )

    # set the figure property
    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_title('Circle packing plot of the clones', fontsize=30)
    # Remove the axes
    ax.axis('off')
    # set the boundaries off axes
    lim = []
    for c in circles:
        lim.append(max(abs(c.x) + c.r, abs(c.y) + c.r))
    lim = max(lim)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)

    # plot
    for c in circles:
        x, y, r = c
        if c.level == 2:
            ax.add_patch(plt.Circle((x, y), r, alpha=0.5,
                         linewidth=2, color="mediumturquoise"))
    os.makedirs("./plot")
    plt.savefig("./plot/circle_all_clones.png")


# %%
def circle_allclones_subclones(clone_circle: list):
    '''
    The circle packing plot shows the full single-cell phylogeny of all clones
    '''
    # create the object
    circles = circ.circlify(
        clone_circle,
        show_enclosure=False,
        target_enclosure=circ.Circle(x=0, y=0, r=1)
    )

    # set the figure property
    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_title(
        'Circle packing plot of the full single-cell phylogeny', fontsize=30)
    # Remove the axes
    ax.axis('off')
    # set the boundaries off axes
    lim = []
    for c in circles:
        lim.append(max(abs(c.x) + c.r, abs(c.y) + c.r))
    lim = max(lim)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)

    # find the highest level of all circles
    levels = set()
    for c in circles:
        levels.add(c.level)
    maxlevel = max(levels)

    # set the color selections
    colors = ["lightblue", "lightcyan", "mediumturquoise", "#69b3a2", "turquoise", "darkcyan",
              "darkslategrey", "skyblue", "deepskyblue", "steelblue", "dodgerblue", "slategrey", "azure"]
    # plot
    for c in circles:
        x, y, r = c
        # plot without the largest circle that wrapped all clones
        if c.level > 1 and c.level < len(colors)-1:
            for i in range(2, maxlevel+1):
                if c.level == i:
                    ax.add_patch(plt.Circle((x, y), r, alpha=0.5,
                                 linewidth=2, color=colors[i-2]))
    plt.savefig("./plot/circle_allclones_subclones.png")

# %%
def circle_single_subclone(subclones: dict, cloneID: int):
    single_clone_circle = clone_circle(subclones, cloneID)
    # create the object
    singlecircle = circ.circlify(
        single_clone_circle,
        show_enclosure=True,
        target_enclosure=circ.Circle(x=0, y=0, r=1)
    )
    # set the figure property
    fig, ax = plt.subplots(figsize=(14, 14))
    ax.set_title('Full circle packing plot of clone ' +
                 str(cloneID), fontsize=30)
    # Remove the axes
    ax.axis('off')
    # set the boundaries off axes
    lim = []
    for c in singlecircle:
        lim.append(max(abs(c.x) + c.r, abs(c.y) + c.r))
    lim = max(lim)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)

    # find the highest level of all circles
    levels = set()
    for c in singlecircle:
        levels.add(c.level)
    maxlevel = max(levels)

    # set the color selections
    colors = ["lightblue", "lightcyan", "mediumturquoise", "#69b3a2", "turquoise", "darkcyan",
              "darkslategrey", "skyblue", "deepskyblue", "steelblue", "dodgerblue", "slategrey", "azure"]
    # plot
    for c in singlecircle:
        x, y, r = c
        # plot with the largest circle that wrapped the clone
        if c.level == 1:
            ax.add_patch(plt.Circle((x, y), r, alpha=0.5,
                         linewidth=2, color=colors[1]))
        elif c.level > 1 and c.level < len(colors)-1:
            for i in range(2, maxlevel+1):
                if c.level == i:
                    ax.add_patch(plt.Circle((x, y), r, alpha=0.5,
                                 linewidth=2, color=colors[i]))
    plt.savefig("./plot/circle_single_subclone.png")



#
# # %%
# # Transformation of clustering results for visualization
#
# # Choice 2: graphs, nodes or trees
# def nodes(subclones: dict, cloneID: int) -> nx.Graph():
#     '''
#     The function that called to create a nx.Graph based on clustering results of one clone.
#     Can be used for several ways of visualization based on networkx or igraph package.
#     '''
#     storage = subclones[cloneID]  # extract out the clone data
#
#     g = nx.Graph()
#
#     # extract out the clustering situations/recordings of this clone
#     clusters = [storage[i][1] for i in storage.keys()]
#
#     if len(clusters) == 1:  # one node and no edges
#         g.add_node("0_0")
#         g.nodes["0_0"]['subset'] = 0
#         return g
#     else:
#         # decide edges to link different clusters
#         for i in range(len(clusters)):  # the times of iteration
#             clustering = clusters[i]
#             for j in range(len(clustering)):  # the name of cluster in this layer
#                 sub_cluster = clustering[j]
#                 for k in range(len(sub_cluster)):  # the name of subcluster contained
#                     cluster = sub_cluster[k]
#                     if i == 0:  # initialization (no edges)
#                         node_name = str(i) + '_' + str(cluster)
#                         g.add_node(node_name)
#                         g.nodes[node_name]['subset'] = i
#                     elif i > 0:  # add edges between previous layer and this layer of iteration
#                         node_name = str(i) + '_' + str(j)
#                         g.add_node(node_name)
#                         g.nodes[node_name]['subset'] = i  # 'subset' used for layout (networkx) package
#                         prev_node_name = str(i-1) + '_' + str(cluster)
#                         g.add_edge(node_name, prev_node_name)
#         return g

#
# # %%
# def Edgelists_Clone(subclones: dict, align: str):
#     '''
#     The function called to visualize edge lists for each clone.
#     '''
#     # set the figure parameter
#     fig_length = len(list(subclones.keys()))
#     fig, ax = plt.subplots(1, fig_length)
#     fig.set_size_inches(20, 15)
#
#
#     for idx, cloneID in enumerate(list(subclones.keys())):
#         # get the nodes lists based on clustering results
#         g = nodes(subclones, cloneID)
#
#         # decide colors
#         d = nx.coloring.greedy_color(g, strategy="largest_first")
#         a = list(d.keys())
#         a.sort()
#         node_color = [d[key] for key in a]
#
#         # decide postions of nodes on the graph
#         # or planar_layout(g), kamada_kawai_layout(g)
#         pos = nx.multipartite_layout(g, align=align)
#
#         # plotting
#         fig = plt.figure()
#         fig.set_size_inches(20, 15)
#         nx.draw_networkx(g, pos=pos, with_labels=False,
#                         node_color=node_color, ax=ax[idx])
#         ax[idx].set_title("Clone" + str(cloneID))

#
# # %%
# def Edgelist_OneClone(cloneID: int, subclones: dict, align: str):
#     '''
#     The function called to visualize the edge list for only one clone.
#     '''
#     # get the nodes lists based on clustering results
#     g = nodes(subclones, cloneID)
#
#     # decide colors
#     d = nx.coloring.greedy_color(g, strategy="largest_first")
#     a = list(d.keys())
#     a.sort()
#     node_color = [d[key] for key in a]
#
#     # decide postions of nodes on the graph
#     # or planar_layout(g), kamada_kawai_layout(g)
#     pos = nx.multipartite_layout(g, align=align)
#
#     # plotting
#     fig = plt.figure()
#     fig.set_size_inches(20, 15)
#     nx.draw_networkx(g, pos=pos, with_labels=False, node_color=node_color)
#



#
# # %%
# # Choice 3 Circle Visualization of Edge Lists
# def circular_tree(subclones: dict, cloneID: int):
#     '''
#     The circular tree shows the lineage tracing in selected cell clone
#     '''
#     # create the networkx object
#     G = nodes(subclones, cloneID)
#     # requires the Graphiz and pygraphiz packages
#     # initiate the plot
#     pos = nx.nx_agraph.graphviz_layout(G, prog="twopi", args="")
#     plt.figure(figsize=(7, 7))
#     nx.draw(G, pos, node_size=20, alpha=0.5,
#             node_color="blue", with_labels=False)
#     plt.axis("equal")
#     # plot
#     plt.title('Circular tree in cell clone ' + str(cloneID), fontsize=30)
#     plt.savefig("./plot/circular_tree.png")
#     plt.show()
#
#
# # %%
# # Choice 4 Lineage Tree Visualization
# def node_tree(subclones: dict, cloneID: int):
#     '''
#     The node tree shows the lineage tracing in the selected cell clone
#     '''
#     # create the networkx object
#     g = nodes(subclones, cloneID)
#     # visualize using igraph
#     # convert the networkx object into igraph object
#     x = ig.Graph.from_networkx(g)
#     # initiate the plot
#     fig, ax = plt.subplots(figsize=(30, 30))
#     x.vs["color"] = "deepskyblue"
#     layout = x.layout("tree", root=[x.vcount()-1])
#     ig.plot(x, bbox=(800, 800), margin=4, layout=layout, target=ax)
#     # plot
#     ax.set_title('Lineage tree in cell clone ' + str(cloneID), fontsize=30)
#     plt.savefig("./plot/node_tree.png")


# %%
