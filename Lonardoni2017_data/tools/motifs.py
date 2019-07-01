import copy
import networkx as nx

import numpy as np
import pylab as plt

import mp2


def partition(G):
    for u, v, d in G.edges(data=True):
        # print d
        d['weight'] = 1  # mapp[u,v]#mapt[u,v]*1./mapp[u,v]
    poss = nx.get_node_attributes(G, 'pos')

    parts = np.load('graph_partition.npy').astype(dict)
    count = 0.
    plt.figure()
    allnodes = []
    for com in set(parts.values()):
        count += 1.
        list_nodes = [nodes for nodes in list(parts.keys())
                      if parts[nodes] == com]
        allnodes.append(list_nodes)
        plt.plot(G.node[list_nodes[0]]['pos'][0], G.node[list_nodes[0]]['pos'][1], 'o',
                 color=plt.get_cmap('jet')(count * .95 / len(set(parts.values()))), label='%s-th community' % count)

        nx.draw_networkx_nodes(G, poss, list_nodes, node_size=20,
                               node_color=plt.get_cmap('jet')(count * .99 / len(set(parts.values()))))

        nx.draw_networkx_edges(G.subgraph(list_nodes), poss, alpha=0.5)

    # plt.legend()
    print("# communities", count)
    return allnodes


def reorder_motifs():
    import pickle
    motifs = pickle.load(open('/home/alonardoni/8192correlation/motifs_DB/motifs_DB_new_new.npy'))
    index = np.array(sorted(list(range(len(motifs))), key=lambda k: len(motifs[k][1].nodes())))
    motifs = np.array(motifs)[index]
    i = 0
    for ele in motifs:
        ele[0] = i
        i += 1
    pickle.dump(motifs, open('/home/alonardoni/8192correlation/motifs_DB/motifs_DB.npy', 'w'))


#    for ele in motifs:
#        print ele[0],len(ele[1].nodes())


def save_motifs(tau):
    import pickle
    pickle.dump(get_mymotifs(tau), open('/home/alonardoni/8192correlation/motifs_DB/motifs_DB_new_new_8.npy', 'w'))


def open_motifs(tau):
    import pickle
    motifs = pickle.load(open('/home/alonardoni/8192correlation/motifs_DB/motifs_DB.npy'))
    i = 0
    ele = motifs[0]
    while len(ele[1].nodes()) <= tau:
        i += 1
        if i == len(motifs):
            break
        else:
            ele = motifs[i]
    return motifs[:i]


def motif_hist(test_CC, ax1=0, ax2=0, show_plot=True):
    graph_comm = 0

    kol = {}
    if show_plot:

        if ax1 == 0:
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_axes([0.1, 0.1, 0.85, 0.15])
            ax2 = fig.add_axes([0.1, 0.3, 0.85, 0.6])

    n_groups = len(test_CC[list(test_CC.keys())[0]])
    bar_width = 0.8 / len(test_CC)
    index = []
    normalizer = np.ones(20)
    for i in sorted(test_CC.keys()):
        normalizer = np.ones(20)
        for ele in test_CC[i]:
            normalizer[len(ele[1].nodes())] += ele[2]
            # print "--------------------->",i,len(test_CC[i]),n_groups
        index = np.array([test_CC[i][k][0] for k in range(n_groups)])
        # test_CC[i][2][2]=test_CC[i][2][2]-test_CC[i][1][2]

        kol[i] = [count[2] * 1. / normalizer[len(count[1].nodes())] for count in test_CC[i]]

        if show_plot:

            if isinstance(i, int):
                ax2.bar(index + (i + 1) * bar_width,
                        [count[2] * 1. / normalizer[len(count[1].nodes())] for count in test_CC[i]], bar_width,
                        color=plt.get_cmap('autumn')(i * .9), label='%s-th component' % i)
            else:
                ax2.bar(index + 0.8 - graph_comm * bar_width,
                        [count[2] * 1. / normalizer[len(count[1].nodes())] for count in test_CC[i]], bar_width,
                        color=plt.get_cmap('autumn')(.9), label='all_graph%s' % graph_comm)
                graph_comm += 1

    if show_plot:
        for i in range(n_groups):
            tmppos = {}
            N = test_CC[0][i][1].nodes()
            k = 0
            NN = len(N)
            for j in N:
                tmppos[j] = [np.cos(2 * np.pi * k * 1. / NN) * .35 + (index[i] + 0.4),
                             np.sin(2 * np.pi * k * 1. / NN) * .35 - 1]
                k += 1
            nx.draw_networkx_edges(test_CC[0][i][1], tmppos, color='k', lw=2, ax=ax1)
        ax1.set_xlim([-1, n_groups + 1])
        ax2.set_xlim([-1, n_groups + 1])
        ax1.set_xlabel('Motif ID')
        ax2.set_ylabel('Occurrences')
        ax2.set_xticks(index + 0.4)
        ax2.set_xticklabels(index)
        ax2.legend(loc='best')
        ax1.set_xticks([])
        ax1.set_yticks([])
    from scipy.stats import ks_2samp
    plt.figure()
    if .9 <= sum(kol[0]) <= 1.1:
        comp = 1
    else:
        comp = len(normalizer) - sum(np.array(normalizer) == 1)
    # print comp
    pvm = np.zeros((len(list(test_CC.keys())) - 1, len(list(test_CC.keys())) - 1))

    if 0:
        control = np.array(kol['graph']) * 1. / comp
        for i in list(test_CC.keys()):
            for j in list(test_CC.keys()):
                if i != 'graph' and j != 'graph' and i != j:
                    pvm[i, j] = \
                        ks_2samp(np.cumsum(np.array(kol[i])) * 1. / comp, np.cumsum(np.array(kol[j])) * 1. / comp)[1]

                    #                if i>j:
                    #                    kol[i]=np.array(kol[i])
                    #                    kol[i][kol[i]==0]=1
                    #                    kol[j][kol[j]==0]=1
                    #                    print i,j,":",pvm[i,j],chisquare(kol[i],kol[j])
        pdiff_index = np.zeros(len(test_CC[0]))
        ndiff_index = np.zeros(len(test_CC[0]))
        for i in list(test_CC.keys()):
            if i != 'graph':
                print("i: ", i, "vs control\t p-value:", \
                    ks_2samp(np.cumsum(np.array(kol[i])) * 1. / comp, np.cumsum(control))[1])
            if i == 'graph':
                plt.plot(np.cumsum(control) + 2, 'r', lw=2, label='all_graph')
            else:
                plt.plot(np.cumsum(kol[i]) * 1. / comp + 2, 'b', lw=2, label='%s-th component' % i)
                diff = np.array(kol[i]) * 1. / comp - control
                for k in range(len(kol[i])):
                    if abs(diff[k]) > 0.5 * control[k]:
                        if diff[k] > 0:
                            pdiff_index[k] += 1
                        else:
                            ndiff_index[k] += -1
                            # print ndiff_index[k]
        # print sum(ndiff_index)
        plt.bar(np.arange(len(pdiff_index)), pdiff_index * 1. / (len(list(test_CC.keys())) - 1), edgecolor='b')
        plt.bar(np.arange(len(ndiff_index)), ndiff_index * 1. / (len(list(test_CC.keys())) - 1), edgecolor='r', color='r')

        plt.ylabel('# of nodes')
        plt.xlabel('motif ID')
        # plt.legend(loc='best')
        plt.suptitle('cumulative distribution')
        if show_plot:
            plt.figure()
            plt.matshow(pvm)
            plt.colorbar()

        for k in range(len(pdiff_index)):
            tmppos = {}
            N = test_CC[0][i][1].nodes()
            k = 0
            NN = len(N)
            for j in N:
                tmppos[j] = [np.cos(2 * np.pi * k * 1. / NN) * .35 + (index[i] + 0.4),
                             np.sin(2 * np.pi * k * 1. / NN) * .35 - 1]
                k += 1
            nx.draw_networkx_edges(test_CC[0][k][1], tmppos, color='k', lw=2, ax=ax1)
            pdiff_index[k] = 0


def motif_counts(tau, base, base_direction=False, mymotifs=True):
    """
    Returns dictionary keyed by graph motifs and values as the number of subgraph isomorphisms
    for the given motif counted in the base structure of the given GMM object.

    Parameters
    ----------
    base : Entire Graph
    base_direction : Consider directed motifs
    mymotifs : Used already saved motifs
    tau : An integer greater than or equal to 2, which designates the number of nodes in the
        largest graph in set of graph motifs used in the given model.

    Returns
    ----------
    subgraph_counts : A list of tuples with the following construction (index,motif,count), where
        "count" is the number of subgraph isomorphisms for the given motif counted in the base
        structure of the given GMM object and index is the motif index.
    """
    #    base=gmm.get_base()
    #    base_direction=base.is_directed()   # Check if GMM base is directed, motifs must match
    ubase = base
    base = nx.line_graph(base)
    if mymotifs:
        # motif_counts=get_mymotifs(tau,base_direction)
        motif_counts = open_motifs(tau)
    else:
        motif_counts = get_motifs(tau, base_direction)
    # Performing the counting of subgraph isomorphism for every motif given the base structure
    for motif in motif_counts:

        index = motif[0]
        if index == 1:
            # print index,"of",len(motif_counts)
            A = nx.adj_matrix(ubase)
            A3 = A.dot(A).dot(A)
            triangles = 0
            for i in range(A3.shape[0]):
                triangles += A3[i, i]

            motif_counts[index] = (index, motif_counts[index][1], triangles / 6)
        else:
            # print index,"of",len(motif_counts),
            # Calculate subgraph ismorphisms
            if base_direction:
                GM = nx.DiGraphMatcher(base, motif[1])
            else:
                GM = nx.GraphMatcher(base, nx.line_graph(motif[1]))
            # Count subgraph isomorphisms and update
            count = 0
            for _ in GM.subgraph_isomorphisms_iter():
                count += 1
            GM = nx.GraphMatcher(nx.line_graph(motif[1]), nx.line_graph(motif[1]))
            acount = 0
            for _ in GM.isomorphisms_iter():
                acount += 1

            if index == 2:
                count = (count - motif_counts[1]) / acount
            else:
                count = count * 1. / acount
            motif_counts[index] = (index, motif_counts[index][1], count)
            # print "elapsed:",time.time()-start
    return motif_counts

def get_mymotifs(tau):
    """
    Returns a list of tuples of all possible single component non-singleton subgraphs from a dyad to
    $G=\{V=tau,E=\frac{tau(tau-1)}{2}\}$, i.e. the K_{tau} complete graph.

    The list represents the set of graph moitifs used in a given model.  This function is
    primarily meant as a helper to motif_counts, but can be used for other user purposes.

    Parameters
    ----------
    tau : An integer greater than or equal to 2, which designates the number of nodes in the
        largest graph in set of graph motifs used in the given model.


    Returns
    ----------
    motifs : A list of tuples of the following construction (index,motif,0).  The final zero is used as
        is a placeholder and will set the counts of subgraph isomorphism for each motif in the GMM's
        base structure.
    """
    motifs = list()
    motif_index = 0
    for g in smart_graphs(tau):
        if len(g.nodes()) >= 3:
            motifs.append((motif_index, g, 0))
            motif_index += 1
    return motifs


def get_motifs(tau, directed_motifs):
    """
    Returns a list of tuples of all possible single component non-singleton subgraphs from a dyad to
    $G=\{V=tau,E=\frac{tau(tau-1)}{2}\}$, i.e. the K_{tau} complete graph.

    The list represents the set of graph moitifs used in a given model.  This function is
    primarily meant as a helper to motif_counts, but can be used for other user purposes.

    Parameters
    ----------
    tau : An integer greater than or equal to 2, which designates the number of nodes in the
        largest graph in set of graph motifs used in the given model.

    directed_motifs : A boolean designating whether the motifs should be directed.

    Returns
    ----------
    motifs : A list of tuples of the following construction (index,motif,0).  The final zero is used as
        is a placeholder and will set the counts of subgraph isomorphism for each motif in the GMM's
        base structure.
    """
    motifs = list()
    motif_index = 0
    for v in range(2, tau + 1):
        graphs = all_graphs(v, directed_motifs)
        for g in graphs:
            motifs.append((motif_index, g, 0))
            motif_index += 1
    return motifs


def smart_graphs(num_nodes):
    graphs = list()
    g = nx.Graph()
    g.add_edge(0, 1)
    graphs.append(g)
    for i in range(3, num_nodes + 1):
        graphs = recursive_worker(graphs, i)

    return graphs


def recursive_worker(graphs, num_nodes):
    def add_correct_edge(g, count):
        count += 1

        tmp = nx.complete_graph(num_nodes)

        if tmp.edges() == g.edges():
            count += -1
            if np.sum([nx.GraphMatcher(g, ele).is_isomorphic() for ele in graphs]) == 0:
                graphs.append(g)
            return
        else:
            for node in g.nodes():
                for pivot in tmp.nodes():
                    if pivot not in g.neighbors(node):
                        if not pivot == node:
                            # print node,pivot
                            complete_copy = copy.deepcopy(g)
                            complete_copy.add_edge(node, pivot)
                            if np.sum([nx.GraphMatcher(complete_copy, ele).is_isomorphic() for ele in graphs]) == 0:
                                graphs.append(complete_copy)
                                add_correct_edge(complete_copy, count)

    count = 0
    for g in graphs:
        if len(g) == num_nodes - 1:
            # print "new graph"
            add_correct_edge(g, 1)
            count += -1
    # print "num_nodes",num_nodes,"----------------->",abs(count),len(graphs)
    return graphs


def very_all_graphs(num_nodes):
    def add_correct_edge(g, count):
        count += 1
        N = len(g.nodes())
        tmp = nx.complete_graph(N)
        # print len(graphs),N,num_nodes,count,len(g.edges())*1./len(g)/(len(g)-1)*2
        if tmp.edges() == g.edges():
            if N + 1 > num_nodes:
                count += -1
                return
            else:
                for h in graphs:
                    for node in h.nodes():
                        if node != N:
                            complete_copy = copy.deepcopy(h)
                            complete_copy.add_edge(node, N)
                            if np.sum([nx.GraphMatcher(complete_copy, ele).is_isomorphic() for ele in graphs]) == 0:
                                graphs.append(complete_copy)
                                add_correct_edge(complete_copy, count)
        else:
            for ele in tmp.edges():
                if ele not in g.edges():
                    if (ele[1], ele[0]) not in g.edges():
                        complete_copy = copy.deepcopy(g)
                        complete_copy.add_edge(ele[0], ele[1])
                        if np.sum([nx.GraphMatcher(complete_copy, ele).is_isomorphic() for ele in graphs]) == 0:
                            graphs.append(complete_copy)
                            add_correct_edge(complete_copy, count)
                        else:
                            add_correct_edge(complete_copy, count)

    graphs = list()
    g = nx.Graph()
    g.add_edge(0, 1)
    graphs.append(g)
    add_correct_edge(g, 1)

    return graphs


def all_graphs(num_nodes, directed_motifs):
    """
Retuns a list of all possible single component graphs given some number of nodes. This function
is meant primarily as a helper function to get_subgraphs, but can be used for other user purposes.

Parameters
----------
num_nodes : The number of nodes in the complete graph

directed_motifs : A boolean designating whether the motifs should be directed.

Returns
----------
graphs : A list of all possible single component graphs given some number of nodes
"""
    graphs = list()
    complete = nx.complete_graph(num_nodes)  # Start with complete graph
    complete_copy = copy.deepcopy(complete)
    if directed_motifs:
        complete = complete.to_directed()
        complete_copy = complete_copy.to_directed()

    edges = complete.edges()
    while complete.number_of_edges() > 0:
        # Iteratively remove edges, and capture single component subgraphs
        e = edges.pop()
        complete.remove_edge(e[0], e[1])
        if directed_motifs:
            if nx.number_weakly_connected_components(complete) == 1:
                graphs.append(copy.deepcopy(complete))
        else:
            if nx.number_connected_components(complete) == 1:
                graphs.append(copy.deepcopy(complete))
    graphs.append(complete_copy)
    # Add recursively produced undirected graphs to directed set
    #    if directed_motifs:
    #        graphs_edges=map(lambda g: g.edges(), graphs)
    #        graphs_edges.sort()
    #        for i in complete_ud:
    #            ud_edges=i.edges()
    #            ud_edges.sort()
    #            if ud_edges not in graphs_edges:
    #                graphs.append(i.to_directed())
    return graphs


def testdict():
    a = dict()
    for count in range(4):
        a[count] = count


def motif_communities(H, th_motif=7):
    count = 0
    test_CC = dict()
    # refined_components=partition(H)
    refined_components = [H.nodes()]
    for allpoints in refined_components:
        if len(allpoints) > 5:
            print("count", count, "# nodes --------------------->", len(allpoints))
            print("count", count, "# edges --------------------->", len(H.subgraph(allpoints).edges()))
            print("price to pay", len(H.subgraph(allpoints).edges()) * 1. / len(
                nx.complete_graph(len(allpoints)).edges()))
            test_CC[count] = mp2.motif_main(th_motif, H.subgraph(allpoints), count)
            count += 1
        else:
            print("---- !SKIPPED! ----- ")
    return test_CC


def clustering_hist(test_CC, ax1=0, ax2=0, show_plot=True):
    cl_th = 0.1
    graph_comm = 0

    kol = {}
    if show_plot:

        if ax1 == 0:
            fig = plt.figure(figsize=(8, 8))
            ax1 = fig.add_axes([0.1, 0.1, 0.85, 0.15])
            ax2 = fig.add_axes([0.1, 0.3, 0.85, 0.6])

    n_groups = len(test_CC[list(test_CC.keys())[0]])
    bar_width = 0.8 / len(test_CC) / 2
    index_nodes = np.ones(20)
    index = np.ones(20)
    for i in sorted(test_CC.keys()):
        normalizer = np.ones(20)
        for ele in test_CC[i]:
            normalizer[len(ele[1].nodes())] += ele[2]
            # print i,"--------------------->",i,len(test_CC[i]),n_groups
        index = np.array([test_CC[i][k][0] for k in range(n_groups)])
        # test_CC[i][2][2]=test_CC[i][2][2]-test_CC[i][1][2]

        kol[i] = [count[2] * 1. / normalizer[len(count[1].nodes())] for count in test_CC[i]]

        tmp = np.array([count[2] * 1. / normalizer[len(count[1].nodes())] for count in test_CC[i]])
        index_nodes = np.array([len(count[1].nodes()) for count in test_CC[i]])
        # print max(index_nodes)
        clustered_hist = np.zeros(max(index_nodes))
        path_hist = np.zeros(max(index_nodes))
        for j in range(max(index_nodes)):
            tmp2 = tmp[index_nodes == j]
            clu_co = np.array([nx.average_clustering(count[1]) for count in test_CC[i] if index_nodes[count[0]] == j])
            clustered_hist[j] = sum(tmp2[clu_co > cl_th])
            path_hist[j] = sum(tmp2[clu_co < cl_th])
            index = np.arange(len(clustered_hist)) + 1
        if show_plot:

            if isinstance(i, int):
                ax2.bar(index + (2 * i) * bar_width, clustered_hist, bar_width,
                        color='r', label='%s-th component' % i)
                ax2.bar(index + (2 * i + 1.) * bar_width, path_hist, bar_width,
                        color='b', label='%s-th component' % i)
            else:
                ax2.bar(index + 0.9 - (2 * graph_comm + 1) * bar_width, clustered_hist, bar_width,
                        color='r', label='all_graph%s' % graph_comm)
                ax2.bar(index + 0.9 - (2 * graph_comm) * bar_width, path_hist, bar_width,
                        color='b', label='all_graph%s' % graph_comm)
                graph_comm += 1

    ax2.set_xlim([2, max(index_nodes) + 2])
    if show_plot:
        kk = 0
        for i in range(n_groups):

            if nx.average_clustering(test_CC[0][i][1]) < cl_th:
                kk += 1
                tmppos = {}
                N = test_CC[0][i][1].nodes()
                NN = len(N)
                if NN > 6:
                    k = 0
                    for j in N:
                        tmppos[j] = [np.cos(2 * np.pi * k * 1. / NN) * .35 + (kk + 0.4),
                                     np.sin(2 * np.pi * k * 1. / NN) * .35 - 1]
                        k += 1
                    nx.draw_networkx_edges(test_CC[0][i][1], tmppos, color='k', lw=2, ax=ax1)
        # print kk
        ax1.set_xlabel('Motif ID')
        ax2.set_ylabel('Occurrences')
        ax2.set_xticks(index + 0.4)
        ax2.set_xticklabels(index)
        ax1.set_xticks([])
        ax1.set_yticks([])
