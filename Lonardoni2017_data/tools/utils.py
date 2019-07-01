__author__ = 'dlonardoni'
#
#   Utilities and statistic procedure
#

import time
import networkx as nx
import pylab as plt
import numpy as np
import numpy.random as rnd
import scipy.stats
from scipy.spatial.distance import pdist, squareform

'''
                ******** Utility procedures ********
'''


def Save_Traces(CellList, groupname, sname):
    trace = []
    plt.figure()
    flags = 0
    for cell in CellList:
        try:
            diction = {'ID_syn': cell.ID_syn, 'ID': cell.num, 'x': cell.pos['x'],
                       'y': cell.pos['y'], 't': np.array(cell.record['time']),
                       'v': np.asarray(cell.record['voltage']),
                       'gAMPA': np.asarray(cell.record['gAMPA']),
                       'gGABA': np.asarray(cell.record['gGABA']),
                       'iAMPA': np.asarray(cell.record['iAMPA']),
                       'iGABA': np.asarray(cell.record['iGABA']),
                       }
            if cell.nmda:
                diction['gNMDA'] = np.asarray(cell.record['gNMDA'])
                diction['iNMDA'] = np.asarray(cell.record['iNMDA'])
            trace.append(diction)
            flags += 1
        except:
            pass
    np.save(open('Results' + groupname + '/' + sname + '_tt.npy', 'w'), trace)


'''
                ******** Graph procedures ********
'''


def CreateRandomGraph(N, p=0.1):
    """
        Procedure to create a random graph G(N,p):
        - a pair of node in the graph is connected with probability p independently
        - the edge orientation is chosen randomly

        Some statistics of the graph are printed:
            Average Shortest Path
            Clustering Coefficient

        If the constructed graph is disconnected return an error

        See "erdos_renyi()" function of Networkx package

        Input:
        N:      Number of nodes
        p:      Probability [0,1] of connection

        Output:
        G:      Generated graph
    """
    G = nx.erdos_renyi_graph(N, p, directed=True)
    for _ in G:
        pass
    H = G.to_undirected()
    try:
        n = nx.average_shortest_path_length(G)
        print("Average Shortest Path : ", n)
    except nx.NetworkXError:
        print('Average Shortest Path : Disconnected')
    print("Clustering Coefficient ", nx.average_clustering(H))
    return G


def random_geometric_gauss(N, dim=2, sigma=1, grad=0, torus=0):
    """
        Return a random geometric graph with gaussian distribution of connection as a function of distance.
        The nodes are selected randomly in a hypercube in the dimension provided.

        The random geometric graph model
        - places n nodes uniformly at random in the unit cube of dimension dim
        - connect two nodes `u,v` with an edge in such a way that `d(u,v)`
                  generates a unilateral gaussian distribution N(0,sigma)

        where `d(u,v)` is the Euclidean distance between node u and v.

        Input:
        N:          Number of nodes
        sigma:      Variance of the Gaussian
        dim :       optional, Dimension of graph

        Output:
        G           Constructed graph

        Position is added to each node
    """
    G = nx.Graph()
    G.name = "Random Geometric Graph"
    G.add_nodes_from(list(range(N)))

    #   sample node position uniformly
    for n in G:
        G.node[n]['pos'] = rnd.random(dim)
    if dim == 3:
        for n in range(len(G.nodes())):
            G.node[n]['pos'][2] = (1 - G.node[n]['pos'][2] ** grad) * .25
    nodes = G.nodes(data=True)
    #   create the connections
    dmax = 0
    i = 0
    s = .5
    prob = rnd.random(N * N / 2).tolist()
    while nodes:

        u, du = nodes.pop()
        print(u)
        pu = du['pos']
        for v, dv in nodes:
            i += 1
            pv = dv['pos']
            d = sum(((a - b) ** 2 for a, b in zip(pu, pv)))
            if dim == 3:
                dxy = sum(((a - b) ** 2 for a, b in zip(pu[:-1], pv[:-1])))
                dz = (pu[-1] - pv[-1]) ** 2
                d = (s * dxy + (1 - s) * dz) * 1. / s
            if torus:
                d = sum(((min(abs(a - b), 1 - abs(a - b))) ** 2 for a, b in zip(pu, pv)))
            if d < .5 ** 2:
                p = scipy.stats.chi2(1).cdf(d / sigma)

                if p <= prob.pop():
                    G.add_edge(u, v)
                    dmax = max(d, dmax)
    return G


def DirectGraph(G):
    """
        Unefficien procedure to convert an undirected graph into a directed one
        The direction of an edges is chosen randomly

        Input:
        G:          networkx undirected graph

        Output:
        H:          directed realization of G
    """
    import time
    start = time.time()
    pivot = rnd.random(len(G.edges()))
    H = G.to_directed()
    import copy
    GG = copy.deepcopy(G)
    edges = np.array(G.edges())
    H.remove_edges_from(edges[pivot > .5])

    GG.remove_edges_from(edges[pivot > .5])
    G2 = GG.to_directed()
    a1 = set(G2.edges())
    a2 = a1.difference(set(GG.edges()))
    H.remove_edges_from(a2)
    print("elapsed:", time.time() - start)

    return H


def CreateGauss(N, dim=2, sigma=1, torus=0):
    """
        Return a random geometric graph with gaussian distribution of connection as a function of distance.
        The nodes are selected randomly in a hypercube in the dimension provided.

        The random geometric graph model
        - places n nodes uniformly at random in the unit cube of dimension dim
        - connect two nodes `u,v` with an edge in such a way that `d(u,v)`
                   generates a unilateral gaussian distribution N(0,sigma)

        where `d(u,v)` is the Euclidean distance between node u and v.

        Some statistics of the graph are printed:
                Average Shortest Path
                Clustering Coefficient

            If the constructed graph is disconnected return an error

        Input:
        N:          Number of nodes
        sigma:      Variance of the Gaussian
        dim :       optional, Dimension of graph

        Output:
        G           Constructed graph

        Position is added to each node
    """

    G, edge_list = random_geometric_gauss2(N=N, dim=dim, sigma=sigma, torus=torus)
    print("done")
    #    nx.draw(G,pos=nx.get_node_attributes(G,'pos'))
    #    edge_list=G.edges()
    import time
    start = time.time()

    def shuffle(xVec):
        if np.random.random() > .5:
            xVec = xVec[::-1]
        return xVec

    edge_list.extend([(ele[1], ele[0]) for ele in edge_list])
    edge_list = [shuffle(edge) for edge in edge_list if edge[0] > edge[1]]

    H = G.to_directed()
    print(time.time() - start)
    H.remove_edges_from(H.edges())
    print(time.time() - start)
    H.add_edges_from(edge_list)
    print(time.time() - start)
    #    H=DirectGraph(G)
    #    try:
    #        n=nx.average_shortest_path_length(H)
    #        print "Average Shortest Path : ",n
    #    except nx.NetworkXError:
    #        print 'Average Shortest Path : Disconnected'
    #    print "Clustering Coefficient ",nx.average_clustering(G)

    # plt.figure()
    # nx.draw(H,pos=nx.get_node_attributes(H,'pos'))
    return H


def CreateRadiusGraph(N, r=0.1):
    """
        Return a radius graph embedeed in a square of side one in the plane.

        The radius graph model
        - places n nodes uniformly at random in the unit cube of dimension two
        - connect two nodes `u,v` with an edge if `d(u,v)` < r

        where `d(u,v)` is the Euclidean distance between node u and v.

        See "random_geometric_graph()" function of Networkx package

        Some statistics of the graph are printed:
                Average Shortest Path
                Clustering Coefficient

            If the constructed graph is disconnected return an error

        Input:
        N:          Number of nodes
        r:          Radius of interaction

        Output:
        G           Directed realization of the constructed graph

        Position is added to each node
    """

    G = nx.random_geometric_graph(N, r, 2)
    H = G.to_directed()
    listed = H.edges()
    k = 0
    for i in range(len(listed) / 2):
        while listed:
            if (listed[k][1], listed[k][0]) in H.edges():
                if rnd.random() < 0.5:
                    H.remove_edge(*listed[k])
                    break
                else:
                    H.remove_edge(listed[k][1], listed[k][0])
                    break
            k += 1
        listed = H.edges()

    try:
        n = nx.average_shortest_path_length(H)
        print("Average Shortest Path : ", n)
    except nx.NetworkXError:
        print('Average Shortest Path : Disconnected')
    print("Clustering Coefficient ", nx.average_clustering(G))
    return H


'''
                ******** Plot procedures ********
'''


def PlotRaster(LS, ts, reo=1):
    """
        Procedure to draw a raster plot of the activity
        A line divides the excitatory neurons (top)
        from the inhibitory ones (bottom)

        Input:
        LS:         List of neuron cells
        ts:         Plot up to time "ts"
        reo:        Reorder neurons
        Output:
        Raster plot

    """

    spike = []
    tmpe = []
    tmpi = []
    N = len(LS)

    #   retrive the spike of each cell splitting into excitatory and inhibitory
    for i in range(N):
        if LS[i].ID_syn.rfind('EXC'):
            tmpe.append(np.array(LS[i].record['spk']))
        else:
            tmpi.append(np.array(LS[i].record['spk']))
    if reo == 1:
        #   reordering of the vector
        for i in range(N):
            if i < len(tmpe):
                spike.append(tmpe[i])
            if i == len(tmpe):
                print(ts)
                spike.append([0, ts])
            if i > len(tmpe):
                spike.append(tmpi[i - len(tmpe)])
    else:
        CellPos = {}
        NCells = len(LS)
        for i in range(NCells):
            # check if it is possible to color differently exc/inh neurons ?!
            CellPos[i] = np.array(list(LS[i].pos.values()))
        index = sorted(list(range(NCells)), key=lambda k: CellPos[k][1])
        for i in range(NCells):
            spike.append(np.array(LS[index[i]].record['spk']))

    # generating the plot
    for i, spk in spike:
        plt.plot(spk, np.ones_like(spike) * i, 'ok')
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron ID')


def PlotVoltage(cell):
    """
        Plot voltage of a single cell "cell" vs time (works when recordAll=TRUE)
    """
    # check before if "recordAll" was TRUE
    t = np.asarray(cell.record['time']) * .001

    v = np.asarray(cell.record['voltage']) * .001
    plt.plot(t, v)
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (mV)')


#    if ShowPlot: plt.show()
#    if HoldPlot: plt.hold(1)

def PlotG(cell):
    """
        Plot synaptic conductance of a single cell "cell" vs time

    """

    t = np.asarray(cell.record['time']) * .001

    gAMPA = np.asarray(cell.record['gAMPA'])
    plt.plot(t, -gAMPA, 'orange', label='gAMPA')
    gGABA = np.asarray(cell.record['gGABA'])
    plt.plot(t, -gGABA, 'b', label='gGABA')
    gNMDA = np.asarray(cell.record['gNMDA'])
    plt.plot(t, -gNMDA, 'r', label='gNMDA')
    plt.ylabel('conductances (nS)')
    plt.xlabel('time (s)')
    plt.legend()


def PlotCurrent(cell):
    """
        Plot synaptic current of a single cell "cell" vs time
    """
    t = np.asarray(cell.record['time']) * .001

    iAMPA = np.asarray(cell.record['iAMPA']) * .001
    iGABA = np.asarray(cell.record['iGABA']) * .001
    iNMDA = np.asarray(cell.record['iNMDA']) * .001
    plt.plot(t, iAMPA, 'orange', lw=2, label='iAMPA')
    plt.plot(t, iGABA, 'b', lw=2, label='iGABA')
    plt.plot(t, iNMDA, 'r', lw=2, label='iNMDA')

    plt.xlabel('time (s)')
    plt.ylabel('currents (nA)')
    plt.legend()


def crop_graph(G, pt=.5):
    H = G.subgraph([i for i in G.nodes() if G.node[i]['pos'][0] < pt and G.node[i]['pos'][1] < pt])
    return nx.convert_node_labels_to_integers(H)


def draw_node(G, node=0):
    pos = nx.get_node_attributes(G, 'pos')
    nx.draw_networkx_edges(G, pos=pos, alpha=.05, edge_color='k')
    nx.draw_networkx_edges(G, pos=pos, alpha=.8, edge_color='r', edgelist=G.in_edges(node), width=2)
    nx.draw_networkx_edges(G, pos=pos, alpha=.8, edge_color='g', edgelist=G.out_edges(node), width=2)
    plt.xlim([0, np.max(list(pos.values()))])
    plt.ylim([0, np.max(list(pos.values()))])


def random_geometric_gauss2(N, dim=2, sigma=1, grad=0, torus=0):
    """
        Return a random geometric graph with gaussian distribution of connection as a function of distance.
        The nodes are selected randomly in a hypercube in the dimension provided.

        The random geometric graph model
        - places n nodes uniformly at random in the unit cube of dimension dim
        - connect two nodes `u,v` with an edge in such a way that `d(u,v)`
                  generates a unilateral gaussian distribution N(0,sigma)

        where `d(u,v)` is the Euclidean distance between node u and v.

        Input:
        N:          Number of nodes
        sigma:      Variance of the Gaussian
        dim :       optional, Dimension of graph

        Output:
        G           Constructed graph

        Position is added to each node
    """
    G = nx.Graph()
    G.name = "Random Geometric Graph"
    G.add_nodes_from(list(range(N)))

    #   sample node position uniformly
    for n in G:
        G.node[n]['pos'] = rnd.random(dim)
    if dim == 3:
        for n in range(len(G.nodes())):
            G.node[n]['pos'][2] = (1 - G.node[n]['pos'][2] ** grad) * .25
    start = time.time()
    pos = np.asarray(list(nx.get_node_attributes(G, 'pos').values()))

    if torus:
        Dmat = squareform(pdist(pos, metric=lambda x, y: (
            min(1 - abs(x[0] - y[0]) % 1, abs(x[0] - y[0])) ** 2 + min(1 - abs(x[1] - y[1]) % 1,
                                                                       abs(x[1] - y[1])) ** 2)))
    else:
        Dmat = squareform(pdist(pos, 'sqeuclidean'))
    print('%d s' % (time.time() - start))
    Dmat = scipy.stats.chi2(1).cdf(Dmat / sigma)
    Dmat += np.eye(len(Dmat)) * 10
    print('%d s' % (time.time() - start))
    Pmat = np.random.random(Dmat.shape)
    print('%d s' % (time.time() - start))

    edges = np.where(Dmat - Pmat < 0.0)
    print('%d s' % (time.time() - start))
    #    G.add_edges_from()

    return G, list(zip(*edges))
