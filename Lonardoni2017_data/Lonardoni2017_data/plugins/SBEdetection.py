import pylab as plt
import numpy as np
import networkx as nx

plt.ion()


def _convertClist(list2xN):
    CellList = []
    for ele in xrange(4097):
        CellList.append({'ID': ele, 'spikes': list2xN[1][list2xN[0] == ele]})
    return CellList


def distance(idA, idB):
    x = abs(idA % 64 - idB % 64)
    y = abs(int(idA) / 64 - int(idB) / 64)
    return (x + y)


def extract_subg(G, cH):
    to_extract = [node for node in G.nodes() if node[0] == cH]
    tmp = []
    for node in to_extract:
        tmp.extend(G.neighbors(node))
    to_extract.extend(tmp)
    return G.subgraph(to_extract)


import time


def compute_distance_graph2(spk, Tdistance, Tdist=10):
    G = nx.Graph()
    N = len(spk[0]) - 1
    tmp = np.argsort(spk[1])
    spk = spk[:, tmp]
    edge2add = []
    clist = _convertClist(spk)
    for i in np.unique(spk[0]):
        tsA = clist[int(i)]['spikes']
        for j in Tdistance[i]:
            tsB = clist[int(j)]['spikes']
            edge2add.extend([((i, A), (j, B)) for A in tsA for B in tsB[(tsB > A) & (tsB - A < Tdist)]])

    G.add_edges_from(edge2add)
    deg = G.degree()
    G.remove_nodes_from([node for node in G.nodes() if deg[node] <= 3])
    subg = list(nx.connected_components(G))[0]
    return subg  # ,nx.connected_component_subgraphs(G)[0]


def compute_distance_graph(spk, Tdistance, Tdist):
    G = nx.Graph()
    N = len(spk[0]) - 1
    tmp = np.argsort(spk[1])
    spk = spk[:, tmp]
    edge2add = []
    for i in xrange(N):
        candidate = np.searchsorted(spk[1][i:], Tdist + spk[1][i]) + i
        Neigh = Tdistance[spk[0][i]]
        i2add = [j for j in xrange(i + 1, candidate) if spk[0][j] in Neigh]
        #            i2add=[ele for ele in spk.T[i+1:candidate] if ele[0] in Neigh]
        edge2add.extend([((spk[0][i], spk[1][i]), (spk[0][j], spk[1][j])) for j in i2add])
    G.add_edges_from(edge2add)

    subg = nx.connected_components(G)[0]
    return subg


def _detect_SMART(spikes, Tbin=10, Tin=0, N_ACTIVE=.05, DIST=5, Tdist=10, MinWind=30):
    if N_ACTIVE > 1:
        N_ACTIVE = N_ACTIVE * 0.01
    G_list = []

    reorder = np.argsort(spikes[1])
    spikes = spikes[:, reorder]
    Tend = np.max(spikes[1])
    rate = np.histogram(spikes[1], bins=np.arange(Tin, Tend, Tbin))[0]

    # Extra parameter definition

    N = len(np.unique(spikes[0]))

    brtmp = []
    burst = []
    ibi_dist = []

    ID = 0
    end = 0
    Smart_Burst = []
    old = 0
    check = np.where(rate > N * N_ACTIVE)[0]
    Tdistance = {}
    IDs = set(np.unique(spikes[0]))
    for i in np.unique(spikes[0]):
        Tdistance[i] = set([j for j in IDs.intersection(set(range(int(i - 64 * DIST), int(i + 64 * DIST + 1)))) if
                            distance(i, j) < DIST])
    print "Tdistance DONE"
    kk = 0
    for pivot in check:
        kk += 1
        if pivot > old:
            #            plt.figure()
            start = time.time()
            to_search = np.where((spikes[1] > max(pivot * Tbin - 50, old * Tbin)) & (spikes[1] < pivot * Tbin + 150))[0]
            spkInBurst = spikes[:, to_search]
            nodes = set([(ele[0], ele[1]) for ele in spkInBurst.T])

            spkInBurst = np.asarray(spkInBurst)
            G = compute_distance_graph2(spkInBurst, Tdistance, Tdist)

            brtmp = []
            brtmp.append([ele[0] for ele in G])
            brtmp.append([ele[1] for ele in G])
            if np.max(brtmp[1]) - np.min(brtmp[1]) > MinWind:
                if len(np.unique(brtmp[0])) > N * N_ACTIVE:
                    Smart_Burst.append(brtmp)
                    #                G_list.append(G2)
            old = np.max(brtmp[1]) / Tbin

            ID += 1
            print "Time: %d, SPK: %d, iteration: %06d of %06d -> %d" % (
            pivot * Tbin, len(spkInBurst[0]), kk, len(check), time.time() - start)
            #        sys.stdout.write("\rTime: %d, SPK: %06d, iteration: %06d of %06d" % (old*Tbin,len(brtmp[0]),kk,len(check)))
            #        sys.stdout.flush()
    return Smart_Burst, G_list

