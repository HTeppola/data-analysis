import igraph as ig
import networkx as nx
import numpy as np
import copy
import pylab as plt
import time


def graph_partition(H, showplot=0, itera=10):
    import copy
    h = copy.deepcopy(H)
    G, mapping = nx_to_igraph(h, pos_flag=1)
    wc = G.community_infomap(trials=itera)
    groups = partition(wc)
    i_mapping = dict(list(zip(list(mapping.values()), list(mapping.keys()))))
    for ele in groups:
        for i in range(len(ele)):
            ele[i] = i_mapping[ele[i]]
    if showplot:
        plot_graph(wc)

    return groups


def relabel(G):
    mapping = {}
    j = 0
    for i in G.nodes():
        mapping[i] = j
        j += 1
    return nx.relabel_nodes(G, mapping), mapping


def my_product(G1, G2):

    G12 = nx.Graph()


    nodes = [(u, v) for u in G1.nodes() for v in G2.nodes()]
    mapping = {}
    j = 0
    for i in nodes:
        mapping[i] = j
        j += 1

    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)

    G12 = nx.cartesian_product(G12, nx.complete_graph(2))

    return G12


def modular_product(g1, g2):
    G1 = g1.to_undirected()
    G2 = g2.to_undirected()
    G12 = nx.Graph()


    nodes = [(u, v) for u in G1.nodes() for v in G2.nodes()]
    mapping = {}
    j = 0
    for i in nodes:
        mapping[i] = j
        j += 1
    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[0]), (edge1[0], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[0], edge2[1]), (edge1[1], edge2[0])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[1]), (edge1[0], edge2[0])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)

    #    UNCOMMENT IF TRUE DETECTION
    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[0]), (edge1[0], edge2[1])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[0], edge2[1]), (edge1[1], edge2[0])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[1]), (edge1[0], edge2[0])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)

    return G12


def nx_to_igraph(H, pos_flag=1, direction=True):
    G = copy.deepcopy(H)
    G.remove_nodes_from(nx.isolates(G))
    G, mapping = relabel(G)
    g = ig.Graph(directed=direction)
    g.add_vertices(len(G.nodes()))
    g.add_edges(G.edges())
    if pos_flag:
        pos = nx.get_node_attributes(G, 'pos')
        if pos == {}:
            pos = nx.spring_layout(G)

        g.vs['x'] = [ele[0] for ele in list(pos.values())]
        g.vs['y'] = [ele[1] for ele in list(pos.values())]
    return g, mapping


def partition(wc, flag=0):
    try:
        tmp = wc.as_clustering()
    except:
        tmp = wc
    lst = np.array(tmp.membership)
    groups = []
    for i in np.unique(lst):
        groups.append(np.where(lst == i)[0].tolist())
    if flag:
        np.save(open('/home/dlonardoni/graph_partition.npy', 'w'), groups)
    return groups


def pairs(lst):
    groups = dict()
    for i in np.unique(lst):
        groups[tuple(np.where(lst == i)[0])] = i
    return groups


def plot_graph(wc):
    try:
        tmp = wc.as_clustering()
    except:
        tmp = wc
    comm = np.array(tmp.membership)
    print(tmp.summary())
    fig = ig.plot(tmp, vertex_size=5, mark_groups=pairs(comm))
    return fig


def all_common(G, groups, pivot=-1):
    import time

    all_k1 = []
    all_k2 = []

    GU = G
    if pivot < 0:
        for i in range(1, len(groups)):
            for j in range(i + 1, len(groups)):
                print("group\t %s\t vs\t %s:" % (i, j))
                print("nodes\t %s\t vs\t %s:" % (len(groups[-i]), len(groups[-j])))
                if len(groups[-i]) > 3 and len(groups[-j]) > 3:
                    start = time.time()
                    # local_k1,local_k2=find_common_subgraph(G.subgraph(groups[-i]),GU.subgraph(groups[-j]))
                    local_k1, local_k2 = find_common2(G.subgraph(groups[-i]), GU.subgraph(groups[-j]))
                    # all_K=add_K(all_K,local_k1,GU)
                    all_k1.append(local_k1)
                    all_k2.append(local_k2)
                    print("----> time elapsed: %f len(all_k1):%d" % (time.time() - start, len(local_k1)))
    else:
        i = pivot
        for j in range(i + 1, len(groups)):
            print("group\t %s\t vs\t %s:" % (i, j))
            print("nodes\t %s\t vs\t %s:" % (len(groups[-i]), len(groups[-j])))
            if len(groups[-i]) > 3 and len(groups[-j]) > 3:
                start = time.time()
                # local_k1,local_k2=find_common_subgraph(G.subgraph(groups[-i]),GU.subgraph(groups[-j]))
                local_k1, local_k2 = find_common2(G.subgraph(groups[-i]), GU.subgraph(groups[-j]))
                # all_K=add_K(all_K,local_k1,GU)
                all_k1.append(local_k1)
                all_k2.append(local_k2)
                print("----> time elapsed: %d len(all_k1):%d" % (time.time() - start, len(local_k1)))

    return all_k1, all_k2


def add_K(all_K, local_K, G2, flag=0):
    """Renumber the values of the dictionary from 0 to n
    """

    for k in local_K:
        new = 1
        for key in list(all_K.keys()):
            graph, occ = all_K[key]

            if nx.is_isomorphic(G2.subgraph(k), G2.subgraph(graph)):
                all_K[key][1] += 1
                new = 0
        if new or flag:
            all_K[len(list(all_K.keys())) + 1] = [k, 1]

    return all_K


def factors(n, result):
    for i in range(2, n + 1):  # test all integers between 2 and n
        if n % i == 0:
            result.append(i)
            result = factors(n / i, result)
            return result
    result.append(1)
    result = np.array(result)
    if len(result) == 3 and np.prod(result) > 4:
        return factors(np.prod(result) + 1, [1])
    else:
        return result


def common_plot(all_K, G, flag=0):
    plt.figure()
    keys = list(all_K.keys())
    if flag:
        keys = [key for key in list(all_K) if all_K[key][1] > 10]
    fac = factors(len(keys), result=[1])

    row = np.prod(fac[np.arange(0, len(fac), 2)])
    col = np.prod(fac[np.arange(1, len(fac), 2)])
    i = 0
    for key in keys:
        i += 1
        plt.subplot(row, col, i)
        ax = plt.gca()
        ax.set_xticks([])
        ax.set_yticks([])

        K = G.subgraph(all_K[key][0])
        #        try:
        #
        #        except:
        pos = nx.circular_layout(K)
        #        pos=nx.get_node_attributes(K,'pos')
        nx.draw_networkx_edges(K, pos=pos)
        #        plt.xlim([0,1])
        #        plt.ylim([0,1])
        plt.xlabel("# of found: %s" % all_K[key][1])


def find_common(G1, G2):
    #    G12=nx.cartesian_product(G1,G2)
    G12 = nx.Graph()
    start = time.time()
    #    print G1.size(),len(G1.nodes()),G2.size(),len(G2.nodes())
    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    #    edge2add=[((edge1[0],edge2[0]),(edge1[1],edge2[1])) for edge1 in nx.complement(G1).edges() for edge2 in nx.complement(G2).edges() ]
    #    G12.add_edges_from(edge2add)

    print("G12 created", time.time() - start)
    C12, mapping = nx_to_igraph(G12, pos_flag=0, direction=False)

    print("C12 created", time.time() - start)

    i_mapping = dict(list(zip(list(mapping.values()), list(mapping.keys()))))

    cliques = C12.largest_cliques()

    del C12
    print("Cliques!", len(cliques), time.time() - start)
    N_cli = 0
    all_maps = []
    for tmp in cliques:
        map2plot = sorted([i_mapping[ele][1] for ele in tmp])

        if len(np.unique(map2plot)) > 2:

            if map2plot not in all_maps:
                flag = 1

                for ele in all_maps:
                    if nx.is_isomorphic(G2.subgraph(ele), G2.subgraph(map2plot)):
                        flag = 0
                        break
                if flag or len(all_maps) == 0:
                    N_cli += 1
                    all_maps.append(map2plot)
                    showplot = 0
                    if showplot:
                        plt.figure()
                        plt.suptitle(
                            'Are Isomorphic: %s' % nx.is_isomorphic(G1.subgraph([i_mapping[ele][0] for ele in tmp]),
                                                                    G2.subgraph(map2plot)))
                        plt.subplot(121)
                        nx.draw(G1.subgraph([i_mapping[ele][0] for ele in tmp]))
                        plt.subplot(122)
                        nx.draw(G2.subgraph([i_mapping[ele][1] for ele in tmp]))

    return all_maps


def find_common2(G1, G2):
    #    G12=nx.cartesian_product(G1,G2)
    G12 = nx.Graph()
    start = time.time()
    #    print G1.size(),len(G1.nodes()),G2.size(),len(G2.nodes())
    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[0]), (edge1[0], edge2[1])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[0], edge2[1]), (edge1[1], edge2[0])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[1]), (edge1[0], edge2[0])) for edge1 in G1.edges() for edge2 in G2.edges()]
    G12.add_edges_from(edge2add)
    #    return G12
    #    UNCOMMENT IF TRUE DETECTION
    edge2add = [((edge1[0], edge2[0]), (edge1[1], edge2[1])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[0]), (edge1[0], edge2[1])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[0], edge2[1]), (edge1[1], edge2[0])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)
    edge2add = [((edge1[1], edge2[1]), (edge1[0], edge2[0])) for edge1 in nx.complement(G1).edges() for edge2 in
                nx.complement(G2).edges()]
    G12.add_edges_from(edge2add)

    print("G12 created", time.time() - start)
    C12, mapping = nx_to_igraph(G12, pos_flag=0, direction=True)
    del G12

    print("C12 created", time.time() - start)

    i_mapping = dict(list(zip(list(mapping.values()), list(mapping.keys()))))
    tstart = time.time()
    cliques = C12.largest_cliques()
    print(time.time() - tstart)
    del C12
    print("Cliques!", len(cliques), time.time() - start)
    N_cli = 0
    all_maps_G1 = set()
    all_maps_G2 = set()
    tl = len(all_maps_G2) + len(all_maps_G1)
    for tmp in cliques:
        map_G1 = sorted([i_mapping[ele][0] for ele in tmp])
        map_G2 = sorted([i_mapping[ele][1] for ele in tmp])
        if len(np.unique(map_G1)) > 2 and len(np.unique(map_G2)) > 2:
            tmp_1 = G1.subgraph(map_G1)
            tmp_2 = G2.subgraph(map_G2)
            if nx.is_isomorphic(tmp_1, tmp_2):
                # print 'Are Directed Isomorphic: True'
                all_maps_G1.add(tuple(map_G1))
                all_maps_G2.add(tuple(map_G2))
                # else:
                # print 'Are Undirected Isomorphic: %s'%nx.is_isomorphic(tmp_1.to_undirected(),tmp_2.to_undirected())
            if len(all_maps_G2) + len(all_maps_G1) > tl:
                tl = len(all_maps_G2) + len(all_maps_G1)
                N_cli += 1
                showplot = 0
                if showplot:
                    plt.figure()
                    tmp_1 = G1.subgraph(map_G1)
                    tmp_2 = G2.subgraph(map_G2)
                    plt.suptitle('Are Directed Isomorphic: %s\nAre Undirected Isomorphic: %s' % (nx.is_isomorphic(tmp_1, tmp_2), nx.is_isomorphic(tmp_1.to_undirected(), tmp_2.to_undirected())))

                    plt.subplot(121)
                    nx.draw(tmp_1)
                    plt.subplot(122)
                    nx.draw(tmp_2)
                    print(" map_G1", map_G1)
                    print(" map_G2", map_G2)

    def track2(path1, path2):

        for up in [ele for i in path1 for ele in G1.to_undirected().neighbors(i)]:
            if up not in path1:
                path1.extend([up])
                flag = 1
                for vp in [ele for j in path2 for ele in G2.to_undirected().neighbors(j)]:
                    if vp not in path2:
                        path2.extend([vp])
                        # print sorted(path2),vp,vp in path2
                        if len(path2) < 3:
                            track2(path1, path2)
                        else:
                            if nx.is_isomorphic(G1.subgraph(path1), G2.subgraph(path2)):
                                flag = 0
                                if len(path2) > 5:
                                    all_path1.append(tuple(path1))
                                    all_path2.append(tuple(path2))
                                track2(path1, path2)
                            else:
                                path2.remove(vp)
                if flag:
                    path1.remove(up)

    all_path1 = []
    all_path2 = []
    i_c = 0
    for path1, path2 in zip(all_maps_G1, all_maps_G2):
        print(i_c, len(all_maps_G1))
        i_c += 1
        if nx.is_isomorphic(G1.subgraph(path1), G2.subgraph(path2)):
            # track2(list(path1),list(path2))
            all_path1.append(path1)
            all_path2.append(path2)
    return all_path1, all_path2


def find_common_subgraph(G1, G2):
    n1xn2 = [(u, v) for u in G1.nodes() for v in G2.nodes()]
    all_path1 = []
    all_path2 = []
    import time

    def track1(path1, path2, g1, g2):

        for up in [ele for ele in g1.to_undirected().neighbors(path1[-1]) if ele not in path1]:
            path1.extend([up])
            for vp in [ele for ele in g2.to_undirected().neighbors(path2[-1]) if ele not in path2]:
                path2.extend([vp])
                if len(path2) < 3:
                    track1(path1, path2, g1, g2)
                else:
                    if nx.is_isomorphic(G1.subgraph(path1), G2.subgraph(path2)):
                        if len(path2) > 5:
                            all_path1.append(tuple(sorted(path1)))
                            all_path2.append(tuple(sorted(path2)))
                        track1(path1, path2, g1, g2)
                    else:
                        a = path2.pop()
                        try:
                            g2.remove_node(a)
                        except:
                            pass

            b = path1.pop()
            try:
                g1.remove_node(b)
            except:
                pass

    i = 0
    start = time.time()
    for ele in n1xn2:
        i += 1

        path1 = [ele[0]]
        path2 = [ele[1]]
        track1(path1, path2, copy.deepcopy(G1), copy.deepcopy(G2))
    print(i, len(n1xn2), time.time() - start)
    return all_path1, all_path2


def unify(all_K):
    new = dict()
    pivot = np.arange(40) * np.arange(-1, 39) / 2
    N_group = np.argmax(pivot == len(all_K[0]))

    for ii in [1, 2]:

        print(N_group)
        ele = all_K[ii - 1]
        count = 0
        j = 0
        i = 0
        for items in ele:

            j += 1

            if (i, j) not in list(new):
                new[(i, j, ii)] = dict()
            for item in items:

                if tuple(sorted(item)) not in list(new[(i, j, ii)]):
                    new[(i, j, ii)][tuple(sorted(item))] = 1
                else:
                    new[(i, j, ii)][tuple(sorted(item))] += 1
            count += 1
            if j >= N_group - 1:
                i += 1
                j = i
        N_group = np.argmax(pivot == len(all_K[0]))
    return new


def tuple_plot(all_K, G):
    N = np.max([ele[1] for ele in list(all_K)]) + 1
    color = ['r', 'b', 'r']
    fig, axs = plt.subplots(N, N)
    for keys in list(all_K):
        ax = axs[keys[0], keys[1]]
        if len(all_K[keys]) > 0:
            for key in all_K[keys]:
                ax.set_xticks([])
                ax.set_yticks([])

                K = G.subgraph(key)
                pos = nx.get_node_attributes(K, 'pos')
                #            pos=nx.circular_layout(K)
                plt.gca()
                nx.draw_networkx_edges(K, pos=pos, ax=ax, edge_color=color[keys[2]])

                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
        else:
            ax.set_xticks([])
            ax.set_yticks([])


def test(G1):
    def track(path1):
        for up in [ele for i in path1 for ele in G1.to_undirected().neighbors(i) if ele not in path1]:
            path1.extend([up])
            print(path1)

    for ele in G1.nodes():
        path1 = [ele]
        track(path1)
