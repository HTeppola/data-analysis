from multiprocessing import Process, Queue

from igraph_tools import *
import os

import os.path
import time

ISO = 0
GGs = []
groups = []
list_motif = []
s_info = 1
only_directed = 0
nmin = 4


def worker(work_queue, w, index, folder):
    #    global list_motif

    while not work_queue.empty():
        start = time.time()
        maps = work_queue.get()
        #        print maps

        #        print "JOB %03d ASSIGNED TO CPU %02d"%(maps[0],w)
        #        print " "
        if not os.path.isfile(folder + '/%s_%s.npy' % (maps[0], maps[1])):
            out = clique_removal2(maps[0], maps[1])
            if len(out[0]) > 0:
                np.save(open(folder + '/%02d_%02d.npy' % (maps[0], maps[1]), 'w'), [maps[0], maps[1], out[0]])
                #            print ">>>>>>>>>>>>>>>>> JOB %03d load: (%03d,%03d) -> %07d nodes"%(maps[0],out[1],out[2],out[1]*out[2])
                #        done_queue.put([maps[1],maps[2],out[0]])
        else:
            print "already processed: %s" % (folder + '/%s_%s.npy' % (maps[0], maps[1]))
        if work_queue.empty():
            print "CPU:%02d process %03d" % (w, maps[0]), "of", index, "\t elapsed: %0ds" % (time.time() - start)
            #        print 'process id:', os.getpid()
            #        print "queuesize %d"%(work_queue.qsize())
            #        print "=========================================="

    return True


def run_all_common(Glist, n_wor=25, gain=3):
    global GGs
    GGs = Glist
    folder = '/home/alonardoni/8192correlation/tmp_%d_%d/' % (len(Glist), n_wor)
    if not os.path.isdir(folder):
        os.mkdir(folder)

    workers = n_wor
    all_K = {}
    N = len(Glist)
    index = 0
    mapping = []
    for i in xrange(0, N):
        for j in xrange(0, N):
            if not i == j:
                mapping.append((i, j))
                #                print i,j,len(Glist[i]),len(Glist[j])
                index += 1
    mapping = mapping[::-1]
    counter = 0
    start = time.time()
    while len(mapping) > 0:
        time.sleep(1)

        print "#######################################"
        print "################# iteration", counter, "e", time.time() - start
        print "#######################################"
        start = time.time()
        work_queue = Queue()

        processes = []

        for k in xrange(int(workers * gain)):
            if len(mapping) == 0:
                break
            match = mapping.pop()
            #                print match
            work_queue.put((match[0], match[1]))
            if work_queue.full():
                print "QUEUE IS FULL"
                print "QUEUE IS FULL"
                print "QUEUE IS FULL"
                print "QUEUE IS FULL"

                time.sleep(20)

        jj = 0

        for w in xrange(workers):
            p = Process(target=worker, args=(work_queue, w, index, folder))
            jj += 1
            p.start()

            processes.append(p)

        for p in processes:
            p.join()

        for ele in os.listdir(folder):
            res = np.load(folder + ele)
            all_K[(res[0], res[1])] = res[2]
            os.remove(folder + ele)

        # for dq in done_queue:
        #            while not dq.empty():
        #                res=dq.get()

        for p in processes:
            del p
        gain *= 1
        counter += 1

    os.rmdir(folder)
    return all_K


def clique_removal2(group1, group2):
    """ Repeatedly remove cliques from the graph.

    Results in a `O(|V|/(\log |V|)^2)` approximation of maximum clique
    & independent set. Returns the largest independent set found, along
    with found maximal cliques.

    Parameters
    ----------
    group1 : NetworkX graph
        Undirected graph
    group2 : NetworkX graph
        Undirected graph

    Returns
    -------
    max_ind_cliques : (set, list) tuple
        Maximal independent set and list of maximal cliques (sets) in the graph.

    References
    ----------
    .. [1] Boppana, R., & Halld0rsson, M. M. (1992).
        Approximating maximum independent sets by excluding subgraphs.
        BIT Numerical Mathematics, 32(2), 180-196. Springer.
    """
    import time
    global nmin
    global GGs
    global groups
    global s_info
    noderemoved = set([])

    if group1.is_directed():
        G1 = group1.to_undirected().copy()
        G2 = group2.to_undirected().copy()
    else:
        G1 = group1.copy()
        G2 = group2.copy()

    nmin = 1
    l1 = len(G1.nodes())
    l2 = len(G2.nodes())


    mc = 0
    deg1 = G1.degree()
    deg2 = G2.degree()
    nnlist = set([(u, v) for u in G1.nodes() for v in G2.nodes() if not (u, v) in noderemoved if u[0] == v[0]])
    nnlist = [(node[0], node[1], deg1[node[0]] * deg2[node[1]]) for node in nnlist]
    stopat = len(nnlist) - 1
    start = time.time()
    c_i, i_i, n = ramseyGlist(G1, G2, ww=0, n_tocheck=nnlist)
    cliques = [c_i]
    isets = [i_i]
    iteration = 0
    counter = 0
    while stopat > len(noderemoved) and counter < 5000000:

        iteration += 1
        if iteration > 100000000:
            break
        # if s_info:
        #            print "iter:",iteration,"node removed:",len(noderemoved),"of:",stopat,len(c_i),"e: %d"%(time.time()-start)


        #        if len(i_i)>=1:
        #         print "----->",i_i,"<-----"
        #         print c_i
        #         print n.intersection(nnlist)
        #         print n
        # time.sleep(1)
        if c_i:
            #           print c_i
            cliques.append(c_i)
        # map1=[item[0] for item in list(c_i)]
        #            map2=[item[1] for item in list(c_i)]
        #            print nnlist
        #            print "clique",nx.is_isomorphic(GG.subgraph(map1).to_undirected(),GG.subgraph(map2).to_undirected())
        #            print c_i
        #            print i_i

        if i_i:
            isets.append(i_i)
        # map1=[item[0] for item in list(i_i)]
        #            map2=[item[1] for item in list(i_i)]
        #            print "indipendent",nx.is_isomorphic(GG.subgraph(map1).to_undirected(),GG.subgraph(map2).to_undirected())
        #            print " "
        i_i = set(i_i)  # [(node[0],node[1]) for node in i_i])
        c_i = set(c_i)  # [(node[0],node[1]) for node in c_i])
        map1 = [item[0] for item in i_i]

        if G1.subgraph(map1).size() <= nmin:
            counter += 1
        else:
            mc = max(counter, mc)
            counter = 0
        for ele in i_i.union(c_i):
            nnlist.remove(ele)
        # print "gain",len(i_i.union(c_i))-len(i_i),"i_i",len(i_i)
        #        print "counter:",counter
        #        g1=G1.subgraph([u[0] for u in nnlist])
        #        g2=G2.subgraph([u[1] for u in nnlist])
        #        print "DMAX",max([g1.degree(u)*g2.degree(v) for (u,v) in nnlist]),mc
        noderemoved = noderemoved.union(i_i.union(c_i))
        #        noderemoved=noderemoved.union(i_i)


        nnlist = [(node[0], node[1], deg1[node[0]] * deg2[node[1]]) for node in nnlist]
        #        if iteration>8:
        #         break
        if nnlist:
            c_i, i_i, n = ramseyGlist(G1, G2, ww=0, n_tocheck=nnlist)

    lenc = len(isets)
    filtered = set()
    isolates = 0

    cc = 0
    while len(isets) > 0:
        cc += 1
        #        print cc
        ele = isets.pop()
        #        print ele
        map1 = [item[0] for item in list(ele)]
        map2 = [item[1] for item in list(ele)]
        g1 = G1.subgraph(map1).copy()
        #        if len(map1)>nmin:
        #            tmp=[(u,v) for u,v in zip(map1,map2)]
        #            filtered.add(tuple(tmp))

        g1.remove_nodes_from(nx.isolates(G1))
        map1 = [item for item in g1.nodes()]
        g2 = G2.subgraph(map2).copy()
        g2.remove_nodes_from(nx.isolates(G2))
        map2 = [item for item in g2.nodes()]
        if len(map1) > nmin:
            tmp = [(u, v) for u, v in zip(map1, map2)]
            filtered.add(tuple(tmp))
        if ISO:
            g1.remove_nodes_from(nx.isolates(G1))
            g2 = G2.subgraph(map2).copy()
            g2.remove_nodes_from(nx.isolates(G2))
            if nx.is_isomorphic(g1, g2):
                isolates += 1
                #          if g1.size()>len(map1)-1:#>(nmin-1)*nx.number_connected_components(g1) or (len(map1)>thx-1):

                #    print "------------------------------------------ MC",mc,counter
                #    print "compression:", lenc,"->",len(filtered)
    if ISO:
        print "CHECK ISO:", not lenc - isolates, isolates, time.time() - start, lenc, isolates

    return filtered, l1, l2


def clique_removal3(group1, group2):
    """ Repeatedly remove cliques from the graph.

    Results in a `O(|V|/(\log |V|)^2)` approximation of maximum clique
    & independent set. Returns the largest independent set found, along
    with found maximal cliques.

    Parameters
    ----------
    group1 : NetworkX graph
        Undirected graph
    group2 : NetworkX graph
        Undirected graph

    Returns
    -------
    max_ind_cliques : (set, list) tuple
        Maximal independent set and list of maximal cliques (sets) in the graph.

    References
    ----------
    .. [1] Boppana, R., & Halld0rsson, M. M. (1992).
        Approximating maximum independent sets by excluding subgraphs.
        BIT Numerical Mathematics, 32(2), 180-196. Springer.
    """

    global nmin
    global GGs
    global groups
    global s_info
    noderemoved = set([])

    if group1.is_directed():
        G1 = group1.to_undirected().copy()
        G2 = group2.to_undirected().copy()
    else:
        G1 = group1.copy()
        G2 = group2.copy()

    nmin = 1

    lM = len(G2)
    eM = G2.size()
    mc = 0
    deg1 = G1.degree()
    deg2 = G2.degree()
    nnlist = set([(u, v) for u in G1.nodes() for v in G2.nodes() if not (u, v) in noderemoved if u[0] == v[0]])
    nnlist = [(node[0], node[1], deg1[node[0]] * deg2[node[1]]) for node in nnlist]
    stopat = len(nnlist) - 1

    c_i, i_i, n = ramseyGlist(G1, G2, ww=0, n_tocheck=nnlist)
    cliques = [c_i]
    isets = [i_i]

    counter = 0
    map1=[]
    while stopat > len(noderemoved):

        if c_i:
            cliques.append(c_i)

        if i_i:
            isets.append(i_i)

        i_i = set(i_i)  # [(node[0],node[1]) for node in i_i])
        c_i = set(c_i)  # [(node[0],node[1]) for node in c_i])
        map1 = [item[0] for item in i_i]

        if G1.subgraph(map1).size() <= nmin:
            counter += 1
        else:
            mc = max(counter, mc)
            counter = 0
        for ele in i_i.union(c_i):
            nnlist.remove(ele)

        noderemoved = noderemoved.union(i_i.union(c_i))
        nnlist = [(node[0], node[1], deg1[node[0]] * deg2[node[1]]) for node in nnlist]
        if nnlist:
            c_i, i_i, n = ramseyGlist(G1, G2, ww=0, n_tocheck=nnlist)


    filtered = []

    #    print "ISOMORPH?",lenc

    cc = 0
    while len(isets) > 0:
        print isets
        cc += 1
        ele = isets.pop()
        map2 = [item[1] for item in list(ele)]
        g2 = G2.subgraph(map2).copy()
        g2.remove_nodes_from(nx.isolates(G2))
        map2 = [item for item in g2.nodes()]
        if len(map2) == lM and g2.size() == eM:
            tmp = [(u, v) for u, v in zip(map1, map2)]
            filtered.append(tuple(tmp))
    return len(filtered)


def ramseyGlist(G1, G2, ww, n_tocheck, AN=set()):
    r"""Approximately computes the Ramsey number `R(2;s,t)` for graph.

    Parameters
    ----------
    G1 : NetworkX graph
        Undirected graph
    G2 : NetworkX graph
        Undirected graph
    ww : Counter for recursion
    n_tocheck : Nodes to check
    AN : Previous Match
    Returns
    -------
    max_pair : (set, set) tuple
        Maximum clique, Maximum independent set.
    """
    #    if ww==0:
    #        n_tocheck=[(node[0],node[1],G1.degree(node[0])*G2.degree(node[1])) for node in n_tocheck]
    ww += 1

    if len(n_tocheck) > 0:
        d = [u[2] for u in n_tocheck]
        pivot = d.index(max(d))

        node = list(n_tocheck)[pivot]
        neigh1 = set(G1.neighbors(node[0]))
        neigh2 = set(G2.neighbors(node[1]))
        N_neigh1 = set([u[0] for u in n_tocheck if not u[0] == node[0]]).difference(neigh1)
        N_neigh2 = set([u[1] for u in n_tocheck if not u[1] == node[1]]).difference(neigh2)

        C_NC = [(u, v, d) for (u, v, d) in n_tocheck if
                (u in neigh1) & (v in neigh2) & (u[0] == v[0])]  # for v in neigh2 if (u,v) in n_tocheck]
        #        C_NC=[(u,v) for u in neigh1 for v in neigh2 if (u,v) in n_tocheck]
        OLD = AN

        tmpN = set(C_NC)
        tmpN.add(node)
        AN = AN.union(tmpN)
        C_NC.extend([(u, v, d) for (u, v, d) in n_tocheck if (u in N_neigh1) & (v in N_neigh2) & (u[0] == v[0])])
        #        C_NC.extend([(u,v) for u in N_neigh1 for v in N_neigh2 if (u,v) in n_tocheck])

        C_NBRS = set(n_tocheck).difference(set(C_NC))  # [nnode for nnode in n_tocheck if nnode not in C_NC]

        C_NBRS.remove(node)

        c_1, i_1, a1 = ramseyGlist(G1, G2, ww, n_tocheck=list(C_NBRS), AN=AN)
        if OLD.intersection(tmpN) or ww == 1:
            c_2, i_2, a2 = ramseyGlist(G1, G2, ww, n_tocheck=C_NC, AN=AN)
            i_2.add(node)
            c_1.add(node)
        else:

            c_2 = set()
            i_2 = set()
            a2 = OLD

        if len(i_2) >= len(i_1):
            i_max = i_2
            a_max = a2
        else:
            i_max = i_1
            a_max = a1
        return max([c_2, c_1], key=len), i_max, a_max  # max([c_2, c_1],key=len),i_max,a_max
    else:
        return set([]), set([]), AN