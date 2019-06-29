from multiprocessing import Process, Queue
import numpy as np
import tools.motifs as m
import pylab as plt
import networkx as nx
import tools.Glist as gl
import igraph_tools as it
import sys
import time

sys.setrecursionlimit(50000)
from multiPP_ST_ref import runFunctionsInParallel


def worker(work_queue, done_queue, base, w):
    while not work_queue.empty():
        nmotif = motif_match(work_queue.get(), base, w, 0)
        done_queue.put((nmotif[0], 0, nmotif[2]))
    # done_queue.put("%s not failed, counts %s"%(nmotif[0],nmotif[2]))

    return True


def worker2(i, motif, base):
    nmotif = motif_match(motif, base, i, 0)
    return nmotif[0], motif[1], nmotif[2]


#        done_queue.put("%s not failed, counts %s"%(nmotif[0],nmotif[2]))

def motif_match(motif, base, w, base_direction):
    start = time.time()
    index = motif[0]
    useNX = 0
    if index == 1:
        # print index,"of",len(motif_counts)
        A = nx.adj_matrix(base)
        A3 = A.dot(A).dot(A)
        triangles = 0
        for i in xrange(A3.shape[0]):
            triangles += A3[i, i]
        count = triangles / 6
        GM = nx.isomorphism.GraphMatcher(motif[1], motif[1])
        aicount = 0
        for _ in GM.isomorphisms_iter():
            aicount += 1

    else:
        # base=nx.line_graph(base)
        # print index,"of",len(motif_counts),
        # Calculate subgraph ismorphisms
        if useNX:
            gl.nmin = 1
            bm = nx.relabel.relabel_nodes(base, {key: (0, key) for key in base.nodes()})
            mm = nx.relabel.relabel_nodes(motif[1], {key: (0, -key - 1) for key in motif[1].nodes()})
            count = gl.clique_removal3(bm, mm)
            GM = nx.isomorphism.GraphMatcher(motif[1], motif[1])
            aicount = 0
            for _ in GM.isomorphisms_iter():
                aicount += 1

        else 1:
            #            a=it.nx_to_igraph(nx.line_graph(motif[1]),direction=0)[0]
            #            b=it.nx_to_igraph(nx.line_graph(base),direction=0)[0]
            a = it.nx_to_igraph(motif[1], direction=base_direction)[0]
            b = it.nx_to_igraph(base, direction=base_direction)[0]

            icount = b.count_subisomorphisms_vf2(a)
            aicount = a.count_isomorphisms_vf2(a)

            count = icount * 1. / aicount
            print "PROCESS:", w, "\t T-index:", len(motif[1]), motif[1].size(), "elapsed: %0ds" % (
                time.time() - start), '\t :', aicount, icount / aicount

    return index, motif[1], count, aicount


def motif_main_IG(tau, base, count, flag=0):
    start = time.time()
    base_direction = False
    mymotifs = True

    if mymotifs:

        motif_counts = m.open_motifs(tau)

    else:
        motif_counts = m.get_motifs(tau, base_direction)

    baseIG = it.nx_to_igraph(base, direction=base_direction)[0]
    task2exec = [[motif_match_IG, [motif[0], it.nx_to_igraph(motif[1], direction=base_direction)[0], baseIG, motif[0]]]
                 for motif in motif_counts[::-1]]

    sRCC = runFunctionsInParallel(listOf_FuncAndArgLists=task2exec, maxAtOnce=40, SZORDAN=0, offsetsSeconds=.5,
                                  sleepTime=np.log(len(base)))

    tmp_test_CC = sorted(sRCC[1], key=lambda x: x[0])
    test_CC = [(ele[0], motif_counts[ele[0]][1], ele[1], base) for ele in tmp_test_CC]
    if flag:
        print "motif computed in %ds" % (time.time() - start)
    return test_CC, count


def motif_match_IG(index, motif2test, base, w, reportInfo=0):
    start = time.time()

    icount = base.count_subisomorphisms_vf2(motif2test)
    aicount = motif2test.count_isomorphisms_vf2(motif2test)

    count = icount * 1. / aicount
    if reportInfo:
        print "PROCESS:", w, "\t T-index:", motif2test.ecount(), motif2test.vcount(), "elapsed: %0ds" % (
            time.time() - start), '\t :', aicount, icount / aicount

    return index, count, aicount


def main():
    import time
    s_times = []
    p_times = []
    errors = []
    test_nodes = np.arange(5, 100, 5)
    for N_node in test_nodes:
        base_direction = False
        mymotifs = True


        graph_list = [nx.erdos_renyi_graph(N_node, 5 * 1. / N_node)]
        init_t = time.time()
        for base in graph_list:
            tau = 6

            if mymotifs:
                # motif_counts=get_mymotifs(tau,base_direction)
                motif_counts = m.open_motifs(tau)
            else:
                motif_counts = m.get_motifs(tau, base_direction)
            workers = 40
            work_queue = Queue()
            done_queue = Queue()
            processes = []

            for motif in motif_counts:
                work_queue.put(motif)

            for w in xrange(workers):
                p = Process(target=worker, args=(work_queue, done_queue, base, w))
                p.start()
                processes.append(p)

            for p in processes:
                p.join()

            test_CC = {len(base.nodes()) * 2: []}
            # print "----------results-----------"
            while not done_queue.empty():
                res = done_queue.get()

                test_CC[len(base.nodes()) * 2].append(res)
                # print test_CC

        ptime = time.time() - init_t
        stime = np.nan

        checker = np.zeros((len(test_CC[test_CC.keys()[0]]), 4))
        for ele in test_CC.keys():
            for item in test_CC[ele]:
                checker[item[0], 0] = item[0]
                if int(ele) % 2 == 0:
                    checker[item[0], 1] = item[2]
                else:
                    checker[item[0], 2] = item[2]
        error = -1
        for i in xrange(len(test_CC[test_CC.keys()[0]])):
            checker[i, 3] = checker[i, 2] - checker[i, 1]
            if not checker[i, 3] == 0:
                error += 1
        # print checker
        print "G serial elapsed:", stime
        print "G parallel elapsed:", ptime
        s_times.append(stime)
        p_times.append(ptime)
        errors.append(error)
        # print "mismatches:",error
    # m.motif_hist(test_CC)
    plt.figure()
    plt.plot(test_nodes, s_times, lw=2, label='serial')
    plt.plot(test_nodes, p_times, 'r', lw=2, label='parallel')
    plt.legend()
    plt.figure()
    plt.plot(test_nodes, errors, lw=2)
    return test_nodes, s_times, p_times


def run(motif_counts, base):
    workers = 40
    work_queue = Queue()
    done_queue = Queue()
    processes = []
    test_CC = []
    for motif in motif_counts:
        work_queue.put(motif)

    for w in xrange(workers):
        p = Process(target=worker, args=(work_queue, done_queue, base, w))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    while not done_queue.empty():
        res = done_queue.get()
        test_CC.append(res)
    return test_CC


def motif_main(tau, base, count):


    base_direction = False
    mymotifs = True

    if mymotifs:
        motif_counts = m.open_motifs(tau)

    else:
        motif_counts = m.get_motifs(tau, base_direction)



    task2exec = [[motif_match, [motif, base, motif[0], base_direction]] for motif in motif_counts[::-1]]

    sRCC = runFunctionsInParallel(listOf_FuncAndArgLists=task2exec, maxAtOnce=10, SZORDAN=0, offsetsSeconds=.5,
                                  sleepTime=np.log(len(base)))


    test_CC = sorted(sRCC[1], key=lambda x: x[0])
    return test_CC, count


if __name__ == '__main__':
    main()
