import os
import time
import igraph_tools as it
import matplotlib
import pylab as plt
import numpy as np
import networkx as nx
import tools.concave_hull as cv
import tools.mp2 as mp2
import copy
import community
from scipy.stats import rankdata
import pickle
from matplotlib.collections import LineCollection


def reduce_community(com, pos, nMax=50):
    com = np.asarray(com)
    posx, posy = list(zip(*list(pos.values())))
    posx = np.asarray(posx)
    posy = np.asarray(posy)
    sset = np.arange(len(com))
    tset = sset[sset >= 0]
    print(posx.shape)
    print(posy.shape)
    print(np.mean(posx[com[tset]]), np.mean(posy[com[tset]]))
    while len(tset) > nMax:
        x = np.mean(posx[com[tset]])
        y = np.mean(posy[com[tset]])
        dist = [abs(posx[node] - x) + abs(posy[node] - y) for node in com[tset]]
        sset[np.where(tset[np.argmax(dist)] == sset)[0]] = -1
        tset = sset[sset >= 0]
    print(np.mean(posx[com[tset]]), np.mean(posy[com[tset]]))
    return np.asarray(com)[tset]


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


plt.rcParams['legend.loc'] = 'best'


def motif_hist2(test_CC):
    fig = plt.figure(figsize=(8, 8))
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.15])
    ax2 = fig.add_axes([0.1, 0.3, 0.8, 0.6])
    simpleaxis(ax1)
    simpleaxis(ax2)
    n_CC = len(test_CC)

    pivot = set()
    for key in list(test_CC.keys()):
        sum_edges = max(sum([count[2] for count in test_CC[key]]), 1)

        counts = np.array([count[2] * 1. / sum_edges for count in test_CC[key]])
        n2add = set(np.argsort(counts)[-30:]).difference(set(np.where(counts == 0)[0].tolist()))
        pivot = pivot.union(n2add)

    pivot = np.array(sorted(list(pivot)))
    n_groups = len(pivot)
    if n_groups > 0:

        index = np.arange(n_groups)
        bar_width = 0.8 / len(test_CC)
        for i, key in enumerate(test_CC.keys()):
            sum_edges = sum([count[2] for count in test_CC[key]])
            counts = np.array([count[2] * 1. / sum_edges for count in test_CC[key]])
            if str(key).rfind('SC') < 0:

                ax2.bar(index + (i - 1) * bar_width, counts[pivot], bar_width,
                        color=plt.get_cmap('Greys')(i * .9 / n_CC),
                        label='%s-th component' % (i + 1))
            else:

                ax2.bar(index + (i - 1) * bar_width, counts[pivot], bar_width,
                        color='r',
                        label=key)
        tmp = 0
        for i in pivot:
            tmppos = {}
            in4g = list(test_CC.keys())[0]
            N = test_CC[in4g][i][1].nodes()
            k = 0
            NN = len(N)
            for j in N:
                tmppos[j] = [np.cos(2 * np.pi * k * 1. / NN) * .35 + (tmp + 0.4),
                             np.sin(2 * np.pi * k * 1. / NN) * .35 - 1]
                k += 1
            nx.draw_networkx_edges(test_CC[in4g][i][1], tmppos, color='k', lw=2, ax=ax1)
            tmp += 1
        ax2.set_xticks(index + 0.3)
        ax2.set_xticklabels(index)
    ax1.set_xlim([-1, n_groups + 1])
    ax2.set_xlim([-1, n_groups + 1])
    ax1.set_xlabel('Motif ID')
    ax2.set_ylabel('Occurrences')

    ax2.legend()
    ax1.set_xticks([])
    ax1.set_yticks([])
    fig.set_size_inches(len(pivot) / 1.5, 8, forward=True)


def extract_functional_graph3(G, hlist, mapp, mapt, D=(), show_plot=False, flag_scatter=False):
    print(G.name)
    link = nx.Graph()

    if show_plot:
        plt.figure()
        plt.xlim([0, 1])
        plt.ylim([0, 1])
    distance = []
    off = []
    distance_in = []
    off_in = []
    IDs = [ele['ID'] for ele in hlist]

    try:
        IDs = [ele[0] for ele in IDs]
    except:
        pass
    mapp[D < 1] = 0
    index = list(zip(*np.where((D > 0) & (mapp > np.percentile(mapp[mapp > 0], 95)))))

    print("percentile", np.percentile(mapp, 95))
    print("percentileAdj", np.percentile(mapp[mapp > 0], 95))

    rData = [rankdata(-mapp[i, :]) for i in range(len(mapp))]
    for i, j in index:
        if rData[i][j] < 10 and rData[j][i] < 10:
            #            if mapp[i,j]>=np.percentile(mapp[i,:][mapp[i,:]>0],90) and mapp[i,j]>=np.percentile(mapp[j,:][mapp[j,:]>0],90):
            #            if rIdata[j]<20 and rJdata[i]<20:
            #            if np.sum(mapp[i,:]>0)>1 and np.sum(mapp[j,:]>0)>1:
            #                if mapp[i,j]>=np.percentile(mapp[i,:][mapp[i,:]>0],99) and mapp[i,j]>=np.percentile(mapp[j,:][mapp[j,:]>0],99) and rIdata[j]<20 and rJdata[i]<20:
            try:
                xi = hlist[i]['x'][0]
                xj = hlist[j]['x'][0]
                yi = hlist[i]['y'][0]
                yj = hlist[j]['y'][0]
            except:
                xi = hlist[i]['x']
                xj = hlist[j]['x']
                yi = hlist[i]['y']
                yj = hlist[j]['y']
            link.add_path([IDs[i], IDs[j]], weight=mapp[i, j], tau=mapt[i, j])
            link.node[IDs[i]]['pos'] = [xi, yi]
            link.node[IDs[j]]['pos'] = [xj, yj]

    if flag_scatter and len(D) == 0 and 0:
        plt.figure()
        plt.plot(distance, off, 'ob', rasterized=True)
        plt.plot(distance_in, off_in, 'or', rasterized=True)
        plt.plot(np.arange(0, np.max(distance), 100), np.arange(0, np.max(distance), 100) * 1. / 400, 'r', lw=3)

    return link


def extract_components_from_functional(link, G, fname='test', show_plot=False):
    all_sub = []
    all_contour = []
    all_nodes = []
    all_refined = []
    all_triangle = []
    RNDG = copy.copy(G)
    H = G.to_undirected()

    allpoints = []
    toPlotList = []
    partition = community.best_partition(link)
    list_nodes = []
    for com in set(partition.values()):
        list_nodes.append([nodes for nodes in list(partition.keys()) if partition[nodes] == com])

    # list_nodes=it.graph_partition(link)

    for iidS, subgL in enumerate(list_nodes):

        subg = link.subgraph(subgL)

        if len(subg) > 10 and subg.size() > 10:  # len(subg)>max(len(link)*0.05,4):
            toPlotList.append(subg)

            all_sub.append(subg)
            pos = nx.get_node_attributes(subg, 'pos')
            points = [(pos[key][0], pos[key][1]) for key in list(pos.keys())]
            if len(np.unique([ele[0] for ele in points])) == 1 or len(np.unique([ele[1] for ele in points])) == 1:
                points.append((points[-1][0] + .01, points[-1][1] + .01))
            try:
                elements, hull_index = cv.hull(0.1, points, plot=1)

                def get_index(point, pos):
                    for index, value in pos.items():
                        if value[0] == points[point][0]: return index

                allpos = nx.get_node_attributes(H, 'pos')

                allpoints = []
                Pdict = np.asarray([key for key in list(allpos.keys())])
                P2test = np.asarray([[allpos[key][0], allpos[key][1]] for key in list(allpos.keys())])
                contained = np.zeros(len(P2test))
                for triangle in elements:
                    Pth = matplotlib.path.Path(triangle)
                    contained += Pth.contains_points(P2test)

                all_triangle.append(elements)

                all_nodes.append(Pdict[contained >= 1].tolist())

                x = [points[k][0] for k in hull_index]
                x.append(x[0])
                y = [points[k][1] for k in hull_index]
                y.append(y[0])

                all_contour.append([x, y])
                hull_index = [get_index(xy, pos) for xy in hull_index]
                all_nodes[-1].extend(hull_index)
                alltmp = allpoints  # +hull_index
                allpoints = list(set(alltmp))
                tmp = nx.isolates(H.subgraph(allpoints))
                if len(tmp) > 0:
                    for ele in tmp:
                        allpoints.remove(ele)

                if H.size() > 0:

                    if len(all_nodes[-1]) > 5000:
                        refined_components = it.graph_partition(H.subgraph(all_nodes[-1]), 0, 500)

                        all_refined.append(refined_components)
                    else:
                        refined_components = [all_nodes[-1]]
                        all_refined.append(refined_components)

                else:
                    all_refined.append([all_nodes[-1]])
            except:
                pass

    if show_plot or 1:

        #        for iidS,subgL in enumerate(it.graph_partition(link)):
        for iidS, subgL in enumerate(toPlotList):
            subg = link.subgraph(subgL)
            functional_graph_plot(subg, link)
    flagI = 0
    all_refined_old = all_refined
    while flagI:
        flagI = 1
        all_refined_new = []

        for ele1 in all_refined_old:
            item1 = set([item for ele in ele1 for item in ele])
            flag2 = 0
            for ele2 in all_refined_old:
                item2 = set([item for ele in ele2 for item in ele])
                if len(set(item1).intersection(set(item2))) > min(len(item1), len(item2)) * .1:
                    all_refined_new.append([list(set(item1).union(set(item2)))])
                    flagI = 0
                    flag2 = 1
            if flag2 == 1:
                if ele1 not in all_refined_new: all_refined_new.append(ele1)
        all_refined_old = all_refined_new
        if flagI:
            all_refined = all_refined_new
            break
    th_motif = 6
    if fname == 'test':
        fname = '%04d' % (np.random.random() * 1000)
    fname = os.path.splitext(fname)[0]
    np.save(open('/home/alonardoni/8192correlation/G_components/GNL_RFCOM_' + fname + '%d.npy' % th_motif, 'w'),
            all_refined)
    pickle.dump(link,
                open('/home/alonardoni/8192correlation/G_components/GNL_link_' + fname + '%d.npy' % th_motif, 'w'))
    pickle.dump(G, open('/home/alonardoni/8192correlation/G_components/HNL_link_' + fname + '%d.npy' % th_motif, 'w'))
    test_CC = dict()
    test_CC_nodes = dict()
    if 0:
        if H.size() > 0:
            test_CC, test_CC_nodes, com2test = count_motifs(th_motif, H, all_refined)
        else:
            test_CC, test_CC_nodes, com2test = count_motifs(th_motif, link, all_refined)
        fCOM = com2test
        np.save(open('/home/alonardoni/8192correlation/G_motifs/GNL_STmot_' + fname + '%d.npy' % th_motif, 'w'),
                test_CC)
        np.save(open('/home/alonardoni/8192correlation/G_motifs/GNL_fCOM_' + fname + '%d.npy' % th_motif, 'w'),
                com2test)
        if H.size() > 0:
            test_CC, com2test = add_entire_graph(test_CC, H, all_nodes, th_motif, fCOM)
            test_CC, com2test = add_entire_graph_random(test_CC, H, all_nodes, th_motif, fCOM)
        else:
            test_CC, com2test = add_entire_graph(test_CC, link, all_nodes, th_motif, fCOM)
            test_CC, com2test = add_entire_graph_random(test_CC, link, all_nodes, th_motif, fCOM)
        np.save(open('/home/alonardoni/8192correlation/G_motifs/GNL_STmot_' + fname + '%d.npy' % th_motif, 'w'),
                test_CC)
        np.save(open('/home/alonardoni/8192correlation/G_motifs/GNL_sCOM_' + fname + '%d.npy' % th_motif, 'w'),
                com2test)
        motif_hist2(test_CC)

    print("---------------------------------------------------------------------")
    print("---------------------------------------------------------------------")
    print("---------------------------------------------------------------------")
    print("---------------------------------------------------------------------")
    return all_refined, allpoints, all_nodes, all_contour, all_sub, RNDG, all_triangle, test_CC, test_CC_nodes


#
def count_motifs(th_motif, H, all_refined):
    start = time.time()
    count = 0
    test_CC = {}
    test_CC_nodes = {}
    all_refined = [item for ele in all_refined for item in ele]

    all_refined = sorted(all_refined, key=lambda x: len(x))
    all_refined = all_refined[::-1]
    com2test = []
    pos = nx.get_node_attributes(H, 'pos')
    for ele in all_refined:
        allpoints = reduce_community(ele, pos)
        K = copy.deepcopy(H.subgraph(allpoints))

        allpoints = [ele for ele in K.nodes() if ele not in nx.isolates(K)]

        if H.subgraph(allpoints).size() > 10:
            com2test.append(H.subgraph(allpoints).copy())
            print("count", count, "# nodes --------------------->", len(allpoints))
            print("count", count, "# edges --------------------->", K.size())
            print("\t price to pay", len(K.edges()) * 1. / len(nx.complete_graph(max(len(allpoints), 2)).edges()))
            #            task2exec.append([mp2.motif_main,[th_motif,H.subgraph(allpoints).to_undirected(),count]])
            # test_CC[count]=test_CC[count]=m.motif_counts(6,H.subgraph(subg.nodes()))
            test_CC[count] = mp2.motif_main_IG(th_motif, H.subgraph(allpoints).to_undirected(), count)[
                0]  # test_CC[count]=m.motif_counts(6,H.subgraph(subg.nodes()))
            test_CC_nodes[len(allpoints)] = count
            count += 1
            if count > 5:
                break

    print("ELAPSED:", time.time() - start)
    return test_CC, test_CC_nodes, com2test


def add_entire_graph(test_CC, H, all_nodes, th_motif, fCOM):
    print("ENTIRE GRAPH")
    pool = H.nodes()
    for ele in all_nodes:
        for item in ele:
            try:
                pool.remove(item)
            except:
                pass

    refined = it.graph_partition(H.subgraph(pool), 0, 20)

    refined = sorted(refined, key=len)[::-1]

    count = -1

    count += 1
    com2test = []
    pos = nx.get_node_attributes(H, 'pos')
    for nMax, ele in enumerate(refined):
        if nMax < len(fCOM):
            ele = reduce_community(ele, pos, len(fCOM[nMax]))
            K = H.subgraph(ele)
            print("Gcount", count, "# nodes --------------------->", len(ele))
            print("Gcount", count, "# edges --------------------->", K.size())
            com2test.append(K.copy())
            allpoints = [ele for ele in K.nodes() if ele not in nx.isolates(K)]
            print("Gcount", count, "# Cnodes --------------------->", len(H.subgraph(allpoints)))
            print("\t price to pay", len(K.edges()) * 1. / len(nx.complete_graph(max(len(allpoints), 2)).edges()))

            test_CC['SC-%d' % count] = mp2.motif_main_IG(th_motif, H.subgraph(allpoints).to_undirected(), count)[
                0]  # test_CC[count]=m.motif_counts(6,H.subgraph(subg.nodes()))

            count += 1
            if count > 5:
                break
    return test_CC, com2test


def add_entire_graph_random(test_CC, H, all_nodes, th_motif, fCOM):
    print("ENTIRE GRAPH")
    flatRemoveList = [item for ele in all_nodes for item in ele]
    pool = [node for node in H.nodes() if node not in flatRemoveList]
    HH = H.subgraph(pool)

    pos = nx.get_node_attributes(HH, 'pos')
    refined = [[] for _ in range(9)]
    for key, value in pos.items():
        count = 0
        for i in range(3):
            for j in range(3):
                if i < value[0] * 3 < i + 1 and j < value[1] * 3 < j + 1:
                    refined[count].append(key)
                count += 1
    refined = sorted(refined, key=len)[::-1]
    count = -1

    count += 1
    com2test = []
    pos = nx.get_node_attributes(H, 'pos')
    for nMax, ele in enumerate(refined):
        if nMax < len(fCOM):
            ele = reduce_community(ele, pos, len(fCOM[nMax]))
            K = H.subgraph(ele)
            print("Gcount", count, "# nodes --------------------->", len(ele))
            print("Gcount", count, "# edges --------------------->", K.size())
            com2test.append(K.copy())
            allpoints = [ele for ele in K.nodes() if ele not in nx.isolates(K)]
            print("Gcount", count, "# Cnodes --------------------->", len(H.subgraph(allpoints)))
            p2pay = len(K.edges()) * 1. / len(nx.complete_graph(max(len(allpoints), 2)).edges())
            print("\t price to pay", p2pay)

            test_CC['RC-%d' % count] = mp2.motif_main_IG(th_motif, H.subgraph(allpoints).to_undirected(), K.size())[
                0]  # test_CC[count]=m.motif_counts(6,H.subgraph(subg.nodes()))
            RR = nx.erdos_renyi_graph(len(K), p2pay)
            test_CC['RR-%d' % count] = mp2.motif_main_IG(th_motif, RR, RR.size())[0]
            count += 1

    return test_CC, com2test


def plot_concave_components(H, all_contour, all_nodes, ISs=[]):
    plt.figure(300)
    allpos = nx.get_node_attributes(H, 'pos')
    performance = []
    colors = []
    CCC = len(all_contour) + 1
    for j in range(len(all_contour)):

        i = j
        ct2analyze = all_contour[i]

        plt.plot(ct2analyze[0], ct2analyze[1], 'k--', lw=2)
        rC = 0
        for elements in all_nodes[i]:
            if len(elements):
                NE = len(all_nodes[i])
                rC += 1
                x = np.array([allpos[ele][0] for ele in elements if ele is not None])
                y = np.array([allpos[ele][1] for ele in elements if ele is not None])
                IS2plot = []
                ISs_new = []
                for points in ISs:
                    flag = 0
                    if len(points[0]):
                        for kkk in range(len(x)):
                            xx = x[kkk]
                            yy = y[kkk]
                            if abs(xx - points[0][0]) ** 2 + abs(yy - points[1][0]) ** 2 < .05 ** 2:

                                if tuple(points) not in IS2plot:
                                    IS2plot.append(points)
                                    flag = 1
                    if not flag:
                        ISs_new.append(points)

                ISs = ISs_new
                if len(IS2plot):
                    CCC -= 1
                    colorFCOM = list(plt.get_cmap('jet')(CCC * 1. / len(all_contour)))
                    colors = []
                    for i in range(20):
                        colors.append(colorFCOM)
                        colors[-1][-1] = 1 - i * 1. / 20
                        colors[-1] = tuple(colors[-1])
                    for item in IS2plot:
                        a = item[0][:20]
                        b = item[1][:20]
                        plot_color_coded_line(item)
                        #                        plt.plot(a,b,zorder=-10,color='k')#plt.get_cmap('hot')(1.-rC*.99/NE+.2))
                        plt.scatter(a, b, zorder=-10, color=colors)
                colors.append(plt.get_cmap('hot')(1. - rC * .99 / NE + .2))
                performance.append(len(IS2plot))

    for points in ISs:
        a = points[0][:20]
        b = points[1][:20]
        plt.plot(a, b, zorder=-10, color='k')
        plt.scatter(a, b, zorder=-10, color=plt.get_cmap('bone')(np.arange(20) * 1. / 20))
    colors.append('k')
    performance.append(len(ISs))
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.figure()
    Norm = max(sum(performance), 1)
    performanceN = [item * 1. / Norm for item in performance]
    plt.bar(np.arange(len(performanceN)), performanceN, color=colors)
    plt.xlabel('FCs')
    plt.suptitle('IS detection performance')
    return performanceN


def plot_concave_components_LW(H, all_contour, all_nodes, ISs=[], composition=[]):
    plt.figure('LW_concave_components')
    allpos = nx.get_node_attributes(H, 'pos')
    performance = []
    colors = []
    CCC = len(all_contour) + 1
    clusterAssigned = np.ones(len(ISs)) * -1

    for j in range(len(all_contour)):

        i = j
        ct2analyze = all_contour[i]
        plt.plot(ct2analyze[0], ct2analyze[1], 'k', lw=2)
        rC = 0
        posComposition = []
        points=[[],[]]
        for cluster in np.unique(composition):
            clustered = np.where(composition == cluster)[0]
            helperXY = [(points[0][0], points[1][0]) for points in np.asarray(ISs)[clustered]]
            posComposition.append(np.mean(helperXY, axis=0))
        for elements in all_nodes[i]:
            if len(elements):

                rC += 1
                x = np.array([allpos[ele][0] for ele in elements if ele is not None])
                y = np.array([allpos[ele][1] for ele in elements if ele is not None])


                for pivot, cluster in enumerate(np.unique(composition)):
                    trjN = 0
                    if len(points[0]):
                        tmpX = x - posComposition[pivot][0]
                        tmpY = y - posComposition[pivot][1]
                        if np.min(tmpX ** 2 + tmpY ** 2) < .05 ** 2:

                            to_plot = np.where(composition == cluster)[0]
                            clusterAssigned[to_plot] += 1
                            CCC -= 1

                            for offset, item in enumerate(np.asarray(ISs)[to_plot]):
                                colorFCOM = list(plt.get_cmap('jet')(
                                    (CCC + .1 * clusterAssigned[to_plot][offset]) * 1. / len(all_contour)))
                                a = item[0][:50]
                                b = item[1][:50]

                                #                            plt.plot(a,b,color=colorFCOM)#plt.get_cmap('hot')(1.-rC*.99/NE+.2))
                                plot_color_coded_line(item)
                                #                            plt.colorbar()
                                plt.plot(a[0], b[0], 'o', color=colorFCOM)
                            trjN += len(to_plot)
                    colors.append(list(plt.get_cmap('jet')(CCC * 1. / len(all_contour))))
                    performance.append(trjN)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.figure('NNcluster_components')
    pivot = np.where(clusterAssigned == -1)[0]
    for points in np.asarray(ISs)[pivot]:
        a = points[0][:50]
        b = points[1][:50]
        #                plt.plot(a,b,'k-')
        #                plt.plot(a,b,'ok--',alpha=.1)
        plot_color_coded_line(points)
        #                plt.plot(a,b,color='k',zorder=-1)#plt.get_cmap('hot')(1.-rC*.99/NE+.2))
        plt.plot(a[0], b[0], 'o', color='k', zorder=-1)
    # plt.scatter(a,b,s=np.arange(len(a))[::-1]*100,c='k')
    #            except:
    #                print "!:",points
    #                plt.scatter(points[0],points[1],alpha=.5,s=10,zorder=-100,c='k')
    colors.append('k')
    performance.append(len(pivot))
    plt.xlim([0, 1])
    plt.ylim([0, 1])


def plot_color_coded_line(item):
    x = item[0][:50]
    y = item[1][:50]
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
                        norm=plt.Normalize(0, len(item[0])))

    lc.set_array(np.arange(len(x)))
    lc.set_linewidth(1)
    plt.gca().add_collection(lc)


def plot_concave_components_NC(H, all_contour, all_nodes, ISs=[]):
    plt.figure(9300)
    allpos = nx.get_node_attributes(H, 'pos')
    performance = []
    colors = []
    CCC = len(all_contour) + 1
    for j in range(len(all_contour)):

        i = j
        ct2analyze = all_contour[i]

        plt.plot(ct2analyze[0], ct2analyze[1], 'k--', lw=2)
        rC = 0
        for elements in all_nodes[i]:
            if len(elements):
                NE = len(all_nodes[i])
                rC += 1
                x = np.array([allpos[ele][0] for ele in elements if ele is not None])
                y = np.array([allpos[ele][1] for ele in elements if ele is not None])
                IS2plot = []
                ISs_new = []
                for points in ISs:
                    flag = 0
                    if len(points[0]):
                        for kkk in range(len(x)):
                            xx = x[kkk]
                            yy = y[kkk]
                            if abs(xx - points[0][0]) ** 2 + abs(yy - points[1][0]) ** 2 < .05 ** 2:

                                if tuple(points) not in IS2plot:
                                    IS2plot.append(points)
                                    flag = 1
                    if not flag:
                        ISs_new.append(points)

                ISs = ISs_new
                if len(IS2plot):
                    CCC -= 1
                    colorFCOM = list(plt.get_cmap('jet')(CCC * 1. / len(all_contour)))
                    colors = []
                    for i in range(20):
                        colors.append(colorFCOM)
                        colors[-1][-1] = 1 - i * 1. / 20
                        colors[-1] = tuple(colors[-1])
                    for item in IS2plot:
                        a = item[0][:20]
                        b = item[1][:20]
                        plot_color_coded_line(item)
                        #                        plt.plot(a,b,zorder=-10,color=colors[0])#plt.get_cmap('hot')(1.-rC*.99/NE+.2))
                        plt.plot(a[0], b[0], 'o', zorder=-10, color=colors[0])
                colors.append(plt.get_cmap('hot')(1. - rC * .99 / NE + .2))
                performance.append(len(IS2plot))


def functional_graph_plot(subg, link, flag=0):
    w = []
    tau = []
    limits = [[], []]
    pos = list(nx.get_node_attributes(link, 'pos').values())
    x, y = list(zip(*pos))
    for ele in link.edges(data=True):
        w.append(ele[2]['weight'])
        tau.append(abs(ele[2]['tau']))

    limits[1].extend([np.min(w), np.max(w)])
    limits[0].extend([np.min(tau), np.max(tau)])
    w = []
    tau = []
    for ele in subg.edges(data=True):
        w.append(ele[2]['weight'])
        tau.append(abs(ele[2]['tau']))

    if 73 in plt.get_fignums():
        flag = 1
    # fig=plt.figure(73)
    # it was
    if len(subg.nodes()) < 10:  # len(subg.nodes())<len(link.nodes())*0.05

        fig = plt.figure(120)
        ax = fig.add_subplot(211)
        im = nx.draw_networkx_edges(subg, pos=nx.get_node_attributes(subg, 'pos'), edge_color=w, edge_vmin=limits[1][0],
                                    edge_vmax=limits[1][1], width=2, ax=ax)
        ax.set_title('sparse CC, weight')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        if flag == 1:
            fig.subplots_adjust(right=0.8)
            pos1 = ax.get_position()
            cbar_ax = fig.add_axes([0.85, (pos1.y0 + pos1.height * .5 - 0.35 * .5), 0.02, 0.35])
            fig.colorbar(im, cax=cbar_ax)

        ax = fig.add_subplot(212)
        im = nx.draw_networkx_edges(subg, pos=nx.get_node_attributes(subg, 'pos'), edge_color=tau,
                                    edge_vmin=limits[0][0], edge_vmax=limits[0][1], width=2, ax=ax)

        plt.plot(x, y, 'ob', zorder=-1000)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_title('sparse CC, tau')
        if flag == 1:
            fig.subplots_adjust(right=0.8)
            pos1 = ax.get_position()
            cbar_ax = fig.add_axes([0.85, (pos1.y0 + pos1.height * .5 - 0.35 * .5), 0.02, 0.35])
            fig.colorbar(im, cax=cbar_ax)

    else:
        if 121 in plt.get_fignums():
            flag = 1
        limits = [[], []]
        w = []
        tau = []
        for ele in subg.edges(data=True):
            w.append(ele[2]['weight'])
            tau.append(abs(ele[2]['tau']))

        limits[1].extend([np.min(w), np.max(w)])
        limits[0].extend([np.min(tau), np.max(tau)])
        fig = plt.figure(121)
        ax = fig.add_subplot(111)
        ax.set_title('dense CC, weight')
        plt.plot(x, y, 'ob', zorder=-1000)
        im = nx.draw_networkx_edges(subg, pos=nx.get_node_attributes(subg, 'pos'), edge_color=w, edge_vmin=limits[1][0],
                                    edge_vmax=limits[1][1], width=2, ax=ax)
        fig.subplots_adjust(right=0.8)
        pos1 = ax.get_position()
        cbar_ax = fig.add_axes([0.85, (pos1.y0 + pos1.height * .5 - 0.35 * .5), 0.02, 0.35])
        if not flag: fig.colorbar(im, cax=cbar_ax)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        fig = plt.figure(122)
        ax = fig.add_subplot(111)
        ax.set_title('dense CC, tau')
        im = nx.draw_networkx_edges(subg, pos=nx.get_node_attributes(subg, 'pos'), edge_color=tau,
                                    edge_vmin=limits[0][0], edge_vmax=limits[0][1], width=2, ax=ax)
        fig.subplots_adjust(right=0.8)
        pos1 = ax.get_position()
        plt.plot(x, y, 'ob', zorder=-1000)
        cbar_ax = fig.add_axes([0.85, (pos1.y0 + pos1.height * .5 - 0.35 * .5), 0.02, 0.35])

        if not flag: fig.colorbar(im, cax=cbar_ax)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
