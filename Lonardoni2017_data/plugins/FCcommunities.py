import matplotlib as mpl

mpl.use('Agg')
import tools.tools as t
import networkx as nx
from difflib import get_close_matches
from multiPP_ST_ref import runFunctionsInParallel
from matplotlib._pylab_helpers import Gcf
import tools.cor2is as c2i
import tools.parallelAnalysis as PA
import tools.utilSpkTrain as UST
import tools.igraph_tools as it
import pylab as plt

import os

import numpy as np


UST.mfrMin = .1
mpl.use('Agg')
plt.ioff()
import community

def extract_SG(G, mapLag, mapPeak, index):
    """
    E.G. index=np.where((abs(posx)>-.5) & (abs(posy)>-.5))[0]
    """
    g = G.subgraph(index)
    mapping = {}
    j = 0
    for i in g.nodes():
        mapping[i] = j
        j += 1

    g = nx.relabel_nodes(g, mapping)
    mapLag = mapLag[np.ix_(index, index)]
    mapPeak = mapPeak[np.ix_(index, index)]
    return g, mapPeak, mapLag


def find_CC_IS(folder, fName):

    fMapLag = 'Results/' + folder + '/' + fName[:-4] + '_mapt2.npy'
    fMapPeak = 'Results/' + folder + '/' + fName[:-4] + '_mapp2.npy'
    fIS = 'Results/' + folder + '/' + fName[:-4] + '_IS2.npy'
    if fName.rfind('mat') > 0:

        fMapLag = folder + '/' + fName[:-4] + '_mapt2.npy'
        fMapPeak = folder + '/' + fName[:-4] + '_mapp2.npy'
        fIS = folder + '/' + fName[:-4] + '_IS2.npy'
    return fMapLag, fMapPeak, fIS


def find_G(folder, fName):
    path2graph = 'Results/' + folder + '/Graph/'
    return path2graph + get_close_matches(fName, os.listdir(path2graph))[0]


def construct_CC_maps(folder, n_work=10):
    if folder.rfind('res') > 0 or folder.rfind('1024') > 0:
        folder += '/Spike/'
        for fName in os.listdir(folder):
            if fName.rfind('.npy') > 0:
                PA.construct_maps(folder, fName, n_work)
    else:
        for fName in os.listdir(folder):
            if fName.rfind('.mat') > 0:
                if fName.rfind('asal') > 0 or fName.rfind('00') > 0 or fName.rfind('pont') > 0:
                    print folder + fName
                    PA.construct_maps(folder, fName, n_work)
                    print "done"


def construct_ISs_maps(folders, n_work=10, destPath=''):
    fNameList = []
    for folder in folders:
        if folder.rfind('res') > 0 or folder.rfind('1024') > 0:
            folder += '/Spike/'
            for fName in os.listdir(folder):
                if fName.rfind('.npy') > 0 > fName.rfind('map') and fName.rfind('STIM') < 0 and fName.rfind(
                        'IS') < 0:
                    fNameList.append(folder + fName)

        else:
            for fName in os.listdir(folder):
                if fName.rfind('.mat') > 0:
                    if fName.rfind('asal') > 0 or fName.rfind('00') > 0 or fName.rfind('pont') > 0:
                        fNameList.append(folder + fName)
    PA.construct_ISs(fNameList, n_work, destPath)
    print "done"


def init(folder, fName):
    if fName.rfind('npy') > 0:
        node2remove = []
        print folder
        print fName
        print find_G(folder, fName)
        G, pos = t.OpenGraph(folder=find_G(folder, fName), flag=1)

        for i in G.nodes():
            if pos[i][0] > 1.55 or pos[i][1] > 1.55:  # or pos[i][0]<.2 or pos[i][1]<.2:
                G.remove_edges_from(G.edges(i))
                G.remove_edges_from(G.in_edges(i))
                #                G.remove_node(i)
                node2remove.append(i)

        hlist = np.load('Results/' + folder + '/Spike/' + fName).tolist()[1:]
    # for i in sorted(node2remove)[::-1]:
    #            hlist[i]['x']=[10**10*np.random.random()]
    else:
        G = nx.Graph()
        pos = dict()
        hlist = UST.LoadSpikeTrains(folder + '/' + fName, ReadHeader=True)

        for ele in hlist:
            ele['ID'] = 64 * (ele['x'] - 1) + ele['y']
            G.add_node(ele['ID'])
            ele['x'] *= 1. / 64
            ele['y'] *= 1. / 64
            pos[ele['ID']] = [ele['x'], ele['y']]
            G.node[ele['ID']]['pos'] = pos[ele['ID']]
            ele['x'] = [ele['x']]
            ele['y'] = [ele['y']]

    fMapt, fMapp, fIS = find_CC_IS(folder, fName)

    mapt = np.load(fMapt)
    mapp = np.load(fMapp)
    posx, posy = np.asarray(zip(*pos.values()))

    for i in xrange(len(mapp)):
        #        print mapp[i,i],
        mapp[i, i] = 0
    mapp[mapp == 1] = 0
    ISs = np.load(fIS).tolist()
    composition = np.load(fIS[:-4] + '_cut.npy')

    H = G.to_undirected()

    lh = int(len(hlist))
    D = np.zeros((lh, lh))
    for i in xrange(lh):
        for j in xrange(lh):
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
            if abs(xi - xj) ** 2 + abs(yi - yj) ** 2 < .25 ** 2:
                D[i, j] = 1


    return hlist, lh, G, H, posx, posy, mapp, mapt, D, ISs, composition


def findMaxFCsPP(folder, fName='', txL=[], n_work=10, destPath='FIG_final/'):
    if fName == '':
        fName = \
            [ele for ele in os.listdir('Results/' + folder + '/Spike/') if fName.rfind('.npy') > 0][
                0]

    hlist, lh, G, H, posx, posy, mapp, mapt, D, _, _ = init(folder, fName)

    def test_tx(tx):

        link = c2i.extract_functional_graph3(G, hlist, mapp, mapt, tx, D)

        NC = sum([1 for subg in it.graph_partition(link) if 100 > len(subg) > 20])
        GC = sum([-1 for subg in it.graph_partition(link) if len(subg) > 100]) * .5

        for subg in nx.connected_component_subgraphs(link):
            c2i.functional_graph_plot(subg, link, flag=0)
        plt.figure(121)
        plt.savefig('testFC/' + 'dense_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
        plt.close(121)
        plt.figure(120)
        plt.savefig('testFC/' + 'sparse_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
        plt.close(120)
        for subg in (it.graph_partition(link)):
            subg = link.subgraph(subg)
            c2i.functional_graph_plot(subg, link, flag=0)
        plt.figure(121)
        plt.savefig('testFC/' + 'denseIT_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
        plt.close(121)
        plt.figure(120)
        plt.savefig('testFC/' + 'sparseIT_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
        plt.close(120)
        return tx, NC + GC + 20

    task2exec = [[test_tx, [tx], ] for tx in txL]
    out = runFunctionsInParallel(listOf_FuncAndArgLists=task2exec, maxAtOnce=n_work)
    plt.close('all')
    plt.figure()
    txs = [ele[-2] for ele in out[1]]
    Nfcs = [ele[-1] for ele in out[1]]
    plt.scatter(txs, Nfcs)
    try:
        new_tx = np.load(destPath + '/tx_all5.npy').item()
    except:
        new_tx = dict()
    new_tx[folder + '/new' + fName] = (txs, Nfcs)
    np.save(open(destPath + '/tx_all5.npy', 'w'), new_tx)
    plt.xlabel('thresholds')
    plt.ylabel('# of FCs > 5 nodes')
    res = sorted(out[1], key=lambda x: x[-1])

    return [fName, res[-1][0], res[-1][1]]


def findMaxFCs(folder, fName='', txL=[]):
    if fName == '':
        fName = \
            [ele for ele in os.listdir('Results/' + folder + '/Spike/') if fName.rfind('.npy') > 0][
                0]

    NC = 0
    ftx = 0
    hlist, lh, G, H, posx, posy, mapp, mapt, D, _, _ = init(folder, fName)

    for tx in txL:

        link = c2i.extract_functional_graph3(G, hlist, mapp, mapt, tx, D)
        if sum([1 for subg in nx.connected_component_subgraphs(link) if 40 > len(subg.nodes()) > 20]) > NC:
            ftx = tx
            NC = sum([1 for subg in nx.connected_component_subgraphs(link) if 40 > len(subg.nodes()) > 20])
    return [fName, ftx, NC]



def finalize(folder, fName, tx, destPath='FIG_final/'):
    hlist, lh, G, H, posx, posy, mapp, mapt, D, ISs, composition = init(folder, fName)
    g = G.copy()
    link = c2i.extract_functional_graph3(g, hlist, mapp, mapt, tx, D, flag_scatter=True)
    Ledges = link.size()
    Nnodes = len([node for node in link.nodes() if any(link.edges(node))])
    print "------------------->", (Ledges * 1. / Nnodes) * 2. / (25 ** 2 * .01)
    print "------------------->", (Ledges * 1. / Nnodes) * 2


    all_refined, allpoints, all_nodes, all_contour, all_sub, RNDG, all_triangle, test_CC, test_CC_nodes = c2i.extract_components_from_functional(
        link, G, folder.replace('/', '_') + '_new' + fName, show_plot=True)


    perf = c2i.plot_concave_components(H, all_contour, all_refined, ISs=ISs)
    c2i.plot_concave_components_NC(H, all_contour, all_refined, ISs=ISs)
    c2i.plot_concave_components_LW(H, all_contour, all_refined, ISs=ISs, composition=composition)
    try:
        new_perf = np.load(destPath + '/perf_all5.npy').item()
    except:
        new_perf = dict()
    new_perf[folder + '/new' + fName] = perf
    np.save(open(destPath + '/perf_all5.npy', 'w'), new_perf)

    figures = [manager.canvas.figure for manager in Gcf.get_all_fig_managers()]
    plt.ioff()
    for i, figure in enumerate(figures):
        figure.savefig(destPath + '/' + fName[:-4] + '_6M_2_%02d_5.svg' % i)
    plt.close('all')
    try:
        all_motifs = np.load(destPath + '/all_motifs5.npy').item()
    except:
        all_motifs = dict()
    all_motifs[folder + '/' + fName[:6]] = [test_CC, test_CC_nodes]
    np.save(open(destPath + '/all_motifs5.npy', 'w'), all_motifs)
    for subg in nx.connected_component_subgraphs(link):
        c2i.functional_graph_plot(subg, link, flag=0)
    plt.figure(121)
    plt.savefig('testFC/' + 'dense_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(121)
    plt.figure(120)
    plt.savefig('testFC/' + 'sparse_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(120)
    for subg in (it.graph_partition(link)):
        subg = link.subgraph(subg)
        c2i.functional_graph_plot(subg, link, flag=0)
    plt.figure(121)
    plt.savefig('testFC/' + 'denseIT_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(121)
    plt.figure(120)
    plt.savefig('testFC/' + 'sparseIT_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(120)
    partition = community.best_partition(link)
    list_nodes = []
    for com in set(partition.values()):
        list_nodes.append([nodes for nodes in partition.keys() if partition[nodes] == com])
    for subg in list_nodes:
        subg = link.subgraph(subg)
        c2i.functional_graph_plot(subg, link, flag=0)
    plt.figure(121)
    plt.savefig('testFC/' + 'denseCOM_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(121)
    plt.figure(120)
    plt.savefig('testFC/' + 'sparseCOM_' + fName[:-4] + '_6M_2_%05d_5.png' % (tx * 100))
    plt.close(120)
    return all_motifs


def build_service(n_work=10, IS=0, CC=0, folders=['Results/testNMDAfinal'],
                  destPath='/home/alonardoni/8192correlation/TESTfigALL/'):
    if IS: construct_ISs_maps(folders, n_work=n_work, destPath=destPath)
    if CC:
        for folder in folders:
            construct_CC_maps(folder, n_work=n_work)


def run_folder(folder='test4FC_sG_dR', destPath='FIG_final/'):
    #    construct_CC_maps(folder,n_work=n_work)
    #    construct_ISs_maps(folder,n_work=n_work)
    plt.ioff()
    if folder.rfind('/') > 0:
        path2search = folder
    else:
        path2search = 'Results/' + folder + '/Spike/'
    all_out = {}
    for fName in os.listdir(path2search):

        if fName.rfind('npy') > 0:
            if fName.rfind('map') < 0 and fName.rfind('STIM') < 0 and fName.rfind('IS') < 0:
                plt.close('all')
                tx = -10
                NC = 10
                print "----> Simulation", fName
                print "----> Folder", folder
                print "TX:", tx
                print "NC:", NC
                print "cutoff: %.02f \t FCs: %02d" % (tx, NC)
                if NC > 0:
                    out = finalize(folder, fName, tx, destPath=destPath)
                    all_out[fName] = out
        elif fName.rfind('mat') > 0:
            if fName.rfind('map') < 0 and fName.rfind('STIM') < 0 and fName.rfind('IS') < 0:
                if fName.rfind('asal') > 0 or fName.rfind('00') > 0 or fName.rfind('pont') > 0:
                    plt.close('all')
                    tx = -10
                    NC = 10
                    print "----> Simulation", fName
                    print "----> Folder", folder
                    print "TX:", tx
                    print "NC:", NC
                    print "cutoff: %.02f \t FCs: %02d" % (tx, NC)
                    if NC > 0:
                        out = finalize(folder, fName, tx, destPath=destPath)
                        all_out[fName] = out

    return all_out


def build_run_save_Model_serv(folder):
    destPath = '/home/alonardoni/8192correlation/result_model/'

    dirName = os.path.dirname('Results/' + folder + '/Spike/').split('alonardoni/')[-1]
    dest = destPath + '/' + dirName.replace('/', '_') + '/'

    if not os.path.isdir(dest):
        build_service(n_work=40, IS=0, CC=1, folders=['Results/' + folder + '/'], destPath=destPath)
    else:
        build_service(n_work=40, IS=0, CC=1, folders=['Results/' + folder + '/'], destPath=destPath)
    out = run_folder(folder=folder, destPath=dest)
    np.save(open('/home/alonardoni/8192correlation/result_model/' + folder + '_positive_3_1.npy', 'w'), out)


def build_run_save_Model(folder):
    destPath = '/home/alonardoni/8192correlation/result_model/'

    build_service(n_work=40, IS=0, CC=0, folders=['Results/' + folder + '/'], destPath=destPath)
    dirName = os.path.dirname('Results/' + folder + '/Spike/').split('alonardoni/')[-1]
    dest = destPath + '/' + dirName.replace('/', '_') + '/'

    out = run_folder(folder=folder, destPath=dest)
    np.save(open('/home/alonardoni/8192correlation/result_model/' + folder + '_positive_3_1.npy', 'w'), out)


def build_run_save_Exp(folder):
    destPath = '/home/alonardoni/8192correlation/result_model/'

    dirName = ('/home/alonardoni/8192correlation/EXPDATA/' + folder).split('alonardoni/')[-1]
    dest = destPath + '/' + dirName.replace('/', '_') + '/'
    run_folder(folder='/home/alonardoni/8192correlation/EXPDATA/' + folder + '/', destPath=dest)


if __name__ == "__main__":

    for ele in os.listdir('EXPDATA/'):
        build_run_save_Exp(ele)

    for ele in ['test4FC_dG_dR', 'test4FC_sG_dR', 'test4FC_dG_dR_2']:
        build_run_save_Model(ele)
