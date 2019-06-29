__author__ = 'dlonardoni'
import networkx as nx
import os
import pylab as plt
import numpy as np
import pickle
import utilSpkTrain as UST
from matplotlib.collections import LineCollection

mfrMax = 15  # (Hz) - maximum acceptable MFR
mfrMin = 0.1  # (Hz) - minumum acceptable MFR


# @profile
def SaveSpikeTrains(CellList, HEADER, fname='SpikeTrain.txt', folder='/', flag=False, random=False):
    """
        NOTE: introduce HDF5 files to save data
        TO DO
        TO DO
        TO DO
    """
    spikelist = list()
    spikelist.append(HEADER)
    for cell in CellList:

        formatName = {'type': cell.ID_syn, "spikes": np.array(cell.record["spk"]), "ID": np.array([cell.num]),
                      "x": np.array([cell.pos["x"]]), "y": np.array([cell.pos["y"]])}
        if random:
            formatName = {'type': cell.ID_syn, "spikes": np.array(cell.record["rspk"]), "ID": np.array([cell.num]),
                          "x": np.array([cell.pos["x"]]), "y": np.array([cell.pos["y"]])}
        if flag:
            param = {'tau_w': cell.AdExp.tau_w, 'G_l': cell.AdExp.G_l, 'a': cell.AdExp.a, 'b': cell.AdExp.b,
                     'C': cell.AdExp.C, 'E_l': cell.AdExp.E_l, 'V_thre': cell.AdExp.V_thre,
                     'Delta_T': cell.AdExp.Delta_T, 'V_reset': cell.AdExp.V_reset,
                     'I_ext': cell.AdExp.I_ext}  # 'type':cell.ID_syn,}
            formatName['parameters'] = param

        spikelist.append(formatName)
    path = os.getcwd()

    np.save(open(path + '/' + folder + fname + '.npy', 'w'), spikelist)


def ReadNpyTrains(fname='SpikeTrain.mat', folder='/', flag=1):
    """
        NOTE: Read .mat recordings in the format 2xN ordered by cellID
        NOTE: introduce HDF5 files to save data
        TO DO
        TO DO
        TO DO
    """
    path = os.getcwd()
    path = os.sep.join(path.split(os.sep)[:-1])
    if flag == 0:
        tmp = np.load(folder + fname + '.npy')
    else:
        tmp = np.load(path + folder + fname + '.npy')

    return tmp


# @profile
def convertposxN(CellList,ax='y'):
    tmpid = []
    tmptime = []
    for cell in CellList:
        tmpid.extend(np.ones(len(cell['spikes'])) * cell[ax])
        tmptime.extend(cell['spikes'])
    tmp = np.array([np.transpose(np.array(tmpid)), np.transpose(np.array(tmptime))])

    return tmp


def convert2xN(CellList):
    tmpid = []
    tmptime = []
    for cell in CellList:
        tmpid.extend(np.ones(len(cell['spikes'])) * cell['ID'])
        tmptime.extend(cell['spikes'])
    tmp = np.array([np.transpose(np.array(tmpid)), np.transpose(np.array(tmptime))])

    return tmp

def SaveGraph(G, fname='graph', folder='/'):
    path = os.getcwd()
    pickle.dump(G, open(path +'/' + folder + fname + '.graph', 'w'))
    # hysto_degree(G,fname,folder)


def OpenGraph(fname='graph', folder='/', flag=False):
    path = os.getcwd()
    path = os.sep.join(path.split(os.sep)[:-1])
    if flag:
        G = pickle.load(open(folder, 'r'))
    else:
        G = pickle.load(open(path +'/' + folder + fname + '.graph', 'r'))
    pos = nx.get_node_attributes(G, 'pos')
    return G, pos


def hysto_degree(G, fname='graph', folder='/'):
    """
        Procedure to plot the distribution of
            - distribution of ingoing & outgoing edges of each node
            - distribution of ingoing edges of each node
        Input:
        G:      graph in networkx format

        Output:
        Plot of the two distributions
    """
    path = os.getcwd()
    degree_sequence = sorted(nx.degree(G).values(), reverse=True)
    plt.figure()
    plt.plot(degree_sequence, 'b-', marker='o')
    plt.ylabel("Degree")
    plt.xlabel("Neuron #")

    degree_sequence = sorted(nx.DiGraph.in_degree(G).values(), reverse=True)
    plt.plot(degree_sequence, 'r-', marker='o')
    plt.savefig(path + folder + fname + 'degree.png', bbox_inches='tight', pad_inches=0)
    plt.close('all')




def detect(fname='nan', folder='/', tbin=5):
    import copy
    path = os.getcwd()
    path = os.sep.join(path.split(os.sep)[:-1])

    CellList = ReadNpyTrains(fname=fname, folder=folder + '/Spike/')
    ts = CellList[0]['ts'] * 1000
    tin = 0
    twind = int(25.0 / tbin)
    brtmp = []
    end = 0
    spikes = convert2xN(CellList[1:])
    rate = np.histogram(spikes[1], bins=np.arange(tin, ts, tbin))
    rate = rate[0]
    N = int(spikes[0][-1])
    burst = []
    ID = 0
    while (end + twind + tin) * tbin < ts:
        start = np.argmax(rate[end:] > 5. / 8192 * N) + end
        end = np.argmax(rate[start:] < 5. / 8192 * N) + start
        if end - start > twind and max(rate[start:end]) > N * 0.01:
            ID += 1
            for h in xrange(1, N):
                pivotstart = np.searchsorted(CellList[h]['spikes'], start * tbin, 'left')
                pivotend = pivotstart + np.searchsorted(CellList[h]['spikes'][pivotstart:], end * tbin, 'left')
                times = copy.deepcopy(CellList[h])
                times['spikes'] = CellList[h]['spikes'][pivotstart:pivotend]
                brtmp.append(times)

            burst.append(brtmp)
            brtmp = []
            # print ID,start*tbin,"-->",end*tbin,"max",(end-start)*tbin
        end += 1
    np.save(open(path + folder + '/Burst/' + fname + 'burst.npy', 'w'), burst)


# @profile
def trj_plot(burst, tbin=2, mswind=10):
    posx = np.zeros(len(burst))
    posy = np.zeros(len(burst))
    for cell in burst:
        posx[cell['ID']] = cell['x']
        posy[cell['ID']] = cell['y']
    burst = convert2xN(burst)
    twind = int(mswind / tbin)
    tstart = np.floor(min(burst[1]))
    tend = np.floor(max(burst[1] + 1))
    rate, t = np.histogram(burst[1], bins=np.arange(tstart, tend, tbin))
    rate[rate == 0] = 1
    meanx = np.histogram(burst[1], bins=np.arange(tstart, tend, tbin), weights=[posx[i] for i in burst[0]])[0]
    meany = np.histogram(burst[1], bins=np.arange(tstart, tend, tbin), weights=[posy[i] for i in burst[0]])[0]
    meanx = np.divide(meanx, rate)
    meany = np.divide(meany, rate)
    # tr=np.array([[np.mean(meanx[d-twind:d]),np.mean(meany[d-twind:d])] for d in xrange(twind,len(meanx)-1)])
    tmptr = np.array(
        [[np.mean(meanx[d - twind:d]), np.mean(meany[d - twind:d])] for d in xrange(twind, len(meanx) - 1)])
    tmptr = tmptr.reshape(-1, 1)
    # print time()-start
    segs = np.concatenate([tmptr[:-1], tmptr[1:]], axis=1)
    lc = LineCollection(segs, cmap=plt.get_cmap('jet'), lw=2)
    lc.set_array(t)
    plt.gca().add_collection(lc)  # add the collection to the plot


def trj_all(fname, gname, folder, mswind=10):
    G, pos = OpenGraph(fname=gname, folder=folder + '/Graph/')
    path = os.getcwd()
    path = os.sep.join(path.split(os.sep)[:-1])
    burst = np.load(path + folder + '/Burst/' + fname + 'burst.npy')
    plt.figure()
    plt.hold(True)
    nx.draw_networkx_edges(G, pos, alpha=0.02)
    for i in xrange(len(burst)):
        burstmp = convert2xN(burst[i])
        # print "Id burst:", i, "\t started at :", int(tin/1000), 's ',end,len(burst)
        # print "tin:", tstart, "tend:",tend
        trj_plot(burstmp, mswind=mswind)
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.savefig(path + folder + '/fig/' + fname + '_all.png')
    plt.close('all')




def SpikeStatWindow(fname, folder, tin, ts=0):
    """

    """
    strict = 1
    try:
        CellList = ReadNpyTrains(fname, folder)
    except:
        CellList = UST.LoadSpikeTrains(folder + fname, ReadHeader=True)
    if strict:
        for ele in CellList[1:]:
            if len(ele['spikes']) < mfrMin * CellList[0]['ts']:
                ele['spikes'] = []
    spikes = convert2xN(CellList[1:])

    if ts == 0:
        ts = CellList[0]['ts']
    del CellList
    spikes = np.array(
        [spikes[0][(spikes[1] > tin) & (spikes[1] < ts)], spikes[1][(spikes[1] > tin) & (spikes[1] < ts)]])
    # print "tin:", int(min(spikes[1]/1000)), "to:", int(max(spikes[1]/1000)),ts,tin
    print "CellList loaded and converted"
    out = UST.SpikeStatistics(spikes, TimeWindow=(ts - tin) / 1000)
    A = UST.ShowSpkStat(out)
    return A


def hist_mbl(path2file):
    """
        Procedure to compute a burst to burst correlation:
        - each synchronized event is isolated with its own duration

    """
    burst = np.load(path2file)
    if len(burst) == 0:
        print "empty burst"
        return
    nburst = len(burst)
    N = len(burst[0])
    tinit = np.ones(nburst) * 1e10
    tend = np.zeros(nburst)
    for j in xrange(nburst):
        for i in xrange(N):
            if len(burst[j][i]['spikes']) > 0:
                tinit[j] = int(min(tinit[j], burst[j][i]['spikes'][0]))
                tend[j] = int(max(tend[j], burst[j][i]['spikes'][-1]))
    H = np.histogram(tend[tend < 300000] - tinit[tend < 300000], bins=np.arange(0, 1000, 10))
    plt.bar(H[1][:-1] - 2, H[0], width=4, color='r')
    H = np.histogram(tend[tinit > 300000] - tinit[tinit > 300000], bins=np.arange(0, 1000, 10))
    plt.bar(H[1][:-1] + 2, H[0], width=4, color='b')

    plt.show()
