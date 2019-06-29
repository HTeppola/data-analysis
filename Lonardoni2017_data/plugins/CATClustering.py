import sys

import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['svg.fonttype'] = 'none'
from matplotlib._pylab_helpers import Gcf
import networkx as nx
import pylab as plt
import numpy as np
import scipy.cluster.hierarchy as sch
from tools import utilSpkTrain as UST
from scipy.interpolate import interp1d
from matplotlib import cm
from matplotlib.widgets import Slider, Button
from scipy.spatial.distance import squareform
from tools.parallelAnalysis import fastburstcorrelation
import SBEdetection as SBED
import os
import time
import tools.SScore as ss

UST.mfrMin = .01


def simpleaxis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


# dictParameter:     Tbin, MinWin,TH_end,Tend,thre,MinWin,mstd_flag
def th_plot(spike2n, parameterList):
    plt.ion()
    fig, ax = plt.subplots()

    plt.subplots_adjust(left=0.25, bottom=0.25)
    Ts = 100 * 1000
    rate = np.histogram(spike2n[1], bins=np.arange(0, Ts, parameterList['tbin (ms)']))
    ax.plot(rate[1][:-1], rate[0] * 100. / len(np.unique(spike2n[0])), 'b')
    TH_end = parameterList['TH_end (%)']
    TH_start = parameterList['TH_start (%)']
    Active = parameterList['N_ACTIVE (%)']
    plt.xlabel('Time (ms)')
    plt.ylabel('SpkCount per %s ms' % parameterList['tbin (ms)'])
    plots = list()

    plots.append(plt.plot([0, Ts], [TH_start, TH_start], 'g', lw=2, label='TH_start (%)')[0])
    plots.append(plt.plot([0, Ts], [TH_end, TH_end], 'r', lw=2, label='TH_end (%)')[0])
    plots.append(plt.plot([0, Ts], [Active, Active], 'c', lw=2, label='ACTIVE (%)')[0])
    plt.legend(loc='best')

    ax_ACTIVE = plt.axes([0.25, 0.12, 0.65, 0.02])
    ax_TH_start = plt.axes([0.25, 0.08, 0.65, 0.02])
    ax_TH_end = plt.axes([0.25, 0.05, 0.65, 0.02])

    SL_TH_start = Slider(ax_TH_start, 'TH_start (%)', 0, 20, valinit=parameterList['TH_start (%)'])
    SL_TH_end = Slider(ax_TH_end, 'TH_end (%)', 0, 20, valinit=parameterList['TH_end (%)'])
    SL_TH_Active = Slider(ax_ACTIVE, 'ACTIVE (%)', 0, 50, valinit=parameterList['N_ACTIVE (%)'])

    def update(val):
        TH_start = SL_TH_start.val
        TH_end = SL_TH_end.val
        Active = SL_TH_Active.val
        plots[0].set_ydata([TH_start, TH_start])
        plots[1].set_ydata([TH_end, TH_end])
        plots[2].set_ydata([Active, Active])
        fig.canvas.draw_idle()

    SL_TH_start.on_changed(update)
    SL_TH_end.on_changed(update)
    SL_TH_Active.on_changed(update)

    resetax = plt.axes([0.05, 0.5, 0.1, 0.04])
    button = Button(resetax, 'Reset', hovercolor='0.975')

    def reset(event):
        SL_TH_start.reset()
        SL_TH_end.reset()

    button.on_clicked(reset)
    plt.show(block=True)
    return SL_TH_start.val, SL_TH_end.val, SL_TH_Active.val


def run():
    folder = 'final4096NMDA_HD'
    fname = '4f_4096NMDA_small_tau_w_15_v4v3_37.5_13.9_sfinal_NMDA_00.npy'
    path2file = '/home/alonardoni/8192correlation/' + folder + '/Spike/' + fname
    Dfile = {'name': path2file}
    #######################
    ## LOADING
    plt.ioff()
    fileName = Dfile['name']
    Tstart = time.time()
    if fileName.rfind('.mat') > 0:
        parameterList = {'Chessboard clustered': True, 'Wind (ms)': 10, 'Save': False, 'Detection': 'Standard',
                         'Cut': 1.3, 'Clustering': 'CAT', 'BC tbin (ms)': 5, 'TH_end (%)': 1, 'Tend (s)': 10000,
                         'CATplot': 'All',
                         'TH_start (%)': 1, 'tbin (ms)': 10, 'MeanSTD': False, 'MinWind (ms)': 20, 'N_ACTIVE (%)': 5,
                         'CatWind (ms)': 10, 'ThPlot': True}
        spike2n = UST.FormatSpikeTrains(UST.LoadSpikeTrains(fileName))
        tmp = np.where((spike2n[1] < parameterList['Tend (s)'] * 1000) & (spike2n[1] > 0))[0]
        #        tmp=np.where((spike2n[1]<30*1000)&(spike2n[1]>0))[0]
        spike2n = np.asarray(spike2n)
        spike2n = spike2n[:, tmp]

        Dfile['_index'] = 0
        Dfile['analysisList'] = [[]]
        CellList = _convertClist(spike2n)


    elif fileName.rfind('.npy') > 0:
        parameterList = {'Chessboard clustered': True, 'Wind (ms)': 10, 'Save': False, 'Detection': 'Standard',
                         'Cut': 1.3, 'Clustering': 'CAT', 'BC tbin (ms)': 5, 'TH_end (%)': 1, 'Tend (s)': 10000,
                         'CATplot': 'All',
                         'TH_start (%)': 1, 'tbin (ms)': 10, 'MeanSTD': False, 'MinWind (ms)': 20, 'N_ACTIVE (%)': 5,
                         'CatWind (ms)': 10, 'ThPlot': True}
        CellList = np.load(fileName).tolist()

        spike2n = _convert2xN(CellList[1:])
        tmp = np.where((spike2n[1] < parameterList['Tend (s)'] * 1000) & (spike2n[1] > 0))[0]
        spike2n[0] = spike2n[0][tmp]
        spike2n[1] = spike2n[1][tmp]
        posx = np.zeros(4096)
        posy = np.zeros(4096)
        for cell in CellList[1:]:
            posx[cell['ID'][0]] = cell['x']
            posy[cell['ID'][0]] = cell['y']

        CellList = _convertClist(spike2n)
    spike2n = np.asarray(spike2n)
    print "TIME------------------>", time.time() - Tstart

    #######################
    ## DETECTION
    burst = []
    out = []
    if parameterList['Detection'].rfind('SMART') >= 0:
        out, G = SBED._detect_SMART(spike2n, Tbin=parameterList['tbin (ms)'], Tin=0,
                                    N_ACTIVE=parameterList['N_ACTIVE (%)'], Tdist=parameterList['Tdist (ms)'],
                                    DIST=parameterList['Max e2e dist'], MinWind=parameterList['MinWind (ms)'])
        IDs = set(range(4096))
        for ele in out:
            IDs = IDs.intersection(set(ele[0]))
        ibi_dist = []
        burst = []
        for ele in out:
            burst.append(_convertClist(np.asarray(ele)))
            ibi_dist.append(np.min(ele[1]))
    elif parameterList['Detection'].rfind('Standard') >= 0:
        if parameterList['ThPlot']:
            parameterList['TH_start (%)'], parameterList['TH_end (%)'], parameterList['N_ACTIVE (%)'] = th_plot(spike2n,
                                                                                                                parameterList)
        burst, ibi_dist = _detect(CellList, spike2n, Tbin=parameterList['tbin (ms)'],
                                  MinWin=parameterList['MinWind (ms)'], TH_end=parameterList['TH_end (%)'] * .01,
                                  TH_start=parameterList['TH_start (%)'] * .01,
                                  N_ACTIVE=parameterList['N_ACTIVE (%)'] * .01)
        out = []
        for ele in burst:
            out.append(_convert2xN(ele))
        IDs = set(range(4096))
        for ele in out:
            IDs = IDs.intersection(set(ele[0]))
    elif parameterList['Detection'].rfind('Reload') >= 0:
        burst = np.load(os.getcwd() + '/plugins/tmpFile/n' + fileName.split('/')[-1][:-5] + '_burst.npy')
        out = []
        ibi_dist = []
        for ele in burst:
            out.append(_convert2xN(ele))
            ibi_dist.append(np.min(ele[1]))
        IDs = set(range(4096))
        for ele in out:
            IDs = IDs.intersection(set(ele[0]))
    print "TIME------------------>", time.time() - Tstart
    #######################
    ## SAVE
    if parameterList['Save'] and parameterList['Detection'].rfind('Reload') < 0:
        np.save(open(os.getcwd() + '/plugins/tmpFile/n' + fileName.split('/')[-1][:-5] + '_burst.npy', 'w'), burst)

    listtrj, alisttrj, Tinit, elaps = _generate_trj(out, MinWin=parameterList['MinWind (ms)'], showplot=0,
                                                    pos=(posx, posy))

    #######################
    ## CLUSTERING
    mstd_flag = 0
    show_plot = 0
    if parameterList['CATplot'].rfind('All') >= 0:
        show_plot = 1
    if parameterList['CATplot'].rfind('Mean') >= 0:
        mstd_flag = 1
        show_plot = 1
    if parameterList['Clustering'].rfind('BC') >= 0:
        #######################

        lburst = []
        for ele in out:
            lburst.append(_convertClist(np.asarray(ele)))
        Dmap = fastburstcorrelation(lburst, tbin=parameterList['Precision (ms)'])
        del lburst
        print " "
        all_times = plot_cat_clustering(burst, listtrj, alisttrj, Tinit, Dmap,
                                        mstd_flag=mstd_flag, show_plot=show_plot,
                                        ChessboardC=parameterList["Chessboard clustered"])
    else:
        Dmap = _dist_mat(listtrj)
        all_times = plot_cat_clustering(burst, listtrj, alisttrj, Tinit, Dmap,
                                        mstd_flag=mstd_flag, show_plot=show_plot,
                                        ChessboardC=parameterList["Chessboard clustered"])
    #######################
    ## OUTPUT
    plist = [key for key in sorted(parameterList.keys())]
    vlist = [parameterList[key] for key in sorted(parameterList.keys())]
    Dfile['analysisList'] = dict()

    Dfile['analysisList']['CATClustering'] = np.array(
        [ele for item in [[plist, vlist], [[' ']], elaps, [[' ']], all_times] for ele in item])
    plt.show(1)
    if (parameterList['CATplot'].rfind('None') < 0 | parameterList['ThPlot']):
        plt.ion()
        plt.show()
    else:
        plt.close('all')
    return Dfile


def init(burst, wind):
    plt.ioff()
    parameterList = {'Cut': 3, 'BC': False, 'BC tbin (ms)': 5, 'TH_end (%)': 1,
                     'Tend (s)': 100, 'CATPlot': True, 'TH_start (%)': 1, 'tbin (ms)': 10, 'MeanSTD': False,
                     'MinWind (ms)': wind, 'N_ACTIVE (%)': 20, 'CatWind (ms)': 10, 'ThPlot': False}
    listtrj, alisttrj, Tinit, elaps = _generate_trj(burst, MinWin=parameterList['MinWind (ms)'], showplot=0)
    return listtrj, alisttrj, Tinit, elaps


def test(listtrj, alisttrj, Tinit, cut):
    parameterList = {'Cut': cut, 'BC': False, 'BC tbin (ms)': 5, 'TH_end (%)': 1,
                     'Tend (s)': 100, 'CATPlot': True, 'TH_start (%)': 1, 'tbin (ms)': 10, 'MeanSTD': False,
                     'MinWind (ms)': 20, 'N_ACTIVE (%)': 20, 'CatWind (ms)': 10, 'ThPlot': False}

    Dmap = _dist_mat(listtrj)
    burst = []
    plot_cat_clustering(burst, listtrj, alisttrj, Tinit, Dmap, mstd_flag=parameterList['MeanSTD'], show_plot=parameterList['CATPlot'])
    plt.ion()
    plt.show()


def _dist_mat(listtrj):
    N_trj = len(listtrj)
    Dmap = np.zeros((N_trj, N_trj))
    for i in xrange(N_trj):
        for j in xrange(i, N_trj):
            sys.stdout.write("\r Burst %04d: %04d of %04d" % (i, j, N_trj))
            sys.stdout.flush()
            l2Vect = abs(listtrj[i]['meanx'] - listtrj[j]['meanx']) ** 2 + abs(
                listtrj[i]['meany'] - listtrj[j]['meany']) ** 2
            #            l2Vect=l2Vect*(np.arange(1,len(l2Vect)+1)[::-1])*1./sum(np.arange(1,len(l2Vect)+1))
            Dmap[i, j] = np.mean(np.sqrt(l2Vect))
            Dmap[i, j] = np.max(np.sqrt(l2Vect))
            Dmap[i, j] = np.sqrt(
                np.sum(l2Vect * np.arange(len(l2Vect), 0, -1) ** 2) * 1. / np.sum(np.arange(len(l2Vect), 0, -1) ** 2))
            #            Dmap[i,j]=np.mean(l2Vect)
            Dmap[j, i] = Dmap[i, j]
    print " "
    return Dmap


def plot_cat_clustering(listtrj, alisttrj, Tinit, Dmap, thre=5, mstd_flag=0, show_plot=1,
                        ChessboardC=0):
    """
    compute clustering of propagations based on cat-distance:
    d(trj[i],trj[j])=sum of pointwise distances across 50 points of the interpolated (cubic-splines) CAT

    INPUT

    lburst:     list of CellList (one per SBE). CellList: list of dictionaries {'ID':int,'spikes':array,'x':[0,1],'y':[0,1]}
    thre:       dendrogram threshold

    CATWin:     CAT time-window to average over
    mstd_flag:  0 to have all the trj, 1 for mean trj +- std

    OUTPUT

    all_times:  per groups time_stamps of SBEs

    FIGURES:

    trajectories divided per group (size ordered)
    unordered/ordered Distance Matrix + box of trajectories
    raster o trj group size-ordered
    """

    Dmap, D0, composition, cut_index, pattern_mat = _reorder_Cmatrix(Dmap, thre, valid_Tinit=10 ** 10)
    _matrix_plot(Dmap, D0, composition, cut_index)
    plt.show(1)
    if show_plot:
        _matrix_plot(Dmap, D0, composition, cut_index)
    if ChessboardC:
        all_times = _plot_trj_Q(listtrj, composition, cut_index, Tinit, mstd_flag=mstd_flag, alisttrj=alisttrj,
                                show_plot=show_plot)
    else:
        all_times = _plot_trj(listtrj, composition, cut_index, Tinit, mstd_flag=mstd_flag, alisttrj=alisttrj,
                              show_plot=show_plot)
    return all_times


def plot_CC_clustering(lburst, BC, thre=2.8, MinWin=20., mstd_flag=0):
    """
    compute clustering of propagations based on BC-distance:
    d(trj[i],trj[j])=sum of pointwise distances across 50 points of the interpolated (cubic-splines) CAT

    INPUT

    lburst:     list of CellList (one per SBE). CellList: list of dictionaries {'ID':int,'spikes':array,'x':[0,1],'y':[0,1]}
    BC:         cross-correlation matrix or another distance matrix
    thre:       dendrogram threshold

    MinWin:     CAT time-window to average over
    mstd_flag:  0 to have all the trj, 1 for mean trj +- std

    OUTPUT

    all_times:  per groups time_stamps of SBEs

    FIGURES:

    trajectories divided per group (size ordered)
    unordered/ordered Distance Matrix + box of trajectories
    raster o trj group size-ordered
    """
    listtrj, alisttrj, Tinit, elaps = _generate_trj(lburst, MinWin)

    Dmap = BC
    Dmap, D0, composition, cut_index, pattern_mat = _reorder_Cmatrix(Dmap, thre, valid_Tinit=10 ** 10)
    _matrix_plot(Dmap, D0, composition, cut_index)
    all_times = _plot_trj_Q(listtrj, composition, cut_index, Tinit, mstd_flag=mstd_flag, alisttrj=alisttrj)
    return all_times


def _generate_trj(lburst, MinWin, showplot=0, pos=(), dirName='TESTfigALL/', figname='test'):
    Tbin = 1
    listtrj = []
    alisttrj = []
    Tinit = np.zeros(len(lburst))
    elaps = [['Time stamp of NB (ms)'], ['Duration (ms)'], ['Spk count'], ['IBI'], ['EChannels']]
    elaps[3].append(0)
    posx = pos[0]
    posy = pos[1]


    for ii, burst in enumerate(lburst):
        if np.max(burst[1]) - np.min(burst[0]) > MinWin:
            burst = np.asarray(burst)
            Tinit[ii] = np.min(burst[1])
            tstart = np.floor(Tinit[ii])
            tend = np.floor(max(burst[1]) + 1)

            # STATS to export into .xlsx file
            elaps[0].append(tstart)
            elaps[1].append(tend - tstart)
            elaps[2].append(len(burst[0]))
            if len(elaps[0]) > 2:
                elaps[3].append(elaps[0][-1] - elaps[0][-2])
            elaps[4].append(len(np.unique(burst[0])))

            t = np.arange(tstart, tend, Tbin)
            mx = np.zeros(int(np.floor(tend - tstart) / Tbin))
            my = np.zeros(int(np.floor(tend - tstart) / Tbin))
            lmx = len(mx)
            colorRasterX = []
            colorRasterY = []
            if lmx > 5:
                for ind_D, d in enumerate(np.arange(tstart - MinWin / 2, tend - MinWin / 2, Tbin)):
                    indices = np.where((burst[1] > d) & (burst[1] < d + MinWin))[0]
                    if indices.any():
                        try:
                            mx[ind_D] = np.mean(posx[burst[0][indices.astype(int)].astype(int)])
                            my[ind_D] = np.mean(posy[burst[0][indices.astype(int)].astype(int)])
                            colorRasterX.append(posx[burst[0][indices.astype(int)].astype(int)])
                            colorRasterY.append(posy[burst[0][indices.astype(int)].astype(int)])
                        except:
                            print "Problem Reported:"
                            print "empty slice?", indices.astype(int)
                            print "empty slice?", len(burst[0]), len(burst[1])
                            print "empty slice?", burst[0][indices.astype(int)].astype(int)
                            print "empty slice?", posx[burst[0][indices.astype(int)].astype(int)]



                    else:
                        mx[ind_D] = mx[ind_D - 1]
                        my[ind_D] = my[ind_D - 1]

                graph_ref = 0
                if graph_ref:
                    txx = np.linspace(.08, .5, len(mx))
                    sg = np.arange(4096)
                    for sm in np.arange(int(len(mx) / 2))[::-1]:
                        NSp = len(colorRasterX[sm])
                        SpX = colorRasterX[sm]
                        SpY = colorRasterY[sm]
                        test = np.zeros((NSp, NSp))
                        for i in xrange(NSp):
                            for j in xrange(i + 1, NSp):
                                if abs(SpX[i] - SpX[j]) < txx[sm] and abs(SpY[i] - SpY[j]) < txx[sm]:
                                    test[i, j] = 1
                                    test[j, i] = 1
                        G = nx.from_numpy_matrix(test, create_using=nx.Graph())
                        sgS = np.array(sorted([subg for subg in nx.connected_components(G)], key=len))[::-1]

                        for ele in sgS:
                            #                    print  "\t - \t - \t ",np.intersect1d(sgS,sg)
                            if any(np.intersect1d(ele, sg)):
                                sg = ele
                                smS = sm
                                break
                                #                print sm,smS, '   ---->   ',np.mean(colorRasterX[smS]),np.mean(colorRasterX[smS][sg])
                        mx[sm] = np.mean(colorRasterX[smS][sg])
                        my[sm] = np.mean(colorRasterY[smS][sg])

                if showplot:
                    plt.plot(mx, my, 'b')
                intx = interp1d(np.arange(0, tend - tstart, Tbin)[:lmx], mx, kind='cubic')
                inty = interp1d(np.arange(0, tend - tstart, Tbin)[:lmx], my, kind='cubic')
                if showplot:
                    plt.plot(intx(np.linspace(0, t[lmx - 1] - t[0], 50)), inty(np.linspace(0, t[lmx - 1] - t[0], 50)),
                             'r')

                listtrj.append({'meanx': intx(np.linspace(0, t[lmx - 1] - t[0], 100)),
                                'meany': inty(np.linspace(0, t[lmx - 1] - t[0], 100))})
                alisttrj.append({'meanx': mx, 'meany': my})
                if len(alisttrj) < 0:
                    plt.ioff()
                    fig = plt.figure(figsize=(18, 9))
                    ax1 = fig.add_axes([0.05, 0.05, 0.4, 0.4])
                    ax2 = fig.add_axes([0.05, 0.5, 0.4, 0.4])
                    #                plt.subplot(121)
                    #            plt.subplot2grid((2,2), (0, 0), colspan=1)
                    ax1.plot(burst[1], posy[burst[0].astype(int)], 'ok', ms=2)
                    ax1.plot(np.linspace(t[0], t[lmx - 1], 100), inty(np.linspace(0, t[lmx - 1] - t[0], 100)), 'r')
                    ax1.plot([tstart, tstart], [0, 1], lw=2)
                    ax1.plot([tstart + MinWin, tstart + MinWin], [0, 1], lw=2)
                    ax1.locator_params(nbins=4)

                    ax2.plot(burst[1], posx[burst[0].astype(int)], 'ok', ms=2)
                    ax2.plot(np.linspace(t[0], t[lmx - 1], 100), intx(np.linspace(0, t[lmx - 1] - t[0], 100)), 'r')
                    ax2.locator_params(nbins=4)

                    axM = fig.add_axes([0.5, 0.05, 0.4, 0.85])
                    axC = fig.add_axes([0.93, 0.05, 0.02, 0.85])
                    NCol = len(colorRasterX)
                    already_plotted = set([])
                    colors2use = []
                    for ncol in xrange(NCol):
                        i2plot = np.asarray([i for i in xrange(len(colorRasterX[ncol])) if
                                             not (colorRasterX[ncol][i], colorRasterY[ncol][i]) in already_plotted])
                        colors2use.append(len(i2plot))
                        already_plotted = already_plotted.union(
                            set([(x, y) for x, y in zip(colorRasterX[ncol], colorRasterY[ncol])]))
                    colors2use = np.max([i for i in xrange(len(colors2use)) if colors2use[i] > 5])
                    if colors2use == 0:
                        colors2use = 1
                    already_plotted = set([])
                    for ncol in xrange(NCol):
                        i2plot = np.asarray([i for i in xrange(len(colorRasterX[ncol])) if
                                             not (colorRasterX[ncol][i], colorRasterY[ncol][i]) in already_plotted])

                        if any(i2plot):
                            if ncol == 0:
                                axM.plot(colorRasterX[ncol][i2plot], colorRasterY[ncol][i2plot], 'o', color='k',
                                         zorder=-ncol, rasterized=True)
                            else:
                                axM.plot(colorRasterX[ncol][i2plot], colorRasterY[ncol][i2plot], 'o',
                                         color=plt.get_cmap('cool')(min(ncol * 1. / colors2use, .99)), zorder=-ncol,
                                         rasterized=True)
                        already_plotted = already_plotted.union(
                            set([(x, y) for x, y in zip(colorRasterX[ncol], colorRasterY[ncol])]))

                    axM.plot(mx, my, 'r', lw=4)
                    axM.plot(mx[0], my[0], 'or', ms=10)
                    import matplotlib as mpl

                    simpleaxis(axM)

                    plt.suptitle('ActiveCh: %04d' % len(np.where(posx > 0)))
                    plt.savefig(dirName + '/propagations/' + figname + '_%06d.png' % ii)

                plt.close('all')
    return listtrj, alisttrj, Tinit, elaps


def get_ISs(fileName, figname='', destPath='TESTfigALL/'):
    parameterList = {'Chessboard clustered': True, 'Wind (ms)': 10, 'Save': False, 'Detection': 'Standard',
                     'Cut': 1.3, 'Clustering': 'CAT', 'BC tbin (ms)': 5, 'TH_end (%)': 5, 'Tend (s)': 10000,
                     'CATplot': 'All',
                     'TH_start (%)': .1, 'tbin (ms)': 10, 'MeanSTD': False, 'MinWind (ms)': 20, 'N_ACTIVE (%)': 10,
                     'CatWind (ms)': 10, 'ThPlot': True}

    if fileName.rfind('.mat') > 0:
        CellList = UST.LoadSpikeTrains(fileName)
        spike2n = UST.FormatSpikeTrains(CellList)
        posx = np.zeros(4097)
        posy = np.zeros(4097)
        for cell in CellList:
            posx[cell['ID'][0]] = cell['x'] * 1. / 64
            posy[cell['ID'][0]] = cell['y'] * 1. / 64


    elif fileName.rfind('.npy') > 0:
        CellList = np.load(fileName).tolist()
        spike2n = _convert2xN(CellList[1:])
        tmp = np.where((spike2n[1] < parameterList['Tend (s)'] * 1000) & (spike2n[1] > 0))[0]
        spike2n[0] = spike2n[0][tmp]
        spike2n[1] = spike2n[1][tmp]
        posx = np.zeros(4096)
        posy = np.zeros(4096)
        for cell in CellList[1:]:
            posx[cell['ID'][0]] = cell['x']
            posy[cell['ID'][0]] = cell['y']

    burst, ibi_dist = _detect(CellList, spike2n, Tbin=parameterList['tbin (ms)'], MinWin=parameterList['MinWind (ms)'],
                              TH_end=parameterList['TH_end (%)'] * .01, TH_start=parameterList['TH_start (%)'] * .01,
                              N_ACTIVE=parameterList['N_ACTIVE (%)'] * .01)

    plt.close('all')
    out = []
    for ele in burst:
        out.append(_convert2xN(ele))
    if figname == '':
        figname = fileName.split('alonardoni/')[-1]
        figname = figname.replace('/', '_')
        figname = figname.replace('.', '_')
    if fileName.rfind('.npy') > 0:
        listtrj, alisttrj, Tinit, elaps = _generate_trj(out, MinWin=30, showplot=0, pos=(posx, posy), dirName=destPath,
                                                        figname=figname)
    else:
        listtrj, alisttrj, Tinit, elaps = _generate_trj(out, MinWin=30, showplot=0, pos=(posx, posy), dirName=destPath,
                                                        figname=figname)

    plt.ioff()
    Dmap = _dist_mat(listtrj)
    Dmap, D0, composition, cut_index, pattern_mat = _reorder_Cmatrix(Dmap, ss.Iterate_HierClust(Dmap, min(20, len(
        burst) - 3))[-1], valid_Tinit=10 ** 10)
    _plot_trj(listtrj, composition, cut_index, Tinit, mstd_flag=0, alisttrj=alisttrj, show_plot=1)
    figures = [manager.canvas.figure for manager in Gcf.get_all_fig_managers()]

    for i, figure in enumerate(figures):
        figure.savefig(destPath + '/summary_' + figname + '_6M_2_%02d.svg' % i)
    plt.close('all')
    return [(ele['meanx'], ele['meany']) for ele in listtrj], cut_index


def _reorder_Cmatrix(D, thre, valid_Tinit=10 ** 10, showplot=0):
    Y = sch.linkage(squareform(D), method='complete')
    Z = sch.dendrogram(Y, no_plot=True)

    index = np.array(Z['leaves'])
    cut_index = sch.fcluster(Y, t=thre, criterion='distance')
    composition = np.array([-sum(cut_index[:valid_Tinit] == i) for i in np.unique(cut_index)])
    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    neworder = []
    group_index = []
    for i in xrange(len(ranks)):
        pivot = ranks.argsort()[i]
        indices = [index[ele] for ele in xrange(len(index)) if cut_index[index[ele]] == pivot + 1]
        group_index.append(indices)
        neworder.extend(indices)
    neworder = np.array(neworder)  # [::-1]

    D0 = D[neworder, :]
    tmp = neworder  # [::-1]
    D0 = D0[:, tmp]

    test = D * 0
    for i in xrange(len(group_index)):
        a = group_index[i]
        b = group_index[i]
        plist = [(x0, y0) for x0 in a for y0 in b]
        for ele in plist: test[ele] = i + 1

    pattern_mat = -np.ones((10, len(test)))
    tmp = test.diagonal(0)
    pattern_mat[0, :] = tmp
    for i in xrange(1, 10):
        pattern_mat[i, i:] = tmp[:-i]
    if showplot:
        plt.matshow(pattern_mat)
        plt.yticks([9, 7, 5, 3, 1], ['SBE', '+2', '+4', '+6', '+8'])

    return D, D0, composition, cut_index[:valid_Tinit], pattern_mat


def _plot_mat(a, composition, cut_index):

    N_groups = len(np.unique(cut_index))
    N_trj = int(np.ceil(np.sqrt(N_groups)))
    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    plt.figure()

    for kk in np.unique(cut_index):
        i = ranks[kk - 1]
        plt.subplot(np.ceil(N_groups * 1. / N_trj), N_trj, i + 1)

        plt.text(0, 0, '%s' % i, fontsize=40, zorder=-1, alpha=.5)

        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        plt.hold(True)

        indices = np.array([j for j in np.arange(10000)[cut_index == kk]])

        ax.matshow(a[:, indices], vmin=0, vmax=6)


def _generate_mat(a, composition, cut_index):
    l_matrix = []
    temp = composition.argsort()
    for kk in np.unique(cut_index):


        indices = np.array([j for j in np.arange(10000)[cut_index == kk]])

        l_matrix.append(a[:, indices])
    return l_matrix


def _matrix_plot(D, D0, composition, cut_index):
    from matplotlib.patches import Rectangle
    cmap_s = cm.gray
    N_groups = len(np.unique(cut_index))
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].matshow(D, aspect='auto', origin='lower', vmin=0, cmap=cmap_s)

    im = axes[1].matshow(D0, aspect='auto', origin='lower', vmin=0, cmap=cmap_s)
    count = 0
    comp = sorted(composition)
    test = 0
    for i in xrange(len(comp)):
        lindex = -comp[i]
        axes[1].add_patch(Rectangle((count - .5, count - .5), lindex, lindex, fill=None, alpha=1,
                                    color=plt.get_cmap('jet')(float(i) / N_groups), lw=4))
        count += lindex
        test += lindex
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    fig.suptitle('Distance Matrices')
    axes[0].set_xlabel('UNORDERED')
    axes[0].xaxis.set_ticks_position('bottom')
    axes[1].set_xlabel('ORDERED')
    axes[1].xaxis.set_ticks_position('bottom')


def _plot_trj(listtrj, composition, cut_index, Tinit, mstd_flag=0, alisttrj=[], show_plot=1):

    N_groups = len(np.unique(cut_index))
    N_trj = int(np.ceil(np.sqrt(N_groups)))
    all_times = []
    count_b = 0
    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    if show_plot:
        plt.figure()
    for kk in np.unique(cut_index):

        i = ranks[kk - 1]
        if show_plot:
            count = 0
            tx_mean = np.zeros((sum(cut_index == kk), 20))
            ty_mean = np.zeros((sum(cut_index == kk), 20))
            for j in np.arange(10000)[cut_index == kk]:
                tx_mean[count, :] = listtrj[j]['meanx'][np.arange(0, 100, 5)]
                ty_mean[count, :] = listtrj[j]['meany'][np.arange(0, 100, 5)]
                count += 1
            count_b += count
            x0 = tx_mean.mean(0)[0]
            y0 = ty_mean.mean(0)[0]

            plt.subplot(np.ceil(N_groups * 1. / N_trj), N_trj, i + 1)
            cgr = plt.get_cmap('jet')(float(i) / N_groups)
            plt.text(.5, .5, '%s' % i, fontsize=40, zorder=-1, alpha=.5, va='center', ha='center')
            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.hold(True)

            def atan(x, y):
                # atan yields angles on the whole [0,2*pi] interval
                theta = np.arctan(y / x) + np.pi / 2
                theta[x < 0] += np.pi
                return theta

            alpha = atan(np.diff(tx_mean.mean(0)), np.diff(ty_mean.mean(
                0)))  # np.pi/2+np.pi/2-np.arctan(np.divide(np.diff(ty_mean.mean(0)),np.diff(tx_mean.mean(0))))
            dMx = np.zeros(len(tx_mean.mean(0)))
            dMy = np.zeros(len(tx_mean.mean(0)))
            dMx[1:] = np.diff(tx_mean.mean(0))
            dMy[1:] = np.diff(ty_mean.mean(0))
            for iii in xrange(1, len(dMx)):
                norm = np.sqrt((dMx[iii] ** 2 + dMy[iii] ** 2))

                dMx[iii] = dMx[iii] * 1. / norm
                dMy[iii] = dMy[iii] * 1. / norm

            tx_meanR = tx_mean * -dMy + ty_mean * dMx

            stdx = [0]
            stdy = [0]

            std = tx_meanR.std(0)  # np.ones(len(tx_mean.std(0)))*.2
            for k in xrange(19):
                stdx.append(std[k + 1] * np.cos(alpha[k]))  # np.cos(alpha[k])*tx_std[k+1]-np.sin(alpha[k])*ty_std[k+1]
                stdy.append(std[k + 1] * np.sin(alpha[k]))  # np.sin(alpha[k])*ty_std[k+1]+np.cos(alpha[k])*ty_std[k+1]

            from pylab import Polygon
            if mstd_flag:
                import matplotlib.colors as colors
                ax = plt.gca()
                #


                pVerts1 = [[u, v] for u, v in zip(tx_mean.mean(0)[:] - stdx[:], ty_mean.mean(0)[:] - stdy[:])]
                pVerts2 = [[u, v] for u, v in zip(tx_mean.mean(0)[:] + stdx[:], ty_mean.mean(0)[:] + stdy[:])]
                plt.plot(tx_mean.mean(0), ty_mean.mean(0), lw=3, color=cgr)
                plt.scatter(tx_mean.mean(0)[0], ty_mean.mean(0)[0], c=cgr, zorder=10000)
                plt.text(x0 + .1, y0 + .1, '%s' % i, fontsize=15, alpha=.5, va='center', ha='center', color='k',
                         zorder=100)
                from matplotlib.collections import PolyCollection
                for i in range(len(pVerts1) - 1):
                    vtx = [[pVerts1[i + 1], pVerts1[i], pVerts2[i], pVerts2[i + 1]]]
                    tri = PolyCollection(vtx)
                    tri.set_color(cgr)
                    tri.set_alpha(.2)
                    tri.set_edgecolor('none')
                    tri.set_zorder(-i)
                    ax.add_collection(tri)

                # ax.plot(tx_mean.mean(0)[:]-stdx[:],ty_mean.mean(0)[:]-stdy[:],'--',color=cgr)
                #                ax.plot(tx_mean.mean(0)[:]+stdx[:],ty_mean.mean(0)[:]+stdy[:],'--',color=cgr)
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])

            # if mstd_flag and len(alisttrj)>0:
            #                plt.plot(tx_mean.mean(0),ty_mean.mean(0),lw=2,color=cgr)
            #                plt.plot(tx_mean.mean(0)[0],ty_mean.mean(0)[0],'o',color=cgr)
            #
            #                plt.plot(tx_mean.mean(0)[1:]-stdx[1:],ty_mean.mean(0)[1:]+stdy[1:],'--',color=cgr)
            #                plt.plot(tx_mean.mean(0)[1:]+stdx[1:],ty_mean.mean(0)[1:]-stdy[1:],'--',color=cgr)
            else:
                for j in np.arange(10000)[cut_index == kk]:
                    plt.plot(alisttrj[j]['meanx'], alisttrj[j]['meany'], lw=2, color=cgr)
                    plt.plot(alisttrj[j]['meanx'][0], alisttrj[j]['meany'][0], 'o', color=cgr)
                    #                    plt.plot(tx_mean.mean(0)[1:]-stdx[1:],ty_mean.mean(0)[1:]+stdy[1:],'--',color='k')
                    #                    plt.plot(tx_mean.mean(0)[1:]+stdx[1:],ty_mean.mean(0)[1:]-stdy[1:],'--',color='k')
            plt.xlim([0, 1])
            plt.ylim([0, 1])

        all_times.append(Tinit[cut_index == kk])
    if show_plot:
        plt.suptitle('Clustered Trajectories')
        fig = plt.figure(figsize=(8, 8))
        fig.suptitle('Raster of Clustered Trajectories')
        ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.8])
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])

    composition = np.array([-sum(all_times[i] < 60 * 60 * 10000) for i in xrange(len(all_times))])
    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    bp = []
    clusterTS = []
    for j in xrange(len(all_times)):
        a = all_times[j]

        i = ranks[j]
        if show_plot:
            cgr = plt.get_cmap('jet')(float(i) / N_groups)
            axmatrix.plot(a * 1. / 1000, i * np.ones(len(a)), 'o', color=cgr)
            ax1.barh(i - 0.4, len(a) * 1. / len(listtrj), color=cgr)

            if len(a) <= 1:
                axmatrix.plot(a * 1. / 1000, i * np.ones(len(a)), 'o', color=cgr)
                ax1.barh(i - 0.4, len(a) * 1. / len(listtrj), color=cgr)
                ax1.text(i - 0.4, len(a) * 1. / len(listtrj) * .5, '%s' % i, fontsize=20, alpha=.5, va='center',
                         ha='center')
            else:
                bp.append(axmatrix.boxplot(a * 1. / 1000, vert=False, positions=[i], widths=[0.8]))
                for box in bp[-1]['boxes']:
                    box.set(color=cgr, lw=2)
                for box in bp[-1]['whiskers']:
                    box.set(color=cgr, lw=2)
                for box in bp[-1]['medians']:
                    box.set(color='k', lw=3)
                for box in bp[-1]['caps']:
                    box.set(color=cgr, lw=2)
        clusterTS.append(a.tolist())
        clusterTS[-1].insert(0, i)
    if show_plot:
        ax1.set_ylim([-1, N_groups])
        ax1.locator_params(axis='x', nbins=4)
        ax1.invert_xaxis()
        axmatrix.set_ylim([-1, N_groups])
        axmatrix.set_yticks([])
    clusterTS = sorted(clusterTS, key=lambda x: x[0])
    clusterTS.insert(0, ['Cluster #', 'Time of event (ms)'])
    return clusterTS


def _plot_trj_Q(listtrj, composition, cut_index, Tinit, mstd_flag=0, alisttrj=[], show_plot=1):
    plt.close('all')
    N_groups = len(np.unique(cut_index))
    all_times = []
    count_b = 0
    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    groups = dict()
    if show_plot:
        plt.figure()
    for kk in np.unique(cut_index):
        i = ranks[kk - 1]

        if show_plot:
            count = 0
            tx_mean = np.zeros((sum(cut_index == kk), 20))
            ty_mean = np.zeros((sum(cut_index == kk), 20))
            for j in np.arange(10000)[cut_index == kk]:
                tx_mean[count, :] = listtrj[j]['meanx'][np.arange(0, 100, 5)]
                ty_mean[count, :] = listtrj[j]['meany'][np.arange(0, 100, 5)]
                count += 1
            count_b += count
            x0 = tx_mean.mean(0)[0]
            y0 = ty_mean.mean(0)[0]
            matrix = np.array([0, .33, 0.66])

            plt.subplot(3, 3, np.searchsorted(matrix, x0) + (3 - np.searchsorted(matrix, y0)) * 3)
            #            plt.subplot(np.ceil(N_groups*1./N_trj),N_trj,i+1)#,projection='3d')
            cgr = plt.get_cmap('afmhot')((float(i) / N_groups) * .5)
            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks([])
            ax.yaxis.set_ticks([])
            plt.hold(True)

            groups[i] = [tx_mean, ty_mean]

            def atan(x, y):
                # atan yields angles on the whole [0,2*pi] interval
                theta = np.arctan(y / x) + np.pi / 2
                theta[x < 0] += np.pi
                return theta

            alpha = atan(np.diff(tx_mean.mean(0)), np.diff(ty_mean.mean(
                0)))  # np.pi/2+np.pi/2-np.arctan(np.divide(np.diff(ty_mean.mean(0)),np.diff(tx_mean.mean(0))))
            dMx = np.zeros(len(tx_mean.mean(0)))
            dMy = np.zeros(len(tx_mean.mean(0)))
            dMx[1:] = np.diff(tx_mean.mean(0))
            dMy[1:] = np.diff(ty_mean.mean(0))
            for iii in xrange(1, len(dMx)):
                norm = np.sqrt((dMx[iii] ** 2 + dMy[iii] ** 2))

                dMx[iii] = dMx[iii] * 1. / norm
                dMy[iii] = dMy[iii] * 1. / norm

            tx_meanR = tx_mean * -dMy + ty_mean * dMx
            stdx = [0]
            stdy = [0]

            std = tx_meanR.std(0)  # np.ones(len(tx_mean.std(0)))*.2
            for k in xrange(19):
                stdx.append(std[k + 1] * np.cos(alpha[k]))  # np.cos(alpha[k])*tx_std[k+1]-np.sin(alpha[k])*ty_std[k+1]
                stdy.append(std[k + 1] * np.sin(alpha[k]))  # np.sin(alpha[k])*ty_std[k+1]+np.cos(alpha[k])*ty_std[k+1]

            from pylab import Polygon
            if mstd_flag > 0:

                import matplotlib.colors as colors
                ax = plt.gca()
                #


                pVerts1 = [[u, v] for u, v in zip(tx_mean.mean(0)[:] - stdx[:], ty_mean.mean(0)[:] - stdy[:])]
                pVerts2 = [[u, v] for u, v in zip(tx_mean.mean(0)[:] + stdx[:], ty_mean.mean(0)[:] + stdy[:])]
                plt.plot(tx_mean.mean(0), ty_mean.mean(0), lw=3, color=cgr)
                plt.scatter(tx_mean.mean(0)[0], ty_mean.mean(0)[0], c=cgr, zorder=10000)
                plt.text(x0 + .1, y0 + .1, '%s' % i, fontsize=15, alpha=.5, va='center', ha='center', color='k',
                         zorder=100)
                from matplotlib.collections import PolyCollection
                for i in range(len(pVerts1) - 1):
                    vtx = [[pVerts1[i + 1], pVerts1[i], pVerts2[i], pVerts2[i + 1]]]
                    tri = PolyCollection(vtx)
                    tri.set_color(cgr)
                    tri.set_alpha(.2)
                    tri.set_edgecolor('none')
                    tri.set_zorder(-i)
                    ax.add_collection(tri)

                # ax.plot(tx_mean.mean(0)[:]-stdx[:],ty_mean.mean(0)[:]-stdy[:],'--',color=cgr)
                #                ax.plot(tx_mean.mean(0)[:]+stdx[:],ty_mean.mean(0)[:]+stdy[:],'--',color=cgr)
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
            else:
                for j in np.arange(10000)[cut_index == kk]:
                    #                    plt.plot(alisttrj[j]['meanx'],alisttrj[j]['meany'],lw=2,color=cgr)
                    plt.plot(alisttrj[j]['meanx'][0], alisttrj[j]['meany'][0], 'o', color=cgr)
                # plt.plot(tx_mean.mean(0)[1:]-stdx[1:],ty_mean.mean(0)[1:]+stdy[1:],'--',color='k')
                #                    plt.plot(tx_mean.mean(0)[1:]+stdx[1:],ty_mean.mean(0)[1:]-stdy[1:],'--',color='k')
                ax.set_xlim([0, 1])
                ax.set_ylim([0, 1])
        all_times.append(Tinit[cut_index == kk])
    if show_plot:
        plt.suptitle('Clustered Trajectories')
        fig = plt.figure(figsize=(8, 8))
        fig.suptitle('Raster of Clustered Trajectories')
        ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.8])
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])

    temp = composition.argsort()
    ranks = np.arange(len(composition))[temp.argsort()]
    bp = []
    clusterTS = []
    for j in xrange(len(all_times)):
        a = all_times[j]

        i = ranks[j]
        if show_plot:
            cgr = plt.get_cmap('jet')(float(i) / N_groups)
            axmatrix.plot(a * 1. / 1000, i * np.ones(len(a)), 'o', color=cgr)
            ax1.barh(i - 0.4, len(a) * 1. / len(listtrj), color=cgr)

            if len(a) <= 1:
                axmatrix.plot(a * 1. / 1000, i * np.ones(len(a)), 'o', color=cgr)
                ax1.barh(i - 0.4, len(a) * 1. / len(listtrj), color=cgr)
                ax1.text(i - 0.4, len(a) * 1. / len(listtrj) * .5, '%s' % i, fontsize=20, alpha=.5, va='center',
                         ha='center')
            else:
                bp.append(axmatrix.boxplot(a * 1. / 1000, vert=False, positions=[i], widths=[0.8]))
                for box in bp[-1]['boxes']:
                    box.set(color=cgr, lw=2)
                for box in bp[-1]['whiskers']:
                    box.set(color=cgr, lw=2)
                for box in bp[-1]['medians']:
                    box.set(color='k', lw=3)
                for box in bp[-1]['caps']:
                    box.set(color=cgr, lw=2)
        clusterTS.append(a.tolist())
        clusterTS[-1].insert(0, i)
    if show_plot:
        ax1.set_ylim([-1, N_groups])
        ax1.locator_params(axis='x', nbins=4)
        ax1.invert_xaxis()
        axmatrix.set_ylim([-1, N_groups])
        axmatrix.set_yticks([])
    clusterTS = sorted(clusterTS, key=lambda x: x[0])
    clusterTS.insert(0, ['Cluster #', 'Time of event (ms)'])
    return clusterTS


def _convert2xN(CellList):
    tmpid = []
    tmptime = []
    for cell in CellList:
        try:
            tmpid.extend(np.ones(len(cell['spikes'])) * cell['ID'].astype(int))
        except:
            tmpid.extend(np.ones(len(cell['spikes'])) * cell['ID'])
        tmptime.extend(cell['spikes'])
    tmp = np.array([np.transpose(np.array(tmpid)), np.transpose(np.array(tmptime))])
    return tmp


def _convertClist(list2xN):
    CellList = []
    for ele in np.unique(list2xN[0]):
        CellList.append({'ID': ele, 'spikes': list2xN[1][list2xN[0] == ele], 'x': (int(ele - 1) / 64) * 1. / 64,
                         'y': (int(ele - 1) % 64) * 1. / 64})

    return CellList


def _detect(CellList, spikes, Tbin=10, MinWin=50.0, TH_end=.1, TH_start=.002, Tin=0, N_ACTIVE=.05):
    Tend = np.max(spikes[1])
    rate = np.histogram(spikes[1], bins=np.arange(Tin, Tend, Tbin))[0]

    # Extra parameter definition

    N = len(np.unique(spikes[0]))
    Twind = int(MinWin / Tbin)

    burst = []
    ibi_dist = []

    ID = 0
    end = 0

    sbN = 0
    while (end + Twind + Tin) * Tbin < Tend - (Twind * 3) * Tbin:
        print (end + Twind + Tin) * Tbin
        start = np.argmax(rate[end:] > TH_start * N) + end
        end = np.argmax(rate[start:] < TH_end * N) + start
        endOld = end
        if end - start > Twind and max(rate[start:end]) > N * N_ACTIVE:

            if (end - start) * Tbin > 100:
                end = int(100 / Tbin) + 1 + start
            kk = 0
            sbN += 1
            brtmp = []
            for h in xrange(1, N):
                pivotstart = np.searchsorted(CellList[h]['spikes'], (start - 1) * Tbin, 'left')
                pivotend = pivotstart + np.searchsorted(CellList[h]['spikes'][pivotstart:], (end + 1) * Tbin, 'left')
                times = CellList[h].copy()
                times['spikes'] = CellList[h]['spikes'][pivotstart:pivotend]
                if len(times['spikes']) > 0:
                    kk += 1
                brtmp.append(times.copy())

            burst.append(brtmp)
            ibi_dist.append(start * Tbin)
            ID += 1
        end = endOld
        end += 1

    INspk_list = [item for ele in burst for cell in ele for item in cell['spikes']]
    rnd_Spk = 1 - len(INspk_list) * 1. / len(spikes[1])
    print "True Random Spikes:", rnd_Spk * 100

    return burst, ibi_dist
