import numpy as np
import pylab as plt
from sklearn.metrics import silhouette_score
import scipy.cluster.hierarchy as hac
from scipy.spatial.distance import squareform
import sys

sys.setrecursionlimit(10000)


def flatlist(x):
    return [ele for item in x for ele in item]


def Iterate_HierClust(FeatureMatrix, nBranches, method='complete', criterion='distance', metric='precomputed'):
    """
        Iterates the "Hierachical clustering" algorithm
        FeatureMatrix       the usual FeatureMatrix (last column is a label)
        nBranches           number of branches to consider
    """
    try:
        zw = hac.linkage(squareform(FeatureMatrix), method=method)
    except:
        return [0, 0, 10 ** 10]
    d = hac.dendrogram(zw)
    AllTH = np.unique(flatlist(d['dcoord']))  # d['dcoord'] returns all needed thresholds!
    dMIN = np.min(np.diff(AllTH)) / 2.0
    plt.hold(1)

    # prepare variables for the output
    LabelsMat = np.zeros((FeatureMatrix.shape[0], nBranches, 1))
    ClustInfo = {'nclust': [], 'params': {}}
    ClustInfo['params']['cut-value'] = []
    # + other info ?
    plt.close('all')
    FmaxL = []
    for ind, th in enumerate(AllTH[-nBranches:]):
        th2use = th - dMIN
        labels = hac.fcluster(zw, t=th2use, criterion=criterion)
        L = len(set(labels))
        sys.stdout.write("\r RUN# %02d out of %d -> #clusters=%02d" % (ind + 1, nBranches, L))
        sys.stdout.flush()
        #        print labels
        LabelsMat[:, ind, 0] = labels - 1
        ClustInfo['nclust'].append(L)
        ClustInfo['params']['cut-value'].append(th2use)
        tmp = silhouette_score(FeatureMatrix, labels - 1, metric=metric)
        FmaxL.append([tmp, L, th])
    print(" ")
    Fmax = sorted(FmaxL, key=lambda x: x[0])[-1]
    Lind = dict()
    Lind[0] = FmaxL.index(Fmax)

    data2return = sorted(FmaxL, key=lambda x: x[0])[-10]
    print(data2return)
    data2return[-1] = .5 * (FmaxL[Lind[0]][-1] + FmaxL[Lind[0] - 1][-1])
    print(data2return)
    return data2return

def line_picker(line, mouseevent):
    """
    find the points within a certain distance from the mouseclick in
    data coords and attach some extra attributes, pickx and picky
    which are the data points that were picked
    """
    if mouseevent.xdata is None: return False, dict()
    xdata = line.get_xdata()
    ydata = line.get_ydata()
    maxd = 0.05
    d = np.sqrt((xdata - mouseevent.xdata) ** 2. + (ydata - mouseevent.ydata) ** 2.)

    ind = np.nonzero(np.less_equal(d, maxd))
    if len(ind):
        pickx = np.take(xdata, ind)
        picky = np.take(ydata, ind)
        props = dict(ind=ind, pickx=pickx, picky=picky)
        return True, props
    else:
        return False, dict()


