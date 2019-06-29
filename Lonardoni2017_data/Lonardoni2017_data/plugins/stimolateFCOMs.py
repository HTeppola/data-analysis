import numpy as np
import pylab as plt
import os
import tools.tools as t
import pickle
import networkx as nx

Rx = .5
Ry = .5


def isRandomST(STspkX, STspkY, i, Rx=Rx, Ry=Ry):
    return abs(
        np.mean(STspkX[0, np.where((STspkX[1] > i * 1000) & (STspkX[1] < i * 1000 + 500))[0]]) - Rx) < .1 and abs(
        np.mean(STspkY[0, np.where((STspkY[1] > i * 1000) & (STspkY[1] < i * 1000 + 500))[0]]) - Ry) < .1


def isRandom(STspkX, STspkY, i, Rx=Rx, Ry=Ry):
    return abs(
        np.mean(STspkX[0, np.where((STspkX[1] > i * 1000 - 50) & (STspkX[1] < i * 1000 + 50))[0]]) - Rx) < .1 and abs(
        np.mean(STspkY[0, np.where((STspkY[1] > i * 1000 - 50) & (STspkY[1] < i * 1000 + 50))[0]]) - Ry) < .1


a = 100
b = 400

clusterA = []
separate = []
import copy

show_plot = 1
for tried in range(20):
    fpath = os.getcwd() + '/Spike/'
    index = [ele for ele in os.listdir(fpath) if ele.rfind('TIM') < 0 if ele.rfind('%02d4' % tried) >= 0]
    print index, '0%d4' % tried
    N = len(index)
    ii = 0
    cluster = []
    for ele in sorted(index, reverse=True):
        ii = 1
        if show_plot:
            plt.figure()
            plt.suptitle(ele)
            plt.subplot(2, 1, ii)

        clist = np.load(fpath + ele).tolist()[1:]
        spk = t.convertposxN(clist)
        spkY = t.convertposxN(clist)
        spkX = t.convertposxN(clist, 'x')
        spk = np.asarray(spk)
        STlist = np.load(fpath + ele[:-4] + '_STIM.npy').tolist()
        STspk = t.convertposxN(STlist)
        STspk = np.asarray(STspk)

        tmp = np.where(abs(spk[0] - .5) < .1)[0]
        y, x = np.histogram(spk[1] * 1. / 1000, bins=np.arange(0, np.max(spk[1]) * 1. / 1000, 0.05))
        #    plt.plot(x[:-1],y,'b')
        #    spk=spk[:,tmp]
        TS = np.max(spk[1])
        plt.plot(STspk[1] * 1. / 1000, STspk[0] * 4096, 'or', zorder=-2)
        plt.plot(spk[1] * 1. / 1000, spk[0] * 4096, 'ok', zorder=-1, ms=2)

        fpath = os.getcwd() + '/Graph/'
        index = [ele for ele in os.listdir(fpath) if ele.rfind('TIM') < 0 if ele.rfind('040') >= 0]

        G = pickle.load(open(fpath + index[0]))
        positions = []
        names = []
        STspkX = t.convertposxN(STlist, 'x')
        STspkY = t.convertposxN(STlist)
        for i in xrange(4):
            tmp = [
                np.mean(STspkX[0, np.where((STspk[1] > (a + i) * 1000 - 500) & (STspk[1] < (a + i) * 1000 + 500))[0]]),
                np.mean(STspkY[0, np.where((STspk[1] > (a + i) * 1000 - 500) & (STspk[1] < (a + i) * 1000 + 500))[0]])]
            print tmp
            names.append(tmp)
            positions.append(tmp[0] * 100 + tmp[1])
        pivot1 = np.arange(4)[np.argsort(positions)]

        pivot = pivot1.argsort()

        offset = .15

        for j in xrange(4):
            for i in np.arange(a + j, b, 4):
                if np.max(y[np.where((x > i) & (x < i + offset))[0]]) > 100 and np.max(
                        y[np.where((x > i - offset) & (x < i))[0]]) < 20:

                    if isRandomST(STspkX, STspkY, i):
                        cluster.append(pivot[j])
                        if show_plot:
                            plt.plot([i, i], [0, 4096], 'y', lw=2)

                            tmp = np.where((STspk[1] > i * 1000 - 500) & (STspk[1] < i * 1000 + 500))[0]
                            plt.plot(STspk[1][tmp] * 1. / 1000, STspk[0][tmp] * 4096, 'or', zorder=-2)
                            tmp = np.where((spk[1] > i * 1000 - 500) & (spk[1] < i * 1000 + 500))[0]
                            plt.plot(spk[1][tmp] * 1. / 1000, spk[0][tmp] * 4096, 'ok', zorder=-1, ms=2)
                    else:
                        cluster.append(pivot[j])

        clusterA.extend(cluster)
        y, x = np.histogram(cluster, bins=np.arange(-1, 5))
        separate.append(y)
        if show_plot:
            ax = plt.gca()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            strength = ele.split('AMPA_test_')[-1].split('_')[0]
            ax.set_ylim([0, 4200])
            plt.ylabel('%d %%' % (float(strength) * 100))

            plt.xlabel('Time [s]')

            ax.yaxis.set_ticks([0, 4096])
            ii = 2
            plt.subplot(2, 1, ii)

            barlist = plt.bar(x[:-1] - .4, y, color='r')
            ax = plt.gca()
            plt.ylabel('# of elicited SBEs')
            barlist[2].set_facecolor('r')
            plt.xlim([-.8, 4 - .4])

            ax.set_xticks(range(4))
            ax.set_xticklabels(['%.02f %.02f' % (names[pivot1[j]][0], names[pivot1[j]][1]) for j in xrange(4)])
plt.figure()
ax = plt.gca()
y, x = np.histogram(clusterA, bins=np.arange(-1, 5))

barlist = plt.bar(x[:-1] - .4, y, color='r')
plt.ylabel('# of elicited SBEs')
barlist[2].set_facecolor('r')
plt.xlim([-.8, 4 - .4])

ax.set_xticks(range(4))
ax.set_xticklabels(['%.02f %.02f' % (names[pivot1[j]][0], names[pivot1[j]][1]) for j in xrange(4)])
# plt.close('all')
# plt.ion()
plt.figure()
separateN = np.asarray(separate) * 1. / ((b - a) / 4)
plt.bar(x[:-1], separateN.mean(0), yerr=separateN.std(0))
for i in xrange(separateN.shape[1]):
    plt.plot(separateN[:, i] * 0 + i - .5, separateN[:, i], 'ok')
plt.ylim([0, .8])
fpath = os.getcwd() + '/Graph/'
index = [ele for ele in os.listdir(fpath) if ele.rfind('TIM') < 0 if ele.rfind('040') >= 0]
plt.figure()
STspk = np.asarray(t.convert2xN(STlist))

XXX = []
for count, i in enumerate(xrange(a, a + 4)):
    g = G.subgraph(np.unique(STspk[0, np.where((STspk[1] > i * 1000 - 500) & (STspk[1] < i * 1000 + 500))[0]]))
    pos = nx.get_node_attributes(g, 'pos').values()
    x = np.array([ele[0] for ele in pos])
    y = np.array([ele[1] for ele in pos])
    d = ((x - np.mean(x)) ** 2 + (y - np.mean(y)) ** 2)
    c = 40 * (1 - (d / (0.08 ** 2)))
    XXX.append([np.mean(x), np.mean(y)])
    plt.text(np.mean(x) + .1, np.mean(y) + .1, '%d-> %.02f,%.02f' % (pivot[count], names[count][0], names[count][1]))
    plt.scatter(x, y, c=c)
plt.colorbar()
plt.xlim([0, 1])
plt.ylim([0, 1])

plt.ion()

from scipy.stats import ttest_rel

test = np.asarray(separate)[:, 1:]
testM = np.zeros((test.shape[1], test.shape[1]))
for i in xrange(len(testM.T)):
    for j in xrange(len(testM.T)):
        testM[i, j] = ttest_rel(test[:, i] - 3, test[:, j])[1]

print testM