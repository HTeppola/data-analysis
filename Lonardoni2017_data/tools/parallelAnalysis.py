import os
import sys
import copy
import tools.utilSpkTrain as UST
from tools.multiPP_ST_ref import runFunctionsInParallel
import pylab as plt
import scipy.cluster.hierarchy as sch
import numpy as np

sys.path.append('/home/alonardoni/8192correlation/')

UST.mfrMin = 0.1
hlist=[]
Tmax=0

def construct_ISs(fNameList, n_work=10, destPath=''):
    task2exec = []
    for fname in fNameList:
        if fname.rfind('.mat') > 0:

            dest = fname[:-4] + '_IS2.npy'
            task2exec.append([workerIS, [fname, dest, destPath], ])
        elif fname.rfind('.npy') > 0 and fname.rfind('map') < 0 and fname.rfind('STIM') < 0 and fname.rfind('IS') < 0:
            folder = fname.split('/Spike/')[0]
            fname = fname.split('/Spike/')[-1]
            path2file = folder + '/Spike/' + fname
            dest = folder + '/' + fname[:-4] + '_IS2.npy'
            if fname[:-4] + '_IS.npy' in os.listdir(folder) and 0:
                print("already analyzed")
                print(dest)

            else:
                task2exec.append([workerIS, [path2file, dest, destPath], ])
    print(task2exec)
    if any(task2exec):
        runFunctionsInParallel(listOf_FuncAndArgLists=task2exec,  maxAtOnce=n_work,SZORDAN=0)


def construct_maps(folder, fname, n_work=20):
    global hlist, Tmax

    if fname.rfind('.mat') > 0:
        hlist = UST.LoadSpikeTrains(folder + fname, ReadHeader=True)
        dest = folder + fname[:-4]
    elif fname.rfind('.npy') > 0 and fname.rfind('map') < 0 and fname.rfind('STIM') < 0 and fname.rfind('IS') < 0:

        path2file = folder + fname
        dest = folder.split('/Spike/')[0] + '/' + fname[:-4]
        hlist = np.load(path2file).tolist()[1:]

        if fname[:-4] + '_mapp.npy' in os.listdir(folder.split('/Spike/')[0]) and 0:
            print("already analyzed")
            print(dest)
            return ['done', 'done']
    else:
        return ['done', 'done']

    lh = len(hlist)
    Tmax = np.max([np.max(cell['spikes']) for cell in hlist if any(cell['spikes'])])
    print("Computing Cross-Correlation of", fname)
    print("%d x %d pairs" % (lh, lh))

    mapp = np.zeros((lh, lh))
    mapt = np.ones((lh, lh)) * -100
    task2exec = []
    for i in range(lh):
        task2exec.append([worker, [i, ], ])

    sRCC = runFunctionsInParallel(listOf_FuncAndArgLists=task2exec, maxAtOnce=n_work)
    for i, vmapp, vmapt in sRCC[1]:
        mapp[i, :] = vmapp
        mapt[i, :] = vmapt
    mapp = mapp + mapp.T
    for i in range(len(mapp)):
        mapp[i, i] = 1
    mapt = mapt - mapt.T

    np.save(open(dest + '_mapp2.npy', 'w'), mapp)
    np.save(open(dest + '_mapt2.npy', 'w'), mapt)

    return mapp, mapt



def worker(i):
    global hlist, Tmax

    #twin = 5
    #tbin = .1
    twin = 20
    tbin = 2
    lh = len(hlist)
    toChangeP = np.zeros(lh)
    toChangeT = np.zeros(lh)
    st1 = hlist[i]['spikes']
    for j in range(i + 1, lh):
        st2 = hlist[j]['spikes']
        if len(st1) > 0 and len(st2) > 0:
            tmp = CrossCorrFastOriginal(np.asarray(st1).copy(), np.asarray(st2).copy(), twin, tbin, Tmax)
        else:
            tmp = np.zeros(10)
        tmp = [np.mean(tmp[max(ic - 1, 0):ic + 2]) for ic in range(len(tmp))]
        NN = len(np.unique(tmp))
        if NN < 15:
            tmp = np.zeros(10)
        else:
            print(j, NN)
        # plt.plot(np.arange(-twin,twin+tbin,tbin),tmp)
        #            plt.savefig('/home/alonardoni/8192correlation/testFC/corr%05d.png'%j)
        #            plt.close()
        pivot = np.argmax(tmp)
        #            pivot=np.asarray([len(tmp)/2,len(tmp)/2+1])
        #            toChangeP[j]=np.sum(tmp[pivot])
        #            if len(tmp[tmp>0]):
        #                if tmp[pivot]>np.mean(tmp[tmp>0])+2*np.std(tmp[tmp>0]):
        #                    toChangeP[j]=tmp[pivot]

        toChangeP[j] = tmp[pivot]
        #            if tmp[pivot]<np.mean(tmp)+3*np.std(tmp):
        #                toChangeP[j]=-tmp[pivot]
        toChangeT[j] = (pivot - len(tmp) / 2) * tbin

    return i, toChangeP, toChangeT


def _convert2xN(CellList):
    tmpid = []
    tmptime = []
    for cell in CellList:
        cell['x'] = np.argmax(np.arange(0, 65) >= cell['x'] * 64)
        cell['y'] = np.argmax(np.arange(0, 65) >= cell['y'] * 64)
        tmpid.extend(np.ones(len(cell['spikes'])) * int(64 * (cell['x'] - 1) + cell['y']))
        tmptime.extend(cell['spikes'])
    tmp = np.array([np.transpose(np.array(tmpid)), np.transpose(np.array(tmptime))])
    return tmp


def workerIS(path2file, dest, destPath):
    import sys
    sys.path.append('/home/alonardoni/FunctionalAnalysis/plugins')
    import plugins.CATClustering as ct
    if path2file.rfind('mat') > 0:
        clist = UST.LoadSpikeTrains(path2file, ReadHeader=True)

        spike2n = UST.FormatSpikeTrains(clist)
        spike2n = np.asarray(spike2n)
    else:
        clist = np.load(path2file).tolist()[1:]
        spike2n = np.asarray(_convert2xN(clist))
    # tmp=np.where(spike2n[1]<20000)[0]
    #    spike2n=spike2n[:,tmp]
    pos = dict()

    for ele in clist:
        pos[int(64 * (ele['x'] - 1) + ele['y'])] = {'x': ele['x'], 'y': ele['y']}

    if len(spike2n):
        #        out=pSBE._detect_SMART(spike2n,pos=pos,Tbin=10,Tin=0,N_ACTIVE=.1,DISTSBE=8,Tdist=15,flagS=1,folder=path2file.split('/')[-1])
        if path2file.rfind('alonardoni') > 0:
            dirName = os.path.dirname(os.path.realpath(path2file)).split('alonardoni/')[-1]
        else:
            dirName = os.path.dirname(os.path.realpath(path2file)).split('data/')[-1]
        if destPath == '':
            destPath = '/home/alonardoni/8192correlation/TESTfigALL/' + dirName.replace('/', '_') + '/'
        else:
            destPath = destPath + '/' + dirName.replace('/', '_') + '/'

        figname = path2file.split('/')[-1]
        figname = os.path.splitext(figname)[0]
        print(figname)
        print(destPath)
        try:
            os.stat(destPath)
        except:
            os.mkdir(destPath)
            os.mkdir(destPath + '/propagations/')

        out = ct.get_ISs(path2file, figname=figname, destPath=destPath + '/')
        #        np.save(open(dest,'w'),out[-1])
        np.save(open(dest, 'w'), out[0])
        np.save(open(dest[:-4] + '_cut.npy', 'w'), out[1])
    print("done")
    return dest, 1

def NumCommon(STS1, STS2):
    """
        count number of common elements
        the STS1 STS2 should be ordered
        credits:
            the algorithm comes from http://www.geeksforgeeks.org/union-and-intersection-of-two-sorted-arrays-2/
            .. with a bug
    """
    L1=len(STS1)
    L2=len(STS2)

    i=0
    j=0
    c=0
    while (i<L1) & (j<L2):
        if STS1[i]<STS2[j]:
            i += 1
        elif STS2[j]<STS1[i]:
            j += 1
        else:
            c += 1
            i += 1
            j += 1  # this was missing in "geeksforgee"
    return c


def CrossCorrFastOriginal(spiketrain1, spiketrain2, TW, dt, Tm):
    nstep = TW / dt
    Nstep = 2 * nstep + 1
    CCvett = np.zeros(Nstep)
    L1 = len(spiketrain1)
    L2 = len(spiketrain2)
    STS1 = np.zeros(L1)
    STS2 = np.zeros(L2)


    sqL1L2 = np.sqrt(L1 * L2)
    Expected = L1 * L2 / Tm * dt
    Term1 = Expected / sqL1L2

    for i in range(L1): STS1[i] = int(spiketrain1[i] / dt)
    for i in range(L2): STS2[i] = int(spiketrain2[i] / dt) - nstep

    for k in range(Nstep):
        STS3 = STS2 + k
        for j in range(L2):
            CCvett[k] = NumCommon(STS1, STS3)
            CCvett[k] = CCvett[k] / sqL1L2 - Term1
    return CCvett[::-1]


def fastburstcorrelation(path2file, tbin=10, thre=0.883):
    """
        Procedure to compute a burst to burst correlation:
        - each syncrhonized event is isolated with its own duration
        - correlation beetween the activity
                of neuron i in burst n
                &
                the same neuron i in burst m
          is performed
        - the maximum of the correlation term is summed up over i=1:# of neurons
        - such value is then normalized

        At the end the matrix of correlation of burst events is reordered with
        a dendrogram algorithm

        Event lasted before 100ms are ignored
        Neuron whose activity ends up before two time bin are ignored

        Input:
        CellList:       List of neurons
        ts:             Ending time of the simulation
        tbin:           Binning time

        Output:
        D:              Reordered matrix of Correlation
        Heatmap plot of D
    """
    unique = False
    burst = np.load(path2file)
    if len(burst) == 0:
        print("empty burst")
        return

    N = len(burst[0])
    print("Nlen: ", N)
    nburst = len(burst)
    print("nburst: ", nburst)
    tinit = np.ones(nburst) * 1e10
    tend = np.zeros(nburst)
    tmpburst = []
    cell = copy.deepcopy(burst[0][0])
    for i in range(4096):
        tmpburst.append(copy.deepcopy(cell))
        cell['ID'] = i
    lburst = []
    for tburst in burst:
        lburst.append(copy.deepcopy(tmpburst))
        for cell in tburst:
            lburst[-1][cell['ID']] = copy.deepcopy(cell)
    burst = lburst
    N = len(burst[0])
    for j in range(nburst):
        for i in range(N):
            if len(burst[j][i]['spikes']) > 0:
                tinit[j] = int(min(tinit[j], min(burst[j][i]['spikes']) / tbin))
                tend[j] = int(max(tend[j], max(burst[j][i]['spikes']) / tbin))
    for j in range(nburst):
        for i in range(N):
            if len(burst[j][i]['spikes']) > 0:
                if unique:
                    burst[j][i]['spikes'] = np.unique(np.floor(burst[j][i]['spikes'] / tbin).astype(int))
                else:
                    burst[j][i]['spikes'] = np.floor(burst[j][i]['spikes'] / tbin).astype(int)
                burst[j][i]['spikes'] = (burst[j][i]['spikes'] - tinit[j])

    validbin = tend - tinit
    D=np.zeros((nburst,nburst))
    if nburst >= 2:

        BC = np.zeros((nburst, nburst))
        for h in range(nburst):

            #            plt.figure()
            for k in range(h, nburst):
                sys.stdout.write("\r Burst %s: %s of %s" % (h, k, nburst))
                sys.stdout.flush()
                # tmpcor=0
                nstep = max(5, np.floor(abs(validbin[h] - validbin[k]) / 2))
                tmpcor2 = np.zeros(nstep * 2 + 1)
                if h == k:
                    BC[h, k] = 1
                else:
                    tot = N + 0.
                    for i in range(N):
                        if len(burst[h][i]['spikes']) <= 2 or len(burst[k][i]['spikes']) <= 2:
                            tot -= 1
                        else:
                            # tmpcor=tmpcor+max(np.correlate(burst[h][i],burst[k][i]))/max(1,sum(burst[h][i]),sum(burst [k][i]))
                            # tmpcor2=tmpcor2+(np.correlate(burst[h][i],burst[k][i]))/max(1,sum(burst[h][i]),sum(burst[k][i]))
                            tmp = CrossCorrFastOriginal(burst[h][i]['spikes'], burst[k][i]['spikes'], nstep, 1, 10 ** 10)
                            tmpcor2 += tmp
                    BC[h, k] = max(tmpcor2) / max(tot, 1)
                    # print "index:",np.argmax(tmpcor2),validbin[h],validbin[k]
                    #                    plt.plot((np.arange(len(tmpcor2))-400)*tbin,tmpcor2/max(tot,1))
                    #                    plt.plot((np.argmax(tmpcor2)-400)*tbin,BC[h,k],'ko',ms=10)
                    BC[k, h] = BC[h, k]
                    # print h,k,BC[h,k]

        # Compute and plot dendrogram.
        fig = plt.figure()
        # axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
        Y = sch.linkage(BC, method='complete')
        Z = sch.dendrogram(Y, no_plot=True)

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
        index = Z['leaves']
        D = BC[index, :]
        tmp = index[::-1]
        D = D[:, tmp]
        im = axmatrix.matshow(D, aspect='auto', origin='lower')
        # Plot colorbar.
        plt.colorbar(im)
        plt.xlabel('SBE ID')
        plt.ylabel('SBE ID')
        # Display and save figure.
        folderpath = os.path.join(*path2file.split('/')[0:-1])
        np.save(
            open('/' + folderpath + '/Analysis/' + (path2file.split('/')[-1]).split('burst')[0] + 'BCmatrix.npy', 'w'),
            BC)
        plt.savefig('/' + folderpath + '/Analysis/' + (path2file.split('/')[-1]).split('burst')[0] + 'BC.png')
        plt.figure()
        index = sch.fcluster(Y, t=thre, criterion='distance')
        for i in range(1, max(index) + 1):
            print("trj", i, "appeared at", tinit[index == i] * tbin / 1000, "s")
            a = tinit[index == i]
            plt.plot(a * tbin / 1000, i * np.ones(len(a)), 'o--', ms=20, lw=5,
                     color=plt.get_cmap('jet')(float(i) / (max(index) + 1)))
            plt.hold(True)
            tmp = np.zeros((sum(index == i), sum(index == i)))
            for j in range(sum(index == i)):
                for k in range(j, sum(index == i)):
                    tmp[j, k] = BC[np.arange(nburst + 1)[index == i][j], np.arange(nburst + 1)[index == i][k]]
            print(tmp)
        plt.ylim([0, max(index) + 1])
        plt.savefig('/' + folderpath + '/Analysis/' + (path2file.split('/')[-1]).split('burst')[0] + 'TBC.png')
        np.save(open(path2file[:-4] + '_index.npy', 'w'), index)
        if 0:
            old = 10
            for j in range(4, nburst):
                ii = 0
                index = sch.fcluster(Y, t=old + ii, criterion='distance')
                while max(index) < j:
                    index = sch.fcluster(Y, t=old + ii, criterion='distance')
                    ii += -0.001
                    # print max(index)
                old = old + ii
                print(max(index), "-cluster optimal", old)
    plt.show()
    return D