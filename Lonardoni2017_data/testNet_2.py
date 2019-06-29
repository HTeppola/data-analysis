import pylab as plt
import cellCultureNet as Net
import networkx as nx
from tools import tools
from tools import utils as ut
import pickle
import numpy as np


G=ut.crop_graph(pickle.load(open("savedGraph/4fgauss_0.005_1024_1_.graph")),.5)

pos=nx.get_node_attributes(G,'pos')
G=G.subgraph([key for key,value in pos.iteritems() if value[0]<1.3 and value[1]<1.3])
G=nx.relabel_nodes(G,{idx:i for i,idx in enumerate(G.nodes())})
pAMPA = 2
ifNMDA = 1
magnitude = 0.05
mAMPA = 1.00
degree = 'control'
inv = 0

_peso = 24 * 1.5
_rpeso = 16
_tot = 1
groupname = "facil_nmda_05_tau_U"
PESO=60
RPESO=16
U=.5
PINH=.2
MODIFIER=0
for index in xrange(1):
    sigma=.005
    ts=1000 * 300


    G.name = "NEW_graph-all_AMPA_test_%s_%s_%d_%s" % (mAMPA, sigma, len(G), inv)
    #    Net.modWeight['U']=.2
    HEADER =Net.NetworkMaker(recAll=0,mAMPA=mAMPA,ifNMDA=ifNMDA,replicate=False,graph=G,stp=150,stpinh=100,N=len(G.nodes()),peso=_peso,rpeso=RPESO,ts=1000*101,tipo=0,sigma=0.0005,percInh=0.2,tot=_tot,r=0,group=groupname,magnitude=magnitude,vblock=-41,kb=6,U=U)
    b = 40.
    Net.INHchange_peso(PESO * 2*(1+MODIFIER))

    # Net.change_ad('b',b)
    i = 4
    j = 3

    peso = PESO  # 37.5
    rpeso = _rpeso * (.6 + .1 * j)

    Net.EXCchange_peso(peso*(1-MODIFIER))
    sname, gname = Net.formatNames(G.name, 'tau_w', _tot, 0, 0, _peso, _rpeso)
    gname += "_%d_NMDA_%0.2d_NOBIC" % (ifNMDA, magnitude)
    sname += "sfinal_%d_NMDA_%0.2f_NOBIC" % (ifNMDA, magnitude)
    Net.SetAttribute(G=G)
    tools.SaveGraph(G, fname=gname, folder="Results/" + groupname + "/Graph/")
    Net.go()
    #ut.PlotRaster(Net.CellList, ts=ts)
    tools.SaveSpikeTrains(Net.CellList, HEADER, fname=sname, folder="Results/" + groupname + "/Spike/")
    spikelist = []
    for cell in Net.CellList:
        try:
            formatData = {'type': cell.ID_syn, "spikes": np.array(cell.record["Stimspk"]), "ID": np.array([cell.num]),
                          "x": np.array([cell.pos["x"]]), "y": np.array([cell.pos["y"]])}
            spikelist.append(formatData)
        except:
            pass
    np.save(open('Results/' + groupname + '/Spike/' + sname + '_STIM.npy', 'w'), spikelist)

    clist = np.load('Results/' + groupname + '/Spike/' + sname + '.npy')

    for cell in clist[1:]:
        plt.plot(cell['spikes'], np.ones_like(cell['spikes']) * cell['y'], 'ok')
    plt.ion()
    tmp = tools.convert2xN(clist[1:])
    for idx in np.unique(tmp[0]):
        spkCount=len(clist[idx + 1]['spikes'])
        print idx, clist[idx + 1]['type'], spkCount, spkCount*1000./Net.h.t
    plt.figure()
    NCells=5
    axV = plt.subplot(NCells, 3, 1)
    axG = plt.subplot(NCells, 3, 2)
    axI = plt.subplot(NCells, 3, 3)
    for i in xrange(NCells):
        plt.subplot(NCells, 3, 3 * i + 1,sharex=axV)
        ut.PlotVoltage(Net.CellList[i])
        plt.subplot(NCells, 3, 3 * i + 2,sharex=axG)
        ut.PlotG(Net.CellList[i])
        plt.subplot(NCells, 3, 3 * i + 3, sharex=axI)
        ut.PlotCurrent(Net.CellList[i])
    plt.show()
    plt.show()
