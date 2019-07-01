__author__ = 'dlonardoni'

#
#   Module to build up a network of neurons of "AdExp.mod" type and
#   run a simulation with noise provided by "pregen.mod" and Short
#   Term Plasticity modeled by "tmgsyn.mod"
#

#   EXTERNAL MODULES
import sys
import os
import pylab as plt
import numpy as np
import numpy.random as rnd
import networkx as nx

from neuron import h
from neuron import gui

#   ADD PATH
sys.path.append('nrnTemplate/')
sys.path.append('nrnMod/')
sys.path.append('tools/')

import utils as ut
import tools.tools as tools
from cellTemplate import AdExpcell

#
#   Definition of parameters of the model
#


#   Population parameters
modParexc = {'tau_w': 280, 'G_l': 10, 'a': 2, 'b': 40, 'C': 281, 'E_l': -70.6, 'V_thre': -50.4, 'Delta_T': 2,
             'V_reset': -70.6, 'type': 'undef', 'iEXT': 50, 'amp': 0}
modParinh = {'tau_w': 280, 'G_l': 10, 'a': 2, 'b': 10, 'C': 200, 'E_l': -70, 'V_thre': -50, 'Delta_T': 2, 'V_reset': -80,
             'type': 'undef', 'iEXT': 75, 'amp': 0}  # it was 75

#   Synaptic parameters
modWeight = {'synexc': 100, 'syninh': 100, 'stp': 500, 'stpinh': 500, 'U': 0.95}
#   Stimolation parameters
stimPar = {'x': 0.2, 'y': 0.8, 'rate': 3000, 'dur': 10, 'amp': 600, 'r': 0.001, 'active': -1, 'start': 60000,
           'end': 120000, 'step': 0}

#   Random input parameters
rndPar = {'parnoise': 1, 'start': 1000, 'end': 10000000, 'fast_invl': 0, 'slow_invl': 40, 'burst_len': 1, 'noise': 1,
          'delay': 4, 'peso': 30}

CellList = []

def go():
    h.run()

def EXCchange_peso(peso):
    for cell in CellList:
     if cell.ID_syn.rfind('EXC')>=0:

        for item in cell.PreList:
            item.weight[0]=peso
def INHchange_peso(peso):
    for cell in CellList:
     if cell.ID_syn.rfind('EXC')<0:
        tmp=-1
        for item in cell.PreList:
            item.weight[0]=tmp*peso

def formatNames(gname, name, _tot, i, j, peso, rpeso):
    gname = "4f_" + gname + "_%s_%s_v%sv%s" % (name, _tot, i, j)
    sname = gname + "_%s_%s_" % (peso, rpeso)
    return sname, gname


def CreateAdExpcells(N, G, percInh=.2, varDt=[False], pos={}, replicate=False, RecAll=2):
    """
	#
	#   Procedure to create a cell of AdExp type coupled with an indipendent generator of random spikes.
	#   Parameters are provided by
	#   - modParexc
	#   - modParinh
	#   - rndPar
	#   - stimPar
	#
	#   Input:
	#   N:          number of cells
	#   varDt:      variable time step
	#   pos:        dictionary containing position of each node in G
	#               empty vector it is accepted
	#
	#   Output:
	#   CellList:   global variable containing the list of all neurons created
	#
	"""

    global CellList
    CellList = []
    coord = {}

    for i in range(N):
        coord['x'] = pos[i][0]
        coord['y'] = pos[i][1]
        #   Append to cell list the new neuron accordingly to its type
        rndPar['seed'] = i

        if replicate:
            if G.node[i]['type'] == 'INH':
                modParCopy=modParinh.copy()
            else:
                modParCopy = modParexc.copy()
            for key,value in G.node[i].items():
                modParCopy[key] = value
            rndPar['parnoise'] = 0
            CellList.append(
                AdExpcell(num=i, modPar=modParCopy, rndParam=rndPar, RecAll=RecAll, coord=coord, varDt=varDt,
                          stimPar=stimPar))

        else:
            if rnd.random() < percInh:
                modParinh['type'] = 'INH'
                CellList.append(
                    AdExpcell(num=i, modPar=modParinh, rndParam=rndPar, RecAll=RecAll, coord=coord, varDt=varDt,
                              stimPar=stimPar))
            else:
                modParexc['type'] = 'EXC'
                CellList.append(
                    AdExpcell(num=i, modPar=modParexc, rndParam=rndPar, RecAll=RecAll, coord=coord, varDt=varDt,
                              stimPar=stimPar))


def ConnectNetwork(G, net_w, mAMPA=1, printFlag=False):
    """
	#
	#   Procedure to connect neurons in CellList. Uses "connect_to()" of CellTemplate module.
	#   The excitatory/inhibitory nature of the connection is determined through the ID variable of the cell.
	#   The connection is made of two steps
	#       First the presynaptic neuron is connected to the STP mechanism with a fixed weight provided by "net_w"
	#       Then the STP mechanism is connected to the postsynaptic neuron via a netcon whose weight is update at each event
	#           involving the edge by STP mechanism
	#
	#       presyn ------> STP mechanism ------> postsyn
	#
	#   Input:
	#   G:              Directed graph G in networkx typical format. It defines who is connected to who.
	#   net_w:          Vector containing the fixed weight of inhibitory and excitatory synapses
	#   magnitudeAMPA:  Efficacy of ampa syn
	#   printFlag:      Boolean, print or no information about who is connected with
	#
	#   Output:
	#   "nothing", this procedure add properties to each neuron in CellList
	#
	"""

    for cell in G.nodes():
        #   retreive outgoing edges for each node
        out_con_from = G.edges([cell])
        for link in out_con_from:
            #   all connections originating from cell
            destLink = link[1]
            CellList[cell].connect_to(dest=CellList[destLink], net_w=net_w, mAMPA=mAMPA)
            if printFlag:
                print("Neuron #", cell, " of nature ", CellList[cell].ID_syn, " connects to: ", destLink, "delay: ", \
                    CellList[cell].PreList[-1].delay)


def SetAttribute(G):
    """
	transfer parameters from AdExp class to graph
	"""
    for n in G:

        G.node[n]['type'] = CellList[n].ID_syn
        G.node[n]['iEXT'] = CellList[n].AdExp.iEXT
        G.node[n]['tau_w'] = CellList[n].AdExp.tau_w
        G.node[n]['G_l'] = CellList[n].AdExp.G_l
        G.node[n]['a'] = CellList[n].AdExp.a
        G.node[n]['C'] = CellList[n].AdExp.C
        G.node[n]['E_l'] = CellList[n].AdExp.E_l
        G.node[n]['V_thre'] = CellList[n].AdExp.V_thre
        G.node[n]['Delta_T'] = CellList[n].AdExp.Delta_T
        G.node[n]['V_reset'] = CellList[n].AdExp.V_reset
        G.node[n]['b'] = CellList[n].AdExp.b
        G.node[n]['num'] = CellList[n].AdExp.num


def NetworkMaker(mAMPA=1, ifNMDA=1, graph=nx.DiGraph(), N=256, peso=110, rpeso=30, ts=60000, r=0.0, tipo=0,
                 sigma=0.0, percInh=0.2, tot=0, group='', replicate=False, stp=500,
                 vardt=[True, 1, 1], stpinh=500, magnitude=0.0, vblock=-55, kb=5, U=.95,recAll=2):
    """
	#
	#   Generate, create, simulate a network of neuron of AdExp type (see AdExp.mod)
	#       coupled with an indipendent generator of random spikes (see pregen.mod).
	#   Synapses are conduptance based integrated with a STP mechanism (see tmgsyn.mod).
	#
	#   Usage:
	#   Provide the number of neurons "N" and specify the type "tipo" of the graph to generate
	#       See "Input" for more details on the parameters to define
	#
	#
	#   Input:
	#   N:          number of neurons in the network
	#   peso:       synaptic weight
	#   rpeso:      strenght of random spike generator
	#   ts:         duration of the simulation
	#   p:          probability [0,1] of an edge to exists in random graph based network
	#   r:          radius (positive) of interaction of each neuron in a radius graph based network
	#   tipo:       specifies which type of graph will be generated. Accepted values
	#                   0 -> provide the graph (networkx format)
	#                   1 -> generate random graph G(N,p)
	#                   2 -> generate radius graph G(N,r)
	#                   3 -> generate a graph where each node has a gaussian N(0,sigma)
	#                           distribution as regard connection distance G(N,sigma)
	#   sigma:      variance of the gaussian for the third type of graph
	#
	#   Output:
	#   List containing all instances for starting the simulation
	#
	#   Debugging info are printed runtime
	#
	"""

    modWeight['U'] = U
    modWeight['synexc'] = peso
    modWeight['syninh'] = peso
    modWeight['stp'] = stp
    modWeight['stpinh'] = stpinh
    modParexc['mNMDA'] = magnitude * ifNMDA
    modParexc['mAMPA'] = 1 - magnitude
    modParinh['mAMPA'] = 1 - magnitude
    modParinh['mNMDA'] = magnitude * ifNMDA
    modParexc['vblock'] = vblock
    modParinh['vblock'] = vblock
    modParexc['kb'] = kb
    modParinh['kb'] = kb
    rndPar['peso'] = rpeso

    #   set tree folder for saving results of the simulation
    path='Results/'
    groupname = group + '/'

    if not os.path.exists(path + groupname):
        os.makedirs(path + groupname)
        os.makedirs(path + groupname + '/fig/')
        os.makedirs(path + groupname + '/Spike/')
        os.makedirs(path + groupname + '/Graph/')
        os.makedirs(path + groupname + '/Burst/')
        os.makedirs(path + groupname + '/Burst/Analysis')
        os.makedirs(path + groupname + '/Burst/Electrode')
        os.makedirs(path + groupname + '/input_res/')
    import time
    start_time = time.time()
    #   generate the graph
    sname=''
    gname=''
    gnameold=''
    pos=dict()
    G=nx.DiGraph()
    if graph == {}:
        if tipo == 0:  # provide the graph in networkx standard
            if sigma == 0:
                gname += 'radius_%s_%s_' % (r, N)
                sname += '4frad_%s_%s_' % (r, peso)
            else:
                gname += 'gauss_%s_%s_' % (sigma, N)
                sname += '4fgauss_%s_%s_' % (sigma, peso)
            G, pos = tools.OpenGraph('/Topologies/' + gname)
            pos = nx.get_node_attributes(G, 'pos')
            gnameold = gname
            gname = '4f' + gname + '%s_' % tot
    else:
        sname += graph.name
        G = graph
        pos = nx.get_node_attributes(graph, 'pos')
        gname = '4f' + sname
        gnameold = gname
    sname += '%s_%s_%s_' % (peso, rpeso, tot)

    print("graph processed in ", time.time() - start_time, " s \n")
    # tools.SaveGraph(G,fname=gname,folder='/res/'+groupname+'/Graph/')
    partial_time = time.time()
    readout = False
    if readout:

        posx = np.array([ele[0] for ele in list(pos.values())])
        posy = np.array([ele[1] for ele in list(pos.values())])
        index = np.where((abs(posx - .5) < .1) & (abs(posy - .5) < .1))[0]
        N = len(G.nodes())
        print(N)
        for j in range(N, N + 10):
            for i in index:
                G.add_path([i, j])

        for i in range(N, N + 10):
            G.node[i]['pos'] = [.5, .5]
        print(len(G.nodes()))
        pos = nx.get_node_attributes(G, 'pos')
    CreateAdExpcells(N=len(G.nodes()), varDt=vardt, pos=pos, replicate=replicate, G=G, percInh=percInh, RecAll=2)
    print("cells created in ", time.time() - partial_time, " s \n")
    if recAll==2:
        for i in range(10):
            CellList[i].RecAll_traces()
    partial_time = time.time()

    #   connect the network according to the generated/provided graph
    ConnectNetwork(G, modWeight, mAMPA=mAMPA, printFlag=False)

    print("cells connected in ", time.time() - partial_time, " s\n")

    h.tstop = ts

    print("Running\n")

    if tipo == 0:
        gtype = 'imported from ../Topologies/' + gnameold
    else:
        gtype = 'new graph created'
    if stimPar['active'] == -1:
        stimuli = {'stim': 'Spontaneous activity'}
    else:
        stimuli = {'stim': 'Evoked activity', 'rate': stimPar['rate'], 'step': stimPar['step'],
                   'start': stimPar['start']}
    HEADER = {'ts': ts / 1000, 'graph': gname, 'gtype': gtype, 'pinh': percInh, 'synexc': modWeight['synexc'],
              'syninh': modWeight['syninh'], 'synnoise': rndPar['peso'], 'stimuli': stimuli}
    save_summary(tipo, gnameold, gname, percInh, path, groupname, ts, N)

    return HEADER


def StartSmallNetTest(peso=37.5, rpeso=18, U=0.5, magnitude=0.05, vblock=-35, kb=6):
    """
	#   Simple test network, connects:
	#
	#        iEXT --> 0 ----> 1 -----> 2
	#
	#   The synapses are all excitatories
	#
	#   Input:
	#   amp:        Current amplitude provided to neuron 0
	#   peso:       Synaptic strength
	#
	#   Output:
	#   Plot of the synaptic current integrated by neuron 1 vs time
	#   Plot of the voltage of the three neurons vs time
	#
	"""

    rndPar['start'] = 1000
    modWeight['U'] = U
    modWeight['synexc'] = peso
    modWeight['syninh'] = peso
    modWeight['stp'] = 200
    modWeight['stpinh'] = 100
    modParexc['vblock'] = vblock
    modParinh['vblock'] = vblock
    modParinh['kb'] = kb
    modParexc['mNMDA'] = magnitude
    modParexc['mAMPA'] = 1 - magnitude
    modParinh['mAMPA'] = 1 - magnitude
    modParinh['mNMDA'] = magnitude
    rndPar['peso'] = rpeso
    modParexc['kb'] = kb

    percInh = 0
    NCells = 3
    stimPar['active'] = -1
    pos={i: np.random.random(2) for i in range(NCells)}
    CreateAdExpcells(N=NCells, G={},  pos=pos, varDt=[True, 1, 1],percInh=percInh, RecAll=2)
    G = nx.DiGraph()
    G.add_path([0, 2])
    G.add_path([0, 1])
    #G.add_path([1, 2])

    ConnectNetwork(G, modWeight, mAMPA=1 - magnitude, printFlag=True)

    '''
	#
	#   Procedure to run a "StartSmallNetTest()" simulation
	#
	'''


    for i in range(len(CellList)):
        CellList[i].AdExp.v0_block = vblock
        CellList[i].AdExp.k_block = kb
        CellList[i].RecAll_traces()
    CellList[2].AdExp.rpeso = 0


    h.tstop = 61000
    h.run()

    plt.figure()
    axV = plt.subplot(NCells, 3, 1)
    axG = plt.subplot(NCells, 3, 2)
    axI = plt.subplot(NCells, 3, 3)
    for i in range(NCells):
        plt.subplot(NCells, 3, 3 * i + 1,sharex=axV)
        ut.PlotVoltage(CellList[i])
        plt.subplot(NCells, 3, 3 * i + 2,sharex=axG)
        ut.PlotG(CellList[i])
        plt.subplot(NCells, 3, 3 * i + 3, sharex=axI)
        ut.PlotCurrent(CellList[i])
    plt.show()
    return CellList


def save_summary(tipo, gnameold, gname, pinh, path, groupname, ts, N):
    if tipo == 0:
        gtype = 'imported from ../Topologies/' + gnameold
    else:
        gtype = 'new graph created'
    if stimPar['active'] == -1:
        stimuli = {'stim': 'Spontaneous activity'}
    else:
        stimuli = {'stim': 'Evoked activity', 'rate': stimPar['rate'], 'step': stimPar['step'],
                   'start': stimPar['start']}
    HEADER = {'ts': ts / 1000, 'graph': gname, 'gtype': gtype, 'pinh': pinh, 'synexc': modWeight['synexc'],
              'syninh': modWeight['syninh'], 'synnoise': rndPar['peso'], 'stimuli': stimuli}

    f1 = open(path + groupname + 'summary.txt', 'w')
    print('Simulation of %s neurons over %s (s)' % (N, ts / 1000), file=f1)
    print('Graph: %s' % gname, end=' ', file=f1)
    if tipo == 0:
        print(' imported from /Topologies', file=f1)
    else:
        print(' created', file=f1)
    print('Percentage of inh neurons: %s' % pinh, file=f1)
    print('Exc syn strength: %s' % modWeight['synexc'], file=f1)
    print('Inh syn strength: %s' % modWeight['syninh'], file=f1)
    print('Noise strength: %s' % rndPar['peso'], file=f1)
    if stimPar['active'] == -1:
        print('Spontaneous activity', file=f1)
    else:
        print('Evoked activity', file=f1)
        print('\t Time interval among stimula: %s (ms)' % stimPar['rate'], file=f1)
        print('\t # of stimula: %s' % stimPar['step'], file=f1)
        print('\t First delivered stimulus stimula: %s (ms)' % stimPar['start'], file=f1)
    f1.close()
    return HEADER


# noinspection PyTypeChecker
def save_voltage(path2save):
    trace2saveALL = []
    for cell in CellList:
        if cell.voltage:
            trace2save = dict()
            trace2save['time'] = np.array(cell.record['timeT'])
            trace2save['voltage'] = np.array(cell.record['voltage'])
            trace2saveALL.append(trace2save)
    np.save(open(path2save, 'w'), trace2saveALL)
    return trace2saveALL


# noinspection PyTypeChecker
def save_traces(path2save):
    trace2saveALL = []
    for cell in CellList:
        if cell.traces:
            trace2save = dict()
            trace2save['gexc'] = np.array(cell.record['gexc'])
            trace2save['ginh'] = np.array(cell.record['ginh'])
            trace2save['g'] = np.array(cell.record['g'])
            trace2save['time'] = np.array(cell.record['timeT'])
            trace2save['voltage'] = np.array(cell.record['voltage'])
            trace2save['Isyn'] = np.array(cell.record['Isyn'])
            trace2save['w'] = np.array(cell.record['w'])
            trace2save['all_c'] = np.array(cell.record['all_c'])
            trace2save['Inmda'] = np.array(cell.record['Inmda'])
            trace2saveALL.append(trace2save)
    np.save(open(path2save, 'w'), trace2saveALL)
    return trace2saveALL


def replicate_sim(G, degree='degree', inv=0):
    # TORCINI FUNCTION
    import copy
    pivot = np.array([i for i in G.nodes() if G.node[i]['type'].rfind('X') > 0])
    H = G.subgraph(pivot)
    pos = nx.get_node_attributes(H, 'pos')
    GR = nx.DiGraph()
    deg=dict()
    ideg=dict()
    exec('deg= np.array(G.%s().values())' % degree)
    exec('ideg= np.array(G.%s().keys())' % degree)
    deg = deg[pivot]
    ideg = ideg[pivot]
    sdeg = deg.argsort()  # [::-1]
    if inv:
        sdeg = sdeg[::-1]
    R = list()
    R.append([H.node[i]['R'] for i in H.nodes()])
    R.append([int(i) for i in H.nodes()])
    R = np.array(R)
    sR = R[0, :].argsort()

    for i in range(len(H.nodes())):
        relabel = sdeg[i]
        tmp = nx.DiGraph()
        tmp = nx.compose(tmp, copy.deepcopy(H.subgraph(R[1, sR[i]])))

        mapping = {R[1, sR[i]]: ideg[relabel]}
        tmp = nx.relabel_nodes(tmp, mapping)
        tmp.node[ideg[relabel]]['pos'] = pos[ideg[relabel]]
        GR = nx.compose(GR, tmp.subgraph(ideg[relabel]))

    for i in G.nodes():
        if G.node[i]['type'].rfind('H') > 0:
            GR = nx.compose(GR, G.subgraph(i))
    GR.add_edges_from(G.edges())

    R = np.array([GR.node[i]['R'] * 1000 for i in GR.nodes() if i in pivot])
    exec('deg=[value for key,value in GR.%s().iteritems() if key in pivot]' % degree)
    plt.scatter(deg, R, color='b')
    return GR

