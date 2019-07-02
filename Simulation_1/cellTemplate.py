__author__ = 'dlonardoni'
import random as rnd
import numpy as np

from neuron import h


#
#   Template of a neuron of AdExp class
#

class AdExpcell:
    def __init__(self, num, modPar={}, rndParam={}, RecAll=2, coord={}, varDt=[False,1,1], stimPar={'active': 0}):
        """
                Adaptive Exponential Integrate and Fire neuron
                see AdExp.mod for details

                The parameter of the model are given by the vector "modPar"
                each value is perturbated with gaussian noise
        """
        # CHECK WHAT APPEN IF THE RUMOR IS something+N(0,0.025)
        # INSTEAD OF something*N(1,0.025)
        #   define which type of neuron it is
        self.PL = []
        self.Pcon = []
        self.P2con = []
        self.play = []


        self.mNMDA = modPar['mNMDA']
        rndforce = 0.0025

        nmda = True
        self.nmda = nmda
        self.RecAll = RecAll
        self.soma = h.Section(cell=self)
        self.voltage = False
        if modPar['type'].rfind('EXC') >= 0:
            self.ID_syn = 'EXC'
            self.AdExp = h.AdExp(0.5, sec=self.soma)
            self.AdExp.rpeso = rndParam['peso']
            self.AdExp.mNMDA = modPar['mNMDA']
            self.AdExp.v0_block = modPar['vblock']
            self.AdExp.k_block = modPar['kb']
                # print "EXC"
        elif modPar['type'].rfind('INH') >= 0:

            self.ID_syn = 'INH'
            self.AdExp = h.AdExp(0.5, sec=self.soma)
            self.AdExp.v0_block = -5
            self.AdExp.k_block = modPar['kb']
            self.AdExp.rpeso = rndParam['peso']
            self.AdExp.mNMDA = modPar['mNMDA']
        else:
            self.ID_syn = 'UNDEF'
            print("Neither excitatory nor inhibitory neuron. The neural type if undefined and no synapse will be formed!")
        # create new section & add AdExp mechanism
        self.AdExp.mAMPA = modPar['mAMPA']
        self.varDt = varDt
        self.ID = "AdExp"
        self.num = num
        self.soma.nseg = 1
        self.soma.diam = 10
        self.soma.cm = 1

        #   define parameters of the model
        self.AdExp.iEXT = modPar['iEXT']
        self.AdExp.tau_w = modPar['tau_w'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.G_l = modPar['G_l'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.a = modPar['a'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.C = modPar['C'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.E_l = modPar['E_l'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.V_thre = modPar['V_thre']  # *(1+rnd.gauss(0,rndforce)*rndParam['parnoise'])
        self.AdExp.Delta_T = modPar['Delta_T'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.V_reset = modPar['V_reset'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])
        self.AdExp.b = modPar['b'] * (1 + rnd.gauss(0, rndforce) * rndParam['parnoise'])

        if stimPar['step'] == 1:
            self.AdExp.num = num / 50
        else:
            self.AdExp.num = int(num)
        # define position of the cell
        self.pos = coord.copy()
        #   PostList contains all netcon from STP to destination neuron
        self.PostList = []
        #   StpList contains all the STPs of the incoming connections
        self.StpList = []
        #   PreList contains all netcon from neuron to STP mechanism
        self.PreList = []

        self.record = dict()
        self.record['stimulated'] = []

        tota = 0
        if stimPar['active'] != -1:
            LP = len(stimPar['positions'])
            for i in range(LP):
                x = stimPar['positions'][i][0]  # 0.8# stimPar[['active']%4
                y = stimPar['positions'][i][1]  # 0.8#(stimPar[['active']-x)/4
                d = ((self.pos['x'] - x) ** 2 + (self.pos['y'] - y) ** 2)
                if d < .08 ** 2:
                    self.record['stimulated'].append(num)
                    self.Stim = h.SpikeGenerator(0.5, sec=self.soma)
                    self.Stim.start = 1000 * 3 + 1000 * i
                    self.Stim.end = 10 ** 10
                    self.Stim.fast_invl = rndParam['fast_invl']
                    self.Stim.slow_invl = 1000 * LP
                    self.Stim.burst_len = rndParam['burst_len']
                    self.Stim.delay = rndParam['delay']
                    self.Stim.noise = 0.001
                    self.Stim.seed(1.)
                    self.StimInput = h.NetCon(self.Stim, self.AdExp, -20, 0.1, 100, sec=self.soma)
                    self.StimInput.weight[0] = 45 * (1 - (d / (0.08 ** 2)))
                    print("stimulation at node %s!!! S:%.2f" % (num, self.StimInput.weight[0]))
                    self.StimInput.weight[1] = 1
                    self.record['Stimspk'] = h.Vector()
                    self.nc_STIMspike = h.NetCon(self.Stim, None, -20, 0, 1, sec=self.soma)
                    self.nc_STIMspike.record(self.record['Stimspk'])
                    tota += 1

        self.playS = dict()
        self.playS['test'] = h.Section()
        #   initialize and inject some noise
        self.RandomStim = h.SpikeGenerator(0.5, sec=self.soma)
        if len(rndParam) > 0:
            self.RandomStim.start = rndParam['start']
            self.RandomStim.end = rndParam['end']
            self.RandomStim.fast_invl = rndParam['fast_invl']
            self.RandomStim.slow_invl = rndParam['slow_invl']
            self.RandomStim.burst_len = rndParam['burst_len']
            self.RandomStim.delay = rndParam['delay']
            self.RandomStim.noise = rndParam['noise']
            self.RandomStim.seed(rndParam['seed'])

        # direct coupling (not via STP mechanism)


        self.RandomInput = h.NetCon(self.RandomStim, self.AdExp, -20, 0.1, 2000, sec=self.soma)

        # record spike time stamps
        self.nc_spike = h.NetCon(self.AdExp, None, sec=self.soma)

        self.record['spk'] = h.Vector()
        self.nc_spike.record(self.record['spk'])
        #        self.record['rspk']=h.Vector()
        #        self.nc_rspike=h.NetCon(self.RandomStim,None,-20,0,1,sec=self.soma)
        #        self.nc_rspike.record(self.record['rspk'])
        if self.RecAll == 2:
            self.record['iAMPA']=h.Vector()
            self.record['iAMPA'].record(self.AdExp._ref_iAMPA, sec=self.soma)
            self.record['iNMDA']=h.Vector()
            self.record['iNMDA'].record(self.AdExp._ref_iNMDA, sec=self.soma)
            self.record['iGABA']=h.Vector()
            self.record['iGABA'].record(self.AdExp._ref_iGABA, sec=self.soma)

        #        self.nc_spike.record(self.record['spk'])
        self.traces = False
        if varDt[0]:
            # setup the "LOCAL TIME STEP" integrator
            self.Hines = h.CVode()
            self.Hines.active(varDt[1])
            self.Hines.use_local_dt(varDt[2])

    def add_spk(self, spk, mAMPA):
        self.play2 = h.VecStim()
        self.Pcon2 = h.NetCon(self.play2, self.AdExp, sec=self.soma)
        peso = 1
        self.playS['test'] = h.Section()
        for ele in spk:
            self.play.append(h.VecStim())
            self.play[-1].play(h.Vector(ele))

        self.RandomStim.start = 10 ** 10
        aSpk = [ele for item in spk for ele in item]


        np.random.seed(int(np.random.random() * 1000))
        self.play2.play(h.Vector(sorted(np.random.random(40 * np.max(aSpk) // 1000) * np.max(aSpk))))
        #        self.Pcon2.weight[1]=.11
        self.Pcon2.weight[0] = 20000
        no_stp = 0
        if no_stp:
            self.Pcon = []
            for conn in self.play:
                self.Pcon.append(h.NetCon(conn, self.AdExp, sec=self.soma))
                self.Pcon[-1].weight[1] = mAMPA
                self.Pcon[-1].weight[0] = peso
        else:
            for _ in range(len(spk)):
                self.PL.append(h.stp_2(self.soma(0.5)))
                self.PL[-1].tau_rec = 100  # net_w['stp']
                self.PL[-1].U = .2  # net_w['U']
                self.PL[-1].z0 = 0
                self.PL[-1].y0 = .2

            for i, conn in enumerate(self.play):
                self.Pcon.append(h.NetCon(conn, self.PL[i], sec=self.soma))
                self.Pcon[-1].weight[1] = mAMPA
                self.Pcon[-1].weight[0] = peso

                self.P2con.append(h.NetCon(self.PL[i], self.AdExp, sec=self.soma))
                self.P2con[-1].delay = 0.5
                self.P2con[-1].weight[1] = mAMPA
                h.setpointer(self.P2con[-1]._ref_weight[0], 'new', self.PL[i])

    def set_param(self, modPar):

        self.AdExp.tau_w = modPar['tau_w'].value
        self.AdExp.a = modPar['a'].value
        self.AdExp.C = modPar['C'].value
        self.AdExp.E_l = modPar['E_l'].value
        self.AdExp.V_thre = modPar['V_thre'].value  #
        self.AdExp.Delta_T = modPar['Delta_T'].value
        self.AdExp.b = modPar['b'].value

        for conn in self.Pcon:
            conn.weight[0] = modPar['gsyn'].value

    def connect_to(self, dest, net_w, mAMPA):  # ,delay=0,weight=1):
        """
            Procedure to connect cell
            Each connection is build up as follows:
            - A new STP mechanism is added to the destination neuron
            - A netcon whose weight is provided by net_w is extablished between
                the starting neuron "self" and the STP mechanism
            - A netcon whose weight is provided by STP mechanism is extablished
                between the STP mechanism and the destination neuron

            The delay is computed as a function of the distance among neurons

            This method take in account excitatory or inhibitory synapses

            Input:
            self:       Starting neuron
            dest:       Destination neuron
            net_w:      Synapse weight
            N:          Number of cells
        """
        #   define new STP mechanism
        dest.StpList.append(h.stp(dest.soma(0.5)))
        dest.StpList[-1].tau_rec = net_w['stp']
        dest.StpList[-1].U = net_w['U']
        dest.StpList[-1].z0 = np.random.random() * .4 + .4
        dest.StpList[-1].y0 = np.random.random() * .2

        if self.ID_syn.rfind('INH') >= 0:
            dest.StpList[-1].tau_rec = net_w['stpinh']
        # create netcon self ----> STP
        self.PreList.append(h.NetCon(self.AdExp, dest.StpList[-1], sec=self.soma))
        self.PreList[-1].weight[1] = mAMPA
        if self.ID_syn.rfind('EXC') >= 0:
            self.PreList[-1].weight[0] = net_w['synexc']
            #   create netcon STP ---> dest
            self.PostList.append(h.NetCon(dest.StpList[-1], dest.AdExp, sec=self.soma))
            self.PostList[-1].delay = 0

            # set a pointer to the adjusted weight "new" modified by STP mechanism
            h.setpointer(self.PostList[-1]._ref_weight[0], 'new', dest.StpList[-1])
            self.PostList[-1].weight[1] = mAMPA
            self.PostList[-1].weight[2] = self.num
            if self.pos['x'] < 1:
                self.PreList[-1].delay = 0.5 + 3.14 * abs(
                    (((self.pos['x'] - dest.pos['x']) ** 2 + (self.pos['y'] - dest.pos['y']) ** 2) ** 0.5))

        elif self.ID_syn.rfind('INH') >= 0:
            self.PreList[-1].weight[0] = -(net_w['syninh'])
            self.PostList.append(h.NetCon(dest.StpList[-1], dest.AdExp, sec=self.soma))

            h.setpointer(self.PostList[-1]._ref_weight[0], 'new', dest.StpList[-1])
            self.PostList[-1].weight[1] = mAMPA
            self.PostList[-1].weight[2] = mAMPA
            if self.pos['x'] < 1:
                self.PreList[-1].delay = 0.5 + 3.14 * abs(
                    (((self.pos['x'] - dest.pos['x']) ** 2 + (self.pos['y'] - dest.pos['y']) ** 2) ** 0.5))

        else:
            print("Since ", self.ID_syn, " is neither an excitatory nor an inhibitory neuron not connection is made!")

    def destroy(self):
        del self.PreList
        del self.PostList
        del self.nc_spike
        del self.StpList
        del self.record
        del self.AdExp
        del self.Hines

    def destroy_fit(self):
        del self.Pcon
        del self.Pcon2
        try:
            del self.PL
        except:
            pass
        del self.nc_spike
        del self.record
        del self.RandomInput
        del self.RandomStim

        del self.AdExp
        del self.Hines

    def IR_RecAll(self):
        self.record['time'] = h.Vector()
        self.record['voltage'] = h.Vector()
        self.Hines.record(self.AdExp._ref_vv, self.record['voltage'], self.record['time'], sec=self.soma)
        self.voltage = True

    def NoRec(self):
        del self.record['time']
        del self.record['voltage']
        self.voltage = False

    def RecAll_traces(self):
        print('I am recording everything')
        self.record['gAMPA'] = h.Vector()
        self.record['gGABA'] = h.Vector()
        self.record['gNMDA'] = h.Vector()
        self.record['time'] = h.Vector(np.arange(1000, 60000))
        self.record['voltage'] = h.Vector()
        self.record['iAMPA'] = h.Vector()
        self.record['iGABA'] = h.Vector()
        self.record['iNMDA'] = h.Vector()
        self.record['w'] = h.Vector()
        self.record['iTotal'] = h.Vector()


        self.traces = True

        if self.varDt[0]:
            self.Hines.record(self.AdExp._ref_vv, self.record['voltage'], self.record['time'], sec=self.soma)
            self.Hines.record(self.AdExp._ref_iAMPA, self.record['iAMPA'], self.record['time'], sec=self.soma)
            self.Hines.record(self.AdExp._ref_iGABA, self.record['iGABA'], self.record['time'], sec=self.soma)

            self.Hines.record(self.AdExp._ref_gAMPA, self.record['gAMPA'], self.record['time'], sec=self.soma)
            self.Hines.record(self.AdExp._ref_gGABA, self.record['gGABA'], self.record['time'], sec=self.soma)
            self.Hines.record(self.AdExp._ref_ww, self.record['w'], self.record['time'], sec=self.soma)
            self.Hines.record(self.AdExp._ref_iTotal, self.record['iTotal'], self.record['time'], sec=self.soma)
            if self.nmda:
                self.Hines.record(self.AdExp._ref_iNMDA, self.record['iNMDA'], self.record['time'], sec=self.soma)
                self.Hines.record(self.AdExp._ref_gNMDA, self.record['gNMDA'], self.record['time'], sec=self.soma)

        else:
            self.record['time'].record(h._ref_t)
            self.record['voltage'].record(self.AdExp._ref_vv)
            self.record['w'].record(self.AdExp._ref_w)
            self.record['iAMPA'].record(self.AdExp._ref_iAMPA)
            self.record['iGABA'].record(self.AdExp._ref_iGABA)
            self.record['gAMPA'].record(self.AdExp._ref_gAMPA)
            self.record['gGABA'].record(self.AdExp._ref_gGABA)
            self.record['iTotal'].record(self.AdExp._ref_iTotal)
            if self.nmda:
                self.record['iNMDA'].record(self.AdExp._ref_iNMDA)
                self.record['gNMDA'].record(self.AdExp._ref_gNMDA)