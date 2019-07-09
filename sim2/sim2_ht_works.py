"""
SIM2

A simple implementation of an AdExp neuronal network.

Version: 2019jul08
"""

from neuron import h, init, run
from adexp import createcell
import pylab as pl
import sciris as sc

# Without specifying the integrator, it dies on every spike
integrator = h.CVode()
integrator.active(1)
integrator.use_local_dt(1)
pl.seed(205738) # Reproducible results (hopefully)

sc.tic()

scale =1
duration = 10000 # Set the duration in ms
n_e = 80*scale # Number of excitatory cells
n_i = 20*scale # Number of inhibitory cells
ncells = n_e + n_i # Number of cells

globalconnweight = 1.2 # Global modulation for all weights (both connectivity and noise)
connweights = {'e->e': [3.0, 5.0], # Set the connecivity weights for each synapse type, 2nd is for AMPA
               'e->i': [3.0, 5.0], 
               'i->e': [-85.0, 0], 
               'i->i': [-85.0, 0]} 
#noiseweights = [8.0, 1.0] # Set the noise stimulation weights for each synapse
noiseweights = 0.5*pl.array([8.0, 1.0]) # Set the noise stimulation weights for each synapse
noiserate = 100 # Rate of stimulation, in Hz
connprob = 0.2 # The connection probability
conndelay = 40 # Average connection delay (ms)
whichcell = 0 # The cell to record example traces from (shouldn't matter)
whichsyns = [0,1] # Which synapses/receptors to stimulate
receptors = ['w0','w1','w2'] # The names of the receptors (see AdExp.mod)

## Create cells
print('Creating...')
cells = []
spikevecs = []
spikerecorders = []
dummy = h.Section()

for c in range(ncells):
    thiscell = createcell(dummy)
    if c<n_e: 
        thiscell.label = 1 # excitatory -- arbitrary convention
        thiscell.a = 2
        thiscell.G_l = 12
        thiscell.tau_w = 300
        thiscell.mNMDA = 0.2
        thiscell.mAMPA = 0
        thiscell.maxcurrent = 6000
    else:     
        thiscell.label = 2 # inhibitory
        thiscell.a = 2
        thiscell.b = 0
        thiscell.C = 200
        thiscell.G_l = 10
        thiscell.tau_w = 30
        thiscell.mNMDA = 0.2
        thiscell.mAMPA = 0
        thiscell.maxcurrent = 6000
    cells.append(thiscell)
    spikevecs.append(h.Vector())
    spikerecorders.append(h.NetCon(cells[c], None))
    spikerecorders[-1].record(spikevecs[-1])
sc.toc(); pl.pause(0.1)

## Connect cells
print('Connecting...')
connections=[]
for c1 in range(ncells):
    for c2 in range(ncells):
        if c1!=c2 and connprob>pl.rand():
            connections.append(h.NetCon(cells[c1], cells[c2])) # Connect them
            connections[-1].delay = conndelay*(0.5+pl.rand())
            if cells[c1].label == 1 and cells[c2].label == 1: weightkey = 'e->e'
            if cells[c1].label == 1 and cells[c2].label == 2: weightkey = 'e->i'
            if cells[c1].label == 2 and cells[c2].label == 1: weightkey = 'i->e'
            if cells[c1].label == 2 and cells[c2].label == 2: weightkey = 'i->i'
            for syn in whichsyns:
                connections[-1].weight[syn] = globalconnweight*connweights[weightkey][syn]*(0.5+pl.rand())
sc.toc(); pl.pause(0.1)

## Add inputs
print('Inputting...')
noiseinputs = [] # Create empty list for storing synapses
noisegens = [] # Create random number generators
noiseconns = [] # Create stimulus connections
for c in range(ncells): 
    noisegen = h.Random()
    noisegen.MCellRan4(c,c*2)
    noisegen.negexp(1)
    noisegens.append(noisegen)
    
    noiseinput = h.NetStim()
    noiseinput.interval = 1.0/noiserate*1000 # Take inverse and then convert from Hz to ms^-1
    noiseinput.number = 1e9
    noiseinput.start = 0
    noiseinput.noise = 1
    noiseinput.noiseFromRandom(noisegen)
    noiseinputs.append(noiseinput)
    
    noiseconn = h.NetCon(noiseinput, cells[c])
    for syn in whichsyns:
        noiseconn.weight[syn] = noiseweights[syn]
    noiseconn.delay=2
    noiseconns.append(noiseconn)
sc.toc(); pl.pause(0.1)

print('Setting up recording...')
tvec = h.Vector()
vvec = h.Vector()
wvec = h.Vector()
#ivec = h.Vector()
#iivec = h.Vector()
ga = h.Vector()
gn = h.Vector()
gg = h.Vector()


tvec.record(h._ref_t)
vvec.record(cells[whichcell]._ref_vv)
wvec.record(cells[whichcell]._ref_ww)
#ivec.record(cells[whichcell]._ref_gEXC)
#iivec.record(cells[whichcell].ref_gINH)
ga.record(cells[whichcell]._ref_gAMPA)
gn.record(cells[whichcell]._ref_gNMDA)
gg.record(cells[whichcell]._ref_gGABA)

sc.toc(); pl.pause(0.1)

print('Running...')
init()
run(duration)
sc.toc(); pl.pause(0.1)

print('Plotting...')
pl.figure()
pl.subplot(6,1,1)
pl.plot(pl.array(tvec), pl.array(vvec))
pl.xlabel('Time(ms)')
pl.ylabel('V(mV)')

pl.subplot(6,1,2)
pl.plot(pl.array(tvec), pl.array(wvec))
pl.xlabel('Time(ms)')
pl.ylabel('I(pA)')

#pl.subplot(7,1,3)
#pl.plot(pl.array(tvec), pl.array(ivec))
#pl.xlabel('Time(ms)')
#pl.ylabel('gExc(nS)')

#pl.subplot(7,1,4)
#pl.plot(pl.array(tvec), pl.array(iivec))
#pl.xlabel('Time(ms)')
#pl.ylabel('gINH(nS)')

pl.subplot(6,1,4)
pl.plot(pl.array(tvec), pl.array(ga))
pl.xlabel('Time(ms)')
pl.ylabel('gAMPA(nS)')

pl.subplot(6,1,5)
pl.plot(pl.array(tvec), pl.array(gn))
pl.xlabel('Time(ms)')
pl.ylabel('gNMDA(nS)')

pl.subplot(6,1,6)
pl.plot(pl.array(tvec), pl.array(gg))
pl.xlabel('Time(ms)')
pl.ylabel('gGABA(nS)')


fig = pl.figure()
#pl.subplot(2,1,1)
#pl.plot(pl.array(tvec), pl.array())
#pl.xlabel('Time(ms)')
#pl.ylabel('GFR(Hz)')
#ax = fig.add_subplot()
for c in range(ncells): 
    ex = pl.array(spikevecs[c])
    if len(ex)>0:
        why = c*pl.ones(len(spikevecs[c]))
        if c<n_e: 
            spikecolor = 'red' # excitatory -- arbitrary convention
        else:
            spikecolor = 'blue' # inhibitory
        pl.scatter(ex, why, c=spikecolor, alpha=1.0)
    else:
        print('No spikes for cell %i' % c)
pl.xlabel('Time (ms)')
pl.ylabel('Cell ID')
pl.xlim(tvec[0], tvec[-1])
spikespercell = [len(pl.array(spikevecs[c])) for c in range(ncells)]
firingrate = sum(spikespercell)/ncells/duration*1000
pl.title('cells=%i syns/cell=%i noise=%s rate=%0.1f Hz' % (ncells,len(connections)/ncells,noiseweights,firingrate),fontsize=12)
pl.show()
sc.toc(); pl.pause(0.1)


print('Done.')
