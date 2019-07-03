"""
TESTIZHI

This script does a very simple test of Izhikevich neurons to see what's going
wrong.

Version: 2013sep13
"""

from neuron import h, init, run
from adexp import createcell
import pylab as pl
import sciris as sc

sc.tic()

duration = 1000 # Set the duration 
connweights = [10.0, 10.0]
noiseweights = [1.0, 10.0]
whichcell = 0
connprob = 0.0
whichsyns = [0,1]
receptors = ['w0','w1','w2']

## Create cells
print('Creating...')
cells = []
spikevecs = []
spikerecorders = []
dummy = h.Section()
ncells = 100
for c in range(ncells):
    thiscell = createcell(dummy, c)
    cells.append(thiscell)
    spikevecs.append(h.Vector())
    spikerecorders.append(h.NetCon(cells[c], None))
    spikerecorders[-1].record(spikevecs[-1])
sc.toc()

## Connect cells
print('Connecting...')
connections=[]
for c1 in range(ncells):
    for c2 in range(ncells):
        if c1!=c2 and connprob>pl.rand():
            connections.append(h.NetCon(cells[c1], cells[c2])) # Connect them
            for syn in whichsyns:
                connections[-1].weight[syn] = connweights[syn]
sc.toc()

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
    noiseinput.interval = 10 # Take inverse and then convert from Hz to ms^-1
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

tvec = h.Vector()
vvec = h.Vector()
wvec = h.Vector()
ivec = h.Vector()
tvec.record(h._ref_t)
vvec.record(cells[whichcell]._ref_vv)
wvec.record(cells[whichcell]._ref_ww)
ivec.record(cells[whichcell]._ref_gEXC)
sc.toc()

print('Running...')
init()
run(duration)
sc.toc()

print('Plotting...')
pl.figure()
pl.subplot(3,1,1)
pl.plot(pl.array(tvec), pl.array(vvec))
pl.subplot(3,1,2)
pl.plot(pl.array(tvec), pl.array(wvec))
pl.subplot(3,1,3)
pl.plot(pl.array(tvec), pl.array(ivec))

pl.figure()
pl.scatter(pl.array(spikevecs[whichcell]), pl.zeros(len(spikevecs[whichcell])))
for c in range(ncells): 
    ex = pl.array(spikevecs[c])
    if len(ex)>0:
        why = 100+c*pl.ones(len(spikevecs[c]))
        pl.scatter(ex, why)
        pl.show()
    else:
        print('No spikes for cell %i' % c)
pl.xlabel('Time (ms)')
pl.ylabel('Voltage & cell ID')
pl.xlim(tvec[0],tvec[-1])
firingrate = float(sum(len(pl.array(spikevecs[c])) for c in range(ncells)))/ncells/duration*1000
pl.title('cells=%i syns/cell=%i weight=%s noise=%s rate=%0.1f Hz' % (ncells,len(connections)/ncells,connweights,noiseweights,firingrate),fontsize=12)
sc.toc()


print('Done.')