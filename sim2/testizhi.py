"""
TESTIZHI

This script does a very simple test of Izhikevich neurons to see what's going
wrong.

Version: 2013sep13
"""

from neuron import h, init, run
from izhi import pyramidal
import pylab as pl

duration=1000
connweight = 0.0
noiseweight = 10
whichcell=0
connprob=0.1
whichsyn=3
receptors=['AMPA','NMDA','GABAA','GABAB']

## Create cells
print('Creating...')
cells=[]
spikevecs=[]
spikerecorders=[]
dummy = h.Section()
ncells=1000
for c in range(ncells):
    cells.append(pyramidal(dummy))
    spikevecs.append(h.Vector())
    spikerecorders.append(h.NetCon(cells[c], None))
    spikerecorders[-1].record(spikevecs[-1])

## Connect cells
print('Connecting...')
connections=[]
for c1 in range(ncells):
    for c2 in range(ncells):
        if c1!=c2 and connprob>pl.rand():
            connections.append(h.NetCon(cells[c1], cells[c2])) # Connect them
            connections[-1].weight[whichsyn] = connweight
    

## Add inputs
print('Inputting...')
noiseinputs=[] # Create empty list for storing synapses
noisegens=[] # Create random number generators
noiseconns=[] # Create stimulus connections
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
    noiseconn.weight[1] = noiseweight
    noiseconn.delay=2
    noiseconns.append(noiseconn)

tvec = h.Vector()
vvec = h.Vector()
uvec = h.Vector()
ivec = h.Vector()
tvec.record(h._ref_t)
vvec.record(cells[whichcell]._ref_V)
uvec.record(cells[whichcell]._ref_u)
uvec.record(cells[whichcell]._ref_I)

print('Running...')
init()
run(duration)

print('Plotting...')
pl.plot(pl.array(tvec),pl.array(vvec))
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
pl.title('cells=%i syns/cell=%i weight=%0.1f syn=%s noise=%0.1f rate=%0.1f Hz' % (ncells,len(connections)/ncells,connweight,receptors[whichsyn],noiseweight,firingrate),fontsize=12)


print('Done.')