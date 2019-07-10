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

scale = 1
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
noiseweights = 1.0*pl.array([8.0, 1.0]) # Set the noise stimulation weights for each synapse
noiserate = 10 # Rate of stimulation, in Hz
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
        #thiscell.mAMPA = 0
        thiscell.maxcurrent = 3000
    else:     
        thiscell.label = 2 # inhibitory
        thiscell.a = 2
        thiscell.b = 0
        thiscell.C = 200
        thiscell.G_l = 10
        thiscell.tau_w = 30
        thiscell.mNMDA = 0.2
        #thiscell.mAMPA = 0
        thiscell.maxcurrent = 3000
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

#Channels
spx_ch = []

#Spike times
spx = []
for c in range(ncells):
    spx.extend(spikevecs[c])
    number_spikes_in_this_cell = len(spikevecs[c])
    spx_ch.extend([c] * number_spikes_in_this_cell)

spx_no= len(spx)

# Arrays
spx = pl.array(spx)
spx_ch = pl.array(spx_ch)

sort_index = pl.argsort(spx) #sort
spx = spx[sort_index]
spx_ch = spx_ch[sort_index]


# Burst detection parameters
min_el = 40 #10
min_spx_in_NB = 40 #10
min_IBI = 10 # 100 ms

# compute inter spike interval
all_isi = pl.diff(spx)

# Find burst start indices, where ISI > 100ms
indices_start = pl.where(all_isi > min_IBI)[0] + 1
indices_end = pl.where(all_isi > min_IBI)[0] 

# add 0 to the beginning
indices_start = pl.concatenate([[0], indices_start])
# add the burst end to the end
indices_end = pl.concatenate([indices_end, [len(all_isi)]])

burst_count = 0

detected_bursts_start_times = []
detected_bursts_end_times = []
detected_bursts_channel_counts = []
detected_bursts_spike_counts = []

for burst_index in range(len(indices_start)):
    print("checking burst index %d" %burst_index)
    start_time = spx[indices_start[burst_index]]
    end_time = spx[indices_end[burst_index]]
    
    # True or False if this one is in burst
    all_bool_indices_in_burst = (spx >= start_time) * (spx <= end_time) 
    # sum all Trues, to find the total spike count in the burst
    spike_count = pl.sum(all_bool_indices_in_burst)
    
    all_channel_indices_in_the_burst = spx_ch[all_bool_indices_in_burst]
    
    channel_count = len(pl.unique(all_channel_indices_in_the_burst))
    
    if spike_count >= min_spx_in_NB and channel_count >= min_el:
        detected_bursts_start_times.append(start_time)
        detected_bursts_end_times.append(end_time)
        detected_bursts_channel_counts.append(channel_count)
        detected_bursts_spike_counts.append(spike_count)
        burst_count += 1
        
    
pl.subplot(2,1,1)

#binsize= duration/20
pl.hist(spx, bins=500)

pl.ylabel('GFR(Hz)')
spikespercell = [len(pl.array(spikevecs[c])) for c in range(ncells)]
firingrate = sum(spikespercell)/ncells/duration*1000
pl.title('cells=%i syns/cell=%i noise=%s rate=%0.1f Hz' % (ncells,len(connections)/ncells,noiseweights,firingrate),fontsize=12)

## Set a clean upper y-axis limit.
#pl.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


#ax = fig.add_subplot()

pl.subplot(2,1,2)
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
        
# start times (green)
for burst_ind in range(burst_count):
    pl.plot([detected_bursts_start_times[burst_ind], detected_bursts_start_times[burst_ind], ], [0, 100], 'g')
# end times (red)
for burst_ind in range(burst_count):
    pl.plot([detected_bursts_end_times[burst_ind], detected_bursts_end_times[burst_ind], ], [0, 100], 'r')

pl.xlabel('Time (ms)')
pl.ylabel('Cell ID')
pl.xlim(tvec[0], tvec[-1])
pl.show()
sc.toc(); pl.pause(0.1)



print('Done.')
