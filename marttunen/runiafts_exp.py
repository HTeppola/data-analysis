# Python code for running NEST simulator with given network structure using Tsodyks' synapses.
#
# The script loads simulation attributes from a mat-file whose name is given on command line. The
# simulation output is saved to a gcf (spike train) and dat (membrane potential) file: The names
# of these files are saved to the mat-file as additional variables by the script.
#
# Input (read from a file given as first argument)
#     Ne - Number of excitatory cells
#     Ni - Number of inhibitory cells
#     M - Connectivity matrix such that M[i][j] denotes the existence of connection FROM i TO j.
#         The size of the matrix is assumed (Ne + Ni) x (Ne + Ni)
#     Ie - The static input drive given to each excitatory cell
#     Ie - The static input drive given to each inhibitory cell
#     wee - The synaptic weights E->E
#     wei - The synaptic weights I->E
#     wie - The synaptic weights E->I
#     wii - The synaptic weights I->I
#     dt - Simulation time step (ms)
#     t_sim - Simulation time (ms)
#     syndelay - the mean synaptic delay
#     syndelaystd - the STD of synaptic delay (for each synapse individually)
#     noisestd - The STD of the white noise applied to the membrane potential of each cell
#     seed - A random number seed used to produce different sample paths
#
# Output:
#     spikefile - Name of the gcf file containing the spike trains of the neurons. The first column
#                 of the data denotes the indices of the neurons that spiked (from N+1 to N+Ne
#                 excitatory indices, N+Ne+1 to 2N inhibitory indices), as the second column denotes
#                 the corresponding spike times.
#     vmfile - Name of the dat file containing the membrane potential of the first excitatory neuron
#              of the network.
#
# Tuomo Maki-Marttunen
# Last modified 9.1.2013

import nest, scipy.io, sys, numpy
import nest.voltage_trace, nest.raster_plot

mat = scipy.io.loadmat(sys.argv[1]+".mat", struct_as_record=True)
Ne = mat['Ne'].item();
Ni = mat['Ni'].item();
M = mat['M'].tolist();
Ie = float(mat['Ie']);
Ii = float(mat['Ii']);
wee = float(mat['wee'])
wei = float(mat['wei'])
wie = float(mat['wie'])
wii = float(mat['wii'])
dt          = float(mat['dt'])
t_sim       = float(mat['t_sim'])
syndelay    = float(mat['syndelay'])
syndelaystd = float(mat['syndelaystd'])
noisestd    = float(mat['noisestd'])
seed = int(mat['seed'])
print "seed = %i" % seed

N = Ne + Ni

nest.SetKernelStatus({"resolution": dt, "print_time": False,"overwrite_files":True, "grng_seed": seed})
neuron_pars = {"C_m"      : 30.0,
               "tau_m"    : 30.0,
               "t_ref"    : 3.0,
               "E_L"      : 0.0,
               "V_th"     : 15.0,
               "V_reset"  : 13.5,
               "V_m"      : 13.5,
               "I_e"      : Ie,
               "tau_syn_ex" : 3.0,
               "tau_syn_in" : 3.0
              }
nest.CopyModel('iaf_psc_exp', 'iaf_e', neuron_pars)
nest.CopyModel('iaf_psc_exp', 'iaf_i', neuron_pars)
nest.SetDefaults('iaf_i', {"t_ref": 2.0, "I_e": Ii}) #faster recovery for inhibitory neurons

noise_e = nest.Create("noise_generator", Ne)
if Ni > 0:
    noise_i = nest.Create("noise_generator", Ni)
else:
    noise_i = 0

for n in range(0,Ne):
    nest.SetStatus([noise_e[n]], {"mean": 0.0, "std": noisestd})
for n in range(0,Ni):
    nest.SetStatus([noise_i[n]], {"mean": 0.0, "std": noisestd})

synapse_pars_e = {"tau_psc": 3.0,
                  "tau_rec": 800.0,
                  "tau_fac": 0.0,
                  "U": 0.5,
                  "delay": 0.2,
                  "weight": wee,
                  "u": 0.5,
                  "x": 1.0,
                  "y": 0.0
                 }
synapse_pars_i = {"tau_psc": 3.0,
                  "tau_rec": 100.0,
                  "tau_fac": 1000.0,
                  "U": 0.04,
                  "delay": 0.2,
                  "weight": wie,
                  "u": 0.04,
                  "x": 1.0,
                  "y": 0.0
                 }
nest.CopyModel('tsodyks_synapse', 'syn_ee', synapse_pars_e)
nest.CopyModel('tsodyks_synapse', 'syn_ei', synapse_pars_e)
nest.SetDefaults('syn_ei', {"weight": wei})
nest.CopyModel('tsodyks_synapse', 'syn_ie', synapse_pars_i)
nest.CopyModel('tsodyks_synapse', 'syn_ii', synapse_pars_i)
nest.SetDefaults('syn_ii', {"weight": wii})
nest.SetDefaults("voltmeter",     {"flush_after_simulate": True, "withtime":True, "to_file":True, "to_file":True})
nest.SetDefaults("spike_detector",{"flush_after_simulate": True, "withgid":True, "withtime":True, "to_file":True})

neuron_e = nest.Create("iaf_e", Ne)
if Ni > 0:
    neuron_i = nest.Create("iaf_i", Ni)
else:
    neuron_i = 0

spiketrain = nest.Create("spike_detector")
vm = nest.Create("voltmeter")

for n in xrange(0,Ne):
    nest.Connect([neuron_e[n]],spiketrain)
for n in xrange(0,Ni):
    nest.Connect([neuron_i[n]],spiketrain)
nest.Connect(vm,[neuron_e[0]])
nest.SetStatus(vm,[{"to_file": True, "withtime": True}])

nest.Connect(noise_e, neuron_e)
if Ni > 0:
    nest.Connect(noise_i, neuron_i)

for ni in xrange(0,Ne):
    for nj in xrange(0,Ne):
        if M[ni][nj]:
            syndel = syndelay + syndelaystd*numpy.random.randn(1).item()
            syndel = round(syndel/dt)*dt
            nest.Connect([neuron_e[ni]],[neuron_e[nj]],{'weight': wee, 'delay': syndel}, model='syn_ee' )
    for nj in xrange(0,Ni):
        if M[ni][Ne+nj]:
            syndel = syndelay + syndelaystd*numpy.random.randn(1).item()
            syndel = round(syndel/dt)*dt
            nest.Connect([neuron_e[ni]],[neuron_i[nj]],{'weight': wie, 'delay': syndel}, model='syn_ie' )
for ni in xrange(0,Ni):
    for nj in xrange(0,Ne):
        if M[Ne+ni][nj]:
            syndel = syndelay + syndelaystd*numpy.random.randn(1).item()
            syndel = round(syndel/dt)*dt
            nest.Connect([neuron_i[ni]],[neuron_e[nj]],{'weight': wei, 'delay': syndel}, model='syn_ei' )
    for nj in xrange(0,Ni):
        if M[Ne+ni][Ne+nj]:
            syndel = syndelay + syndelaystd*numpy.random.randn(1).item()
            syndel = round(syndel/dt)*dt
            nest.Connect([neuron_i[ni]],[neuron_i[nj]],{'weight': wii, 'delay': syndel}, model='syn_ii' )

nest.SetKernelStatus({"rng_seeds": range(seed+1, seed+2)})
nest.Simulate(t_sim)
spikefile = nest.GetStatus(spiketrain,'filenames')  #in some version "filename", in other "filenames"
vmfile = nest.GetStatus(vm,'filenames')
mat['spikefile'] = spikefile[0][0]
mat['vmfile'] = vmfile[0][0]
scipy.io.savemat(sys.argv[1]+".mat",mat,oned_as='row')






