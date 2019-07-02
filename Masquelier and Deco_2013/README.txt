This readme and the code were contributed by Timothee Masquelier
timothee.masquelier@alum.mit.edu
Aug 2013

This code was used in:

Masquelier T and Deco G (2013) Network bursting dynamics in excitatory
cortical neuron cultures results from the combination of different
adaptive mechanisms. PLoS ONE

Feel free to use/modify but please cite us if appropriate.

We use the Brian simulator described in:
Goodman D, Brette R (2008) Brian: a simulator for spiking neural
networks in python. Front Neuroinformatics 2:5
and available at:
http://www.briansimulator.org/

This code has been tested with:
   - Brian 1.4.0
   - Python 2.6
   - Mac OS and Linux

Main file: main.py (should be called like that "python -i main.py")
Calls param.py to set the parameters (see comments there), launches
the simulation and plots the results.

The current values in param.py corresponds to the baseline simulation
in the paper (1 min of simulated time).

Many more options are provided (with no guaranty), eg sparse
connectivity, inhomogeneous weights, STDP, adaptive thresholds,
inhibitory neurons etc.

Note: the AMPA and NMDA inputs are computed manually at each time step
using a custom @network_operation (NOT using Brian built in
Connection).  For the simplified case of full connectivity with
homogeneous, non-plastic, weights, these contributions are the same
for all neurons, which enormously simplifies computation.

--------------------
Python files:
--------------------
main.py                        main script
param.py                       contains all the parameters
customrefractoriness.py        Brian file to handle both a refractory
                               period and a user-defined reset
                               function
selfuxconnection.py            handles Short Term Plasticity
plotResult.py                  graphical output

---------------------------
Output files: (in ./data/)
---------------------------

spike.%random seed%.%dump number%.mat Output spikes. Format nx2
                                      matrix, first column is neuron
                                      indexes, second is spike times.

dump.%random seed%.mat                Contains g_a ("adapt"), membrane
                                      potentials ("pot_e"),firing
                                      rates ("rate"), facilitation
                                      ("u") and depression ("x")






