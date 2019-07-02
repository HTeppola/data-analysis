# === Define parameters ========================================================
set_global_preferences(useweave=True) # to compile C code

def iround(value):
    return int(value+0.5)

# simulation-related parameters   
#randomSeed = 3
simTime     = 60*second # 5*60*60*second     # total simulation time [ms]  
defaultclock.dt = 0.025*ms # simulation step length [ms]
monitorTime = simTime#4*second
splitTime = 10*second # periodic log / dump
onlyRecord60 = True
dump_a = False # dump adaptation variable

condDelay = 3*ms# experimental estimation 1.5-1.9md (Nakanishi 1998).

# synaptic weights
factor = 1.0 # mean activity increases with factor
w_ee = 8.75*factor # 1.6 - 1.7
w_ei = w_ee
w_ii = 1.0 
w_ie = 1.0

# sparseness
factorS = 1.0 
sparseness_ee = factorS*1.0 
sparseness_ei = factorS*1.0 
sparseness_ii = factorS*1.0
sparseness_ie = factorS*1.0 

matrixStructure = 'dense'# 'dense' or 'sparse'
binaryWeight = True 
removeSelfConnection = False
useLogNormalWeightDist = False

sigma = 6.5*mV # gaussian white noise (Langevin equation). Note external 2.4kHz Poisson AMPA imput with Brunel and Wang's parameters corresponds to ~ 1.4 mV
delta = 0.08  # increase to get oscillations (eg 0.12), see Buehlmann & Deco 2008
              # Albantakis and Deco 2011 used 0.08
             
compensation = 10.0 # Buehlmann, Albantakis etc use 10

#_______________________________________________________________
# Spike Timing-Dependent Plasticity
useSTDP = False
if useSTDP:
    dA_pre = 1e-3 # with 1e-3, about 2000s needed for convergence
    dA_post = - 0.55*dA_pre
    tau_post=34.0*ms #(source: Bi & Poo 2001)
    tau_pre=17.0*ms #(source: Bi & Poo 2001) 
    wmax = 2*w_ee
#_________________________________________
# Adaptation 
useAdaptation = True
if useAdaptation:
    alpha_Ca = 0.0145 
    tau_Ca = 4000*ms 
    g_AHP = 10*nS
    E_K = -80*mV
#_______________________________________________________________
# Adaptive threshold (Kobayashi et al 2009)
useAdaptiveThr = False
if useAdaptiveThr:
    alpha_1 = 0*mvolt
    alpha_2 = 5*mvolt # Tim Feb 2012: ~3mV seems to lead to 0.5s long NS
    tau_1 = 200.0*ms
    tau_2 = 4000.0*ms    
#_______________________________________________________________
# Short Term Plasticity (seems necessary to get network spikes)
useSTP = True
if useSTP: # use tau = 0*ms to unplug facilitation / depression
    tau_f = 1600*ms # Marina uses 2000ms
    tau_d = 800*ms # shoud be > 200ms (Tim, Sept 2011)
    U = 0.025 # Tim Feb 2012: ~0.8. If too big, no trough between NS peak and residual activity
              # Tim Feb 2012: ~0.8. If too small, huge NS peak

#_______________________________________________________________
# Kick (to trigger a NS)
useKick = False
if useKick:
    kickTime = 10*second
    kickStrength = .1
       
# === Parameters determining model dynamics =====

# pyramidal cells 
Cm_e        = 0.5*nF     # [nF] total capacitance 
gl_e        = 25.0*nS    # [ns] total leak conductance 
El_e        = -70.0*mV   # [mV] leak reversal potential 
Vth_e       = -50.0*mV   # [mV] threshold potential
Vr_e   = -55.0*mV   # [mV] reset potential 
tr_e   = 2.0*ms   # [ms] refractory time 
 
# interneuron cells 
Cm_i        = 0.2*nF   # [nF] total capacitance 
gl_i        = 20.0*nS   # [ns] total leak conductance 
El_i        = -70.0*mV    # [mV] leak reversal potential 
Vth_i       = -50.0*mV    # [mV] threshold potential 
Vr_i   = -55.0*mV    # [mV] reset potential 
tr_i   = 1.0*ms     # [ms] refractory time

# network structure
N = 1000
#scale = 0.5 # w.r.t Wang 2002
N_e   = iround(.8*N)
N_i   = 0*iround(.2*N) # use 0 if there is no inhibitory neurons

# AMPA receptor (AMPAR) 
E_ampa   = 0.0*mV     # [mV] synaptic reversial potential 
t_ampa   = 2.0*ms     # [ms] exponential decay time constant  
g_ampa_e = 104.0*nS*(1+compensation*delta) / N # ~ x3.5 to compensate NMDA blocker
g_ampa_i = 81.0*nS*(1+compensation*delta) / N # ~ x3.5 to compensate NMDA blocker

g_ext_e   = 2.08*nS     # [nS] maximum conductance from external to pyramidal cells 
g_ext_i   = 1.62*nS    # [nS] maximum conductance from external to interneuron cells
wext_e=g_ext_e/g_ampa_e
wext_i=g_ext_i/g_ampa_i
 
# GABA receptor (GABAR) 
E_gaba   = -70.0*mV      # [mV] synaptic reversial potential 
t_gaba   = 10.0*ms   # [ms] exponential decay time constant
g_gaba_e = 1250*nS / N 
g_gaba_i = 973*nS / N 
 
# NMDA receptor (NMDAR) 
E_nmda   = 0.0*mV     # [mV] synaptic reversial potential 
ts_nmda   = 100.0*ms   # [ms] decay time of NMDA currents 
tx_nmda   = 2.0*ms   # [ms] controls the rise time of NMDAR channels 
alfa_nmda   = 0.5*kHz      # [kHz] controls the saturation properties of NMDAR channels
g_nmda_e = 327*nS*(1-delta) / N # ~ x0.1 if adding APV
g_nmda_i = 258*nS*(1-delta) / N # ~ x0.1 if adding APV

def stimulus(t):
    return 2400*Hz

def weakStimulus(t):
    if(t>2000*ms):
        stim=0*Hz
    else:
        stim=200*Hz
    return stim

def shortStimulus(t):
    if(t>=0*second):
        stim=0*Hz
    else:
        stim=1800*Hz
    return stim

def timedLog(txt):
    t=time.localtime()
    print '%02d'  % t[3] +  ':%02d' % t[4] + ':%02d' % t[5] + ' - ' + txt
    






