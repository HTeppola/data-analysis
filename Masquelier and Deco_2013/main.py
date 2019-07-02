"""
"""
from scipy.io import *
from scipy.sparse import *
from scipy import *
from scipy import weave
import os
from numpy import *
import time 

debugg   = False      # set True for debugging (checks the units) 
if not debugg: 
  import brian_no_units 
  
from brian import *
from selfuxconnection import SelfUXConnection

import sys
import getopt
#letters = 'r:'
#keywords = ['randomSeed=']
try:
    opts, args = getopt.getopt(sys.argv[1:], 'r:w:W:a:u:n:t:d:c:s:D:f:', ['randomSeed=', 'weightFile=','recurrentWeight=','adaptation=','utilization=','number=','adaptThreshold=','taud_d=','condDelay=','sigma=','delta=','taud_f='])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
#    usage()
    sys.exit(2)
#opts, extraparams = getopt.getopt(sys.argv[1:],'r:w:',['randomSeed='])
# starts at the second element of argv since the first one is the script name
# extraparms are extra arguments passed after all option/keywords are assigned
# opts is a list containing the pair "option"/"value"

execfile('param.py') # set parameters

# eventually override them:
randomSeed=0
weightFile = None
for o,p in opts:
  if o in ['-r','--randomSeed']:
     randomSeed = int(p)
  elif o in ['-w','--weightFile']:
     weightFile = p
  elif o in ['-W','--recurrentWeight']:
     w_ee = double(p)
  elif o in ['-a','--adaptation']:
     alpha_Ca = double(p)
  elif o in ['-u','--utilization']:
     U = double(p)
  elif o in ['-n','--number']:
     N = int(p)
     N_e   = iround(.8*N)
     N_i   = 0*iround(.2*N)
     g_ampa_e = 104.0*nS*(1+compensation*delta) / N
     g_ampa_i = 81.0*nS*(1+compensation*delta) / N 
     g_gaba_e = 1250*nS / N
     g_gaba_i = 973*nS / N 
     g_nmda_e = 327*nS*(1-delta) / N 
     g_nmda_i = 258*nS*(1-delta) / N  
  elif o in ['-t','--adaptThreshold']:
     alpha_2 = double(p)*ms
  elif o in ['-d','--tau_d']:
     tau_d = double(p)*ms
  elif o in ['-c','--condDelay']:
     condDelay = double(p)*ms
  elif o in ['-s','--sigma']:
     sigma = double(p)*mV
  elif o in ['-D','--delta']:
     delta = double(p)
     g_ampa_e = 104.0*nS*(1+compensation*delta) / N # ~ x3.5 to compensate NMDA blocker
     g_ampa_i = 81.0*nS*(1+compensation*delta) / N # ~ x3.5 to compensate NMDA blocker
     g_nmda_e = 327*nS*(1-delta) / N # ~ x0.1 if adding APV
     g_nmda_i = 258*nS*(1-delta) / N # ~ x0.1 if adding APV
  elif o in ['-f','--tau_f']:
     tau_f = double(p)*ms


print '\nRandom seed = %d' % randomSeed
print 'w_ee = %f' % w_ee
if useAdaptation:
    print 'alpha_Ca = %f' % alpha_Ca
if useSTP:
    print 'U = %f' % U
    print 'tau_d = %f' % tau_d
    print 'tau_f = %f' % tau_f    
if useAdaptiveThr:
    print 'alpha_2 = %f' % alpha_2
print 'N = %d' % N
print 'condDelay = %f' % condDelay
print 'sigma = %f' % sigma
print 'delta = %f' % delta

seed(randomSeed)


################################# 

# common parts
eqs_e = ''
if N_i>0:
    eqs_e += 'ds_gaba/dt = -s_gaba/(t_gaba) : 1\n'
eqs_e += ''' 
ds_ampa/dt = -s_ampa/(t_ampa) : 1
ds_nmda/dt = -s_nmda/(ts_nmda)+alfa_nmda*x*(1-s_nmda) : 1 
dx/dt = -x/(tx_nmda) : 1
s_nmda_in : 1
s_ampa_in : 1
tmp : 1
dv/dt = ''' 
eqs_i = ''' 
ds_gaba/dt = -s_gaba/(t_gaba) : 1
s_nmda_in : 1
s_ampa_in : 1
tmp: 1
dv/dt = '''
#reset_e = 'v=Vr_e'
#reset_i = 'v=Vr_i'

if sigma>0: # noise terms
    eqs_e += 'sigma*xi*(2*gl_e/Cm_e)**.5'
    eqs_i += 'sigma*xi*(2*gl_i/Cm_i)**.5'
 
# common parts 
eqs_e += ' + (-gl_e*(v-El_e)-g_ampa_e*(v-E_ampa)*s_ampa_in-s_nmda_in*g_nmda_e*(v-E_nmda)/(1+exp(-0.062*v/(1*mV))/3.57)'
if N_i>0:
    eqs_e += '-g_gaba_e*(v-E_gaba)*s_gaba'
eqs_i += ' + (-gl_i*(v-El_i)-g_gaba_i*(v-E_gaba)*s_gaba-g_ampa_i*(v-E_ampa)*s_ampa_in-s_nmda_in*g_nmda_i*(v-E_nmda)/(1+exp(-0.062*v/(1*mV))/3.57)'

if useAdaptation: # adapation terms
    eqs_e += '-g_AHP*Ca*(v-E_K))/Cm_e : mvolt \ndCa/dt = - Ca / tau_Ca : 1 \n'
    eqs_i += '-g_AHP*Ca*(v-E_K))/Cm_i : mvolt \ndCa/dt = - Ca / tau_Ca : 1 \n'
#    reset_e += ';Ca+=alpha_Ca'
#    reset_i += ';Ca+=alpha_Ca'
else:
    eqs_e += ')/Cm_e : mvolt \n'
    eqs_i += ')/Cm_i : mvolt \n'

if useSTP: # STP equations
    if tau_d > 0: # depression is at work
        eqs_e += 'dx_d/dt = (1-x_d)/tau_d : 1\n'
    if tau_f > 0: # facilitation is at work
        eqs_e += 'du_f/dt = (U-u_f)/tau_f : 1 \n'
#    eqs_e += 'dlogx_d/dt = (exp(-logx_d)-1)/tau_d : 1\n'
    # TO DO: facilitation with delays
#    if tau_f < Inf:
#        eqs_e += 'dlogu_f/dt = (U*exp(-logu_f)-1)/tau_f : 1 \n'
#    reset_e += ';x_d*=(1-u_f);u_f+=U*(1-u_f)'
if useAdaptiveThr:
    eqs_e += 'dtheta_1/dt = -theta_1/tau_1 : mvolt \ndtheta_2/dt = -theta_2/tau_2 : mvolt \ndtimer/dt = 1.0 : ms'
    eqs_i += 'dtheta_1/dt = -theta_1/tau_1 : mvolt \ndtheta_2/dt = -theta_2/tau_2 : mvolt \ndtimer/dt = 1.0 : ms'
    

#print eqs_e
#print eqs_i
    
# Setting up the populations 
print '\nSetting up the populations ...' 

# dynamically define reset functions
if N_i > 0:
    populations = ['e', 'i']
else:
    populations = ['e']
for pop in populations:
    reset_str = 'def reset_' + pop + '(P,spikes):';
    if useAdaptiveThr:
        reset_str = reset_str + '\n    P.theta_1_[spikes] += alpha_1'
        reset_str = reset_str + '\n    P.theta_2_[spikes] += alpha_2'
        reset_str = reset_str + '\n    P.timer_[spikes] = 0'
        reset_str = reset_str + '\n    P.v_[spikes] = Vr_' + pop
    else:
        reset_str = reset_str + '\n    P.v_[spikes] = Vr_' + pop
    if useAdaptation:
        reset_str = reset_str + '\n    P.Ca_[spikes] += alpha_Ca'        
    exec(reset_str)
        
if useAdaptiveThr:
    scr_e = SimpleCustomRefractoriness(reset_e,period=tr_e,state='v')
    Pe = NeuronGroup(N_e, eqs_e, threshold='(v>theta_1+theta_2+Vth_e)*(timer>tr_e)>0', reset=scr_e, order=1, compile=True, freeze=True)
#        Pe = NeuronGroup(N_e, eqs_e, threshold="(v>theta_1+theta_2+Vth_e)*(timer>tr_e)>0", reset=reset_e, order=1, compile=True, freeze=True)
    Pe.theta_1 = 0
    Pe.theta_2 = 0
else:
    scr_e = SimpleCustomRefractoriness(reset_e,period=tr_e,state='v')
    Pe = NeuronGroup(N_e, eqs_e, threshold=Vth_e, reset=scr_e, order=1, compile=True, freeze=True)
    
Pe.v = El_e # Vr_e
Pe.s_ampa = 0 
if N_i>0:
    Pe.s_gaba = 0
Pe.s_nmda = 0
Pe.x = 0

 
if N_i>0:
    #Pi = NeuronGroup(N_i, eqs_i, threshold=Vth_i, reset=reset_i, refractory=tr_i, order=1, compile=True, freeze=True) 
    if useAdaptiveThr:
        Pi = NeuronGroup(N_i, eqs_i, threshold="(v>theta_1+theta_2+Vth_i)*(timer>tr_i)>0", reset=reset_i, order=1, compile=True, freeze=True)
        Pi.theta_1 = 0
        Pi.theta_2 = 0
    else:
        scr_i = SimpleCustomRefractoriness(reset_i,period=tr_i,state='v')
        Pi = NeuronGroup(N_i, eqs_i, threshold=Vth_i, reset=scr_i, order=1, compile=True, freeze=True)
    #Pi.v = El_i + 1.0*rand(N_i)*(Vth_i-El_i)#El_i 
    Pi.v = El_i # Vr_i
    #Pi.s_ampa = 0 
    Pi.s_gaba = 0

if useAdaptation:
    Pe.Ca = 0.0
    if N_i>0:
        Pi.Ca = 0.0
if useSTP:
    if tau_f>0:
        Pe.u_f = U
    if tau_d>0:
        Pe.x_d = 1
#    Pe.logx_d = 0

## 'Starter' (external stimulation)
#PG=PoissonGroup(N_e+N_i, shortStimulus)
#PGe = PG.subgroup(N_e)
#PGi = PG.subgroup(N_i) 
##PGe = PoissonGroup(N_e,shortStimulus)
##PGi = PoissonGroup(N_i,stimulus) 


# Creating connections 
print 'Creating static connections ...'
 
self_nmda=IdentityConnection(Pe, Pe, 'x', weight=1.0, delay=condDelay)
self_ampa=IdentityConnection(Pe, Pe, 's_ampa', weight=1.0, delay=condDelay)
#self_logx=IdentityConnection(Pe, Pe, 'logx_d', weight=log(1-U), delay=condDelay)
if useSTP:
    self_ux=SelfUXConnection(Pe, tau_f>0, tau_d >0, delay=condDelay, U=U)


#Cpe = IdentityConnection(PGe, Pe, 's_ampa', weight=wext_e)
#Cpi = IdentityConnection(PGi, Pi, 's_ampa', weight=wext_i) 

#from scipy import sparse
#Cei = Connection(Pe, Pi, 's_ampa',delay=condDelay,structure=matrixStructure)
#tmp = 2*sparse.rand(N_e, N_i, density=sparseness_ei, format='csr', dtype=None)
#if binaryWeight:
#    tmp = csr_matrix(tmp.todense()>0)
#Cei.connect(Pe,Pi,tmp*w_ei/sparseness_ei)


# Returns random nR x nC connection matrix. Average non-zero weights is 1
def getRandomConnectionMatrix(sparseness,binaryWeight,nR,nC):
    mat = 1.0*(rand(nR,nC)<=sparseness)
    if not binaryWeight:
        mat *= 2.0*rand(nR,nC)
    return csr_matrix(mat) # switching to sparse matrix does not change anything Connection is dense, but is more efficient if Connection is sparse

Cee = Connection(Pe, Pe, 'tmp',delay=condDelay,structure=matrixStructure) # in fact this connection is fake; AMPA and NMDA contributions are computed manually at each time step. This is useful essentially for STDP.
if N_i>0:
    Cie = Connection(Pi, Pe, 's_gaba',delay=condDelay,structure=matrixStructure)
    Cei = Connection(Pe, Pi, 'tmp',delay=condDelay,structure=matrixStructure)
    Cii = Connection(Pi, Pi, 's_gaba',delay=condDelay,structure=matrixStructure)

if weightFile!=None: #loadWeightMatrix and os.path.exists(os.path.join('..','data','Cee.STDP.mat')):
    print 'Loading weight matrix ' + weightFile
#    tmp=loadmat(os.path.join('..','data','Cee.STDP.mat'))
    tmp=loadmat(weightFile)
    Cee.connect(Pe,Pe,csr_matrix(w_ee*tmp['Cee'].todense()))
    if N_i>0:
        Cei.connect(Pe,Pi,csr_matrix(w_ei*tmp['Cei'].todense()))
        Cii.connect(Pi,Pi,csr_matrix(w_ii*tmp['Cii'].todense()))
        Cie.connect(Pi,Pe,csr_matrix(w_ie*tmp['Cie'].todense()))
else:
    tmp  = getRandomConnectionMatrix(min(1,sparseness_ee+removeSelfConnection*1.0/N_e),binaryWeight,N_e,N_e)
    if useLogNormalWeightDist: # set to True to get log normal weight distribution
        from random import lognormvariate
        for i in arange(N_e):
            for j in arange(N_e):
                if tmp[i,j]>0:
                    tmp[i,j] = lognormvariate(0,.6)
    Cee.connect(Pe,Pe,w_ee/sparseness_ee*tmp)
    if N_i>0:
        Cie.connect(Pi,Pe,w_ie/sparseness_ie*getRandomConnectionMatrix(sparseness_ie,binaryWeight,N_i,N_e))
        Cei.connect(Pe,Pi,w_ei/sparseness_ei*getRandomConnectionMatrix(sparseness_ei,binaryWeight,N_e,N_i)) 
        Cii.connect(Pi,Pi,w_ii/sparseness_ii*getRandomConnectionMatrix(min(1,sparseness_ii+removeSelfConnection*1.0/N_i),binaryWeight,N_i,N_i))

##Cii = Connection(Pi, Pi, 's_gaba', weight=N_i/(N_i-1), delay=condDelay, structure='dense') # note: a full i-i connectivity seems to prevent network spikes...
##Cie = Connection(Pi, Pe, 's_gaba', weight=.8*w_inh, delay=condDelay, structure='dense') # note: NS still possible with full i-e connectivity
##Cei = Connection(Pe, Pi, 's_ampa', weight=correctionSTP, delay=condDelay,structure='dense') # note: NS still possible with full i-e connectivity
##Cee = Connection(Pe, Pe, 's_ampa', weight=correctionSTP, delay=condDelay,structure='dense') # note: NS still possible with full i-e connectivity

if removeSelfConnection:    #Remove self-connections
    if N_i>0:
        for i in xrange(len(Pi)):
            Cii[i,i]=0
    for i in xrange(len(Pe)):
        Cee[i,i]=0


Cee.compress() 
if N_i>0:
    Cei.compress() 
    Cie.compress() 
    Cii.compress() 


if useSTDP:
    stdp = ExponentialSTDP(Cee, tau_pre, tau_post, dA_pre, dA_post, wmin=0, wmax=wmax, update='additive',interactions='all')

### Short term plasticity
if useSTP:
#    stp_ee = STP(Cee, taud=tau_d, tauf=tau_f, U=U)
#    stp_ei = STP(Cei, taud=tau_d, tauf=tau_f, U=U)
#    _ux = asarray(stp_ee._contained_objects[0].ux)
#    _ux = asarray(stp_ee.vars._var_array['ux'])
    if tau_f>0:
        _u = asarray(Pe.u_f)
    if tau_d>0:
        _x = asarray(Pe.x_d)
#    _logx = asarray(Pe.logx_d)


_Pe_s_nmda_in = asarray(Pe.s_nmda_in)
_Pe_s_nmda = asarray(Pe.s_nmda)

_Pe_s_ampa_in = asarray(Pe.s_ampa_in)
_Pe_s_ampa = asarray(Pe.s_ampa)

if N_i>0:
    _Pi_s_nmda_in = asarray(Pi.s_nmda_in)
    _Pi_s_ampa_in = asarray(Pi.s_ampa_in)

_Pe_v = asarray(Pe.v)

if matrixStructure=='sparse':
    #_Wee = csr.csr_matrix(Cee.W.todense().transpose())
#    _Cee_W_rowj = asarray(Cee.W.rowj)
#    _Cei_W_rowj = asarray(Cei.W.rowj)
#    _Cee_W_rowdata =  asarray(Cee.W.rowdata)
#    _Cei_W_rowdata =  asarray(Cei.W.rowdata)
#    _Wei = asarray(Cei.W)
    _sparseWee = csr_matrix((Cee.W.alldata, Cee.W.allj, Cee.W.rowind), (N_e, N_e)).transpose()
    if N_i>0:
        _sparseWei = csr_matrix((Cei.W.alldata, Cei.W.allj, Cei.W.rowind), (N_e, N_i)).transpose() 
else:
    _Wee = asarray(Cee.W.transpose())
    if N_i>0:
        _Wei = asarray(Cei.W.transpose())
    
#tmp = loadmat(os.path.join('..','data','Cee.mat'))
#tmp = tmp['Cee']
#connect(P,Q,W)
        

#u = arange(0,simtime,defaultclock.dt)
#x = arange(0,simtime,defaultclock.dt)

str_f =''

# special case: no need to compute dot product
if binaryWeight and (not useSTDP) and sparseness_ee==1.0 and (N_i==0 or sparseness_ei==1.0):
#        str_f += 'x = exp(_logx)\n'
    for receptor in ['nmda','ampa']:
        str_f += 'tmp = sum(_Pe_s_'+receptor
        if useSTP:
            if tau_f>0:
                str_f += '*_u'
            if tau_d>0:
                str_f += '*_x'
        str_f += ')\n'
        for pop in populations:
            str_f += '_P'+pop+'_s_'+receptor+'_in[:] = w_e'+pop+' * tmp\n'
else:
#        str_f += 'x = exp(_logx)\n'
    for receptor in ['nmda','ampa']:
        str_f += 'tmp = _Pe_s_'+receptor
        if useSTP:
            if tau_f>0:
                str_f += '*_u'
            if tau_d>0:
                str_f += '*_x'            
        str_f += '\n'
        for pop in populations:
            if matrixStructure=='sparse':
                str_f += '_P'+pop+'_s_'+receptor+'_in[:] = _sparseWe'+pop+' * tmp\n'
            else:
                str_f += '_P'+pop+'_s_'+receptor+'_in[:] = dot(_We'+pop+',tmp)\n'
  
# kick
if useKick:
    str_f += 'if defaultclock.t==kickTime:\n'
    str_f += '    _Pe_v[rand(N_e)<kickStrength] = Vth_e + 1.0*mV'

print str_f   

@network_operation(when='start') 
def f():
    exec(str_f)    
                   


timedLog('Running ...') 

time_start  = time.time() 
NScount = 0
#g_nmda_e = 0*nS

if onlyRecord60:
    recordedGroup_e = Pe.subgroup(min(60,N_e))
    if N_i>0:
        recordedGroup_i = Pi.subgroup(min(60,N_i))
else:
    recordedGroup_e = Pe
    if N_i>0:
        recordedGroup_i = Pi

for i in arange(0,ceil(simTime/splitTime)):
    
    # reset spike monitors for dumping
    dump_e = SpikeMonitor(recordedGroup_e) 
    if N_i>0:
        dump_i = SpikeMonitor(recordedGroup_i)
    if dump_a:
        dump_a_monitor = StateMonitor(Pe,'Ca',record=True,timestep=iround(10*ms/defaultclock.dt))
    
#    if i==5:
#        print 'Changing dt'
#        get_default_clock().set_dt(0.1*ms/2/2/2/2/2)

    if useSTDP:
        tmp = Cee.W.todense()
        print 'Convergence index = %.3f' % mean(abs( tmp/wmax - (tmp>.5*wmax) ))
        del tmp
        outputdata={'Cee':Cee.W.todense()}
        savemat(os.path.join('..','data','Cee.seed.%09d' % randomSeed + '.t=%06.1f' % (splitTime*i)   +'s.mat'),outputdata,oned_as='row')
        del outputdata

    defaultclock.reinit(splitTime*i) # avoids drifting
    endTime = min(simTime,splitTime*(i+1))
    if endTime > simTime-monitorTime and (i==0 or splitTime*i <= simTime-monitorTime) :
        print 'Start monitoring'
        # initiating monitors 
        Me = SpikeMonitor(Pe) 
        if N_i>0:
            Mi = SpikeMonitor(Pi)
            rate_i=PopulationRateMonitor(Pi,bin=4*ms)
            pot_i = StateMonitor(Pi,'v',record=[0],timestep=iround(.5*ms/defaultclock.dt))
            s_gaba = StateMonitor(Pe,'s_gaba',record=[0,1],timestep=iround(1*ms/defaultclock.dt))
        rate_e=PopulationRateMonitor(Pe,bin=4*ms)
        
        pot_e = StateMonitor(Pe,'v',record=True,timestep=iround(10*ms/defaultclock.dt))
        s_nmda_in = StateMonitor(Pe,'s_nmda_in',record=[0,1],timestep=iround(1*ms/defaultclock.dt))
        s_ampa = StateMonitor(Pe,'s_ampa_in',record=[0,1],timestep=iround(1*ms/defaultclock.dt))
        if useAdaptation:
            sm = StateMonitor(Pe,'Ca',record=range(N_e),timestep=iround(50*ms/defaultclock.dt))
        if useSTP:
            if tau_f>0:
                u = StateMonitor(Pe,'u_f',record=range(N_e),timestep=iround(5*ms/defaultclock.dt))
            if tau_d>0:
                x = StateMonitor(Pe,'x_d',record=range(N_e),timestep=iround(5*ms/defaultclock.dt))
        if useAdaptiveThr:
            th1 = StateMonitor(Pe,'theta_1',record=True,timestep=iround(1*ms/defaultclock.dt))
            th2 = StateMonitor(Pe,'theta_2',record=True,timestep=iround(1*ms/defaultclock.dt))
            
#    run(endTime-splitTime*i,report='text')
    run(endTime-splitTime*i)
    
    timedLog('simulated time = %.1f' % (endTime) + ' s')
    print 'Avg E firing rate = %.2f' % (dump_e.nspikes/(endTime-splitTime*i)*(onlyRecord60*1.0/min(N_e,60)+(1-onlyRecord60)*1.0/N_e))
    if N_i>0:
        print 'Avg I firing rate = %.2f' % (dump_i.nspikes/(endTime-splitTime*i)*(onlyRecord60*1.0/min(N_i,60)+(1-onlyRecord60)*1.0/N_i))
    
    filepath = os.path.join('..','data','spike.%09d.' % randomSeed + '%05d.mat'  % i )
    print 'Saving spikes in ' + filepath
    if N_i>0:
        outputdata={'spike_e':dump_e.spikes,'spike_i':dump_i.spikes}
    else:
        outputdata={'spike_e':dump_e.spikes}
    if dump_a:
        outputdata['adaptation']= mean(dump_a_monitor.values[:,:],0)
        
    savemat(filepath,outputdata,oned_as='row')
    del outputdata 
    
#    # Exit conditions

#    if local_Me.nspikes==0:
#        print 'It seems that activity stopped. Exiting.'
#        break
    if endTime >= 0*second:
        meanRate = dump_e.nspikes*1.0/(endTime-splitTime*i)*(onlyRecord60*1.0/60+(1-onlyRecord60)*1.0/N_e)
        if meanRate > 5*Hz: # some activity
            tmp = array(dump_e.spikes);
            h = histogram(tmp[:,1], bins=arange(0,tmp[-1,1],10e-3))
            maxRate = h[0].max()/60.0/(10e-3*second)
            if maxRate<5*meanRate: # no "real" NS
                print 'Mean rate = %f Hz' % meanRate
                print 'Max rate = %f Hz' % maxRate
                print 'It seems that we have reached a stationary regime with high activity but no NS. Exiting.'
                break
            
    # count NS
    thr = 10*Hz
    timeBin = 50e-3*second
    tmp = array(dump_e.spikes);
    h = histogram(tmp[:,1], bins=arange(tmp[1,1],tmp[-1,1],timeBin))
    rate = h[0]/60.0/timeBin
    i_start = find( (rate[1:] >= thr) * (rate[0:-1] < thr) )
    #i_end = find( (rate[1:] <= thr) * (rate[0:-1] > thr) )
    if len(i_start>0): # at least one NS
        if len(i_start)==1:
            NScount+=1;
        else: # len > 1
            isi = ( i_start[1:] - i_start[0:-1] )*timeBin
            NScount+= 1 + sum(isi>6*second)
        print 'NScount = %d ' % NScount
        
        if NScount>=300: # enough data
            print 'We have enough NS data. Exiting.'
            break


#g_nmda_e = 0.165*nS
#run(simtime)
#g_nmda_e = 0*nS
#run(simtime)

time_end  = time.time()
running_time = time_end-time_start
print 'Simulation finished in %d seconds' % running_time

#print 'Saving weight matrices...'
#if N_i>0:
#    outputdata={'Cee':Cee.W.todense(),'Cii':Cii.W.todense(),'Cei':Cei.W.todense(),'Cie':Cie.W.todense()}
#else:
#    outputdata={'Cee':Cee.W.todense()}    
#savemat(os.path.join('..','data','weights.mat'),outputdata,oned_as='row')
#del outputdata




print 'For excitatory neurons:'
print 'Mean NMDA current = %.3f nA' % (-mean(Pe.s_nmda_in*g_nmda_e*(Pe.v-E_nmda)/(1+exp(-0.062*Pe.v/(1*mV))/3.57))/nA)
print 'Mean AMPA current = %.3f nA' % (-mean(Pe.s_ampa_in*g_ampa_e*(Pe.v-E_ampa))/nA)
if N_i>0:
    print 'Mean GABA current = %.3f nA' % (-mean(Pe.s_gaba*g_gaba_e*(Pe.v-E_gaba))/nA)
if useAdaptation:
    print 'Mean Ca current = %.3f nA' % (-mean(g_AHP*Pe.Ca*(Pe.v-E_K))/nA)
print 'Mean leak current = %.3f nA' % (-mean(gl_e*(Pe.v-El_e))/nA)
print 'Mean Contrib NMDA / Contrib AMPA %.2f ' % (mean(Pe.s_nmda_in*g_nmda_e*(Pe.v-E_nmda)/(1+exp(-0.062*Pe.v/(1*mV))/3.57))/mean(Pe.s_ampa_in*g_ampa_e*(Pe.v-E_ampa)))

if N_i>0:
    print 'For inhibitory neurons:'
    print 'Mean NMDA current = %.1f nA' % (-mean(Pi.s_nmda_in*g_nmda_i*(Pi.v-E_nmda)/(1+exp(-0.062*Pi.v/(1*mV))/3.57))/nA)
    print 'Mean AMPA current = %.1f nA' % (-mean(Pi.s_ampa_in*g_ampa_i*(Pi.v-E_ampa))/nA)
    print 'Mean GABA current = %.1f nA' % (-mean(Pi.s_gaba*g_gaba_i*(Pi.v-E_gaba))/nA)
    if useAdaptation:
        print 'Mean Ca current = %.3f nA' % (-mean(g_AHP*Pi.Ca*(Pi.v-E_K))/nA)
    print 'Mean leak current = %.1f nA' % (-mean(gl_i*(Pi.v-El_e))/nA)

timedLog('Done.')
t=time.localtime()

if monitorTime>0*second:
    execfile('plotResult.py')






