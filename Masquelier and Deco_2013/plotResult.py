period = array([max(0,simTime-monitorTime),simTime])
#period = array([1.86*second,1.96*second])


fig = plt.figure()
fig.canvas.manager.set_window_title('Spikes - rand %d' % randomSeed +  ' - %02d'  % t[3] +  ':%02d' % t[4] + ':%02d' % t[5])


fig.add_subplot(1+(N_i>0),1,1)
raster_plot(Me)
title('Excitatory')
grid()
#axis([0, 1*simtime/ms, -1, 3])
#xlim(1/(1*ms)*period)

if N_i>0:
    fig.add_subplot(1+(N_i>0),1,2)
    raster_plot(Mi)
    title('Inhibitory')
    grid()
#xlim(1/(1*ms)*period)
#plot(s_nmda_in.times/ms, s_nmda_in[0],'r')
#plot(pot.times/ms, pot[0]/mV)
#axis([0, 1*simtime/ms, -75, -50])


fig = plt.figure()
fig.canvas.manager.set_window_title('Rates  - rand %d' % randomSeed +  ' - %02d' % t[3] +  ':%02d' % t[4] + ':%02d' % t[5])

nPlots = 1 + useAdaptation + useSTP + useAdaptiveThr + 1 + 1 -2
i = 0

i += 1
fig.add_subplot(nPlots,1,i)
plot(rate_e.times,rate_e.rate,'-g',label='Excitatory')
if N_i>0:
    plot(rate_e.times,rate_i.rate,'-r',label='Inhibitory')
    leg = legend(loc='upper right')
title('Firing rates')
#xlabel('t (s)')
ylabel('Hz')
xlim(period)
grid()

if useAdaptation:
    i += 1
    fig.add_subplot(nPlots,1,i)
    plot(sm.times,mean(sm.values[:,:],0),'-r')
    title('a')
    xlim(period)
    grid()
if useSTP:
    i += 1
    fig.add_subplot(nPlots,1,i)
    if tau_f>0:
        plot(u.times,mean(u.values[:,:],0),'-g',label='u')
    if tau_d>0:
        plot(x.times,mean(x.values[:,:],0),'-r',label='x')
    if tau_f * tau_d >0:
        plot(x.times,mean(x.values[:,:]*u.values[:,:],0),'-b',label='ux')    
    leg = legend(loc='upper right')
    title('STP')
    xlabel('t (s)')
    ylim(0,1)
    xlim(period)
    grid()
if useAdaptiveThr:
    i += 1
    fig.add_subplot(nPlots,1,i)
    plot(th1.times,Vth_e+mean(th1.values[:,:]+th2.values[:,:],0),'-b')
    title('Avg threshold')
    xlim(period)
    grid()

#i += 1
#fig.add_subplot(nPlots,1,i)
#plot(s_nmda_in.times,s_nmda_in.values.transpose(),'-',label='nmda')
#plot(s_ampa.times,s_ampa.values.transpose(),'--',label='ampa')
#if N_i>0:
#    plot(s_gaba.times,s_gaba.values.transpose(),':',label='gaba')
#leg = legend(loc='upper right')
#title('Input')
#xlabel('t (s)')
##ylim(150,600)
#xlim(period)
#grid()
#
#i += 1
#fig.add_subplot(nPlots,1,i)
#plot(pot_e.times,mean(pot_e.values.transpose(),1),'-',label='Exc')
##plot(pot_i.times,pot_i.values.transpose(),'--',label='Inh')
#title('Pot')
#xlabel('t (s)')
#xlim(period)
##ylim(150,600)
#leg = legend(loc='upper right')
#grid()

outputdata={'pot_e':mean(pot_e.values[:,:],0),'rate_time':rate_e.times,'rate':rate_e.rate}
if N_i>0:
    outputdata['pot_i']=pot_i.values.transpose()
if useAdaptation:
    outputdata['adapt_time'] = sm.times
    outputdata['adapt'] = mean(sm.values[:,:],0)
if useSTP:
    if tau_f>0:
        outputdata['u_time'] = u.times
        outputdata['u'] = mean(u.values[:,:],0)        
    if tau_d>0:
        outputdata['x_time'] = x.times
        outputdata['x'] = mean(x.values[:,:],0)
savemat(os.path.join('..','data','dump.%09d' % randomSeed + '.mat'),outputdata,oned_as='row')


fig = plt.figure()
fig.canvas.manager.set_window_title('Sample spikes  - rand %d' % randomSeed +  ' - %02d' % t[3] +  ':%02d' % t[4] + ':%02d' % t[5])
import random
sample_size = 60
sample_e = random.sample(range(N_e),iround(sample_size*.8))
if N_i>0:
    sample_i = random.sample(range(N_i),iround(sample_size*.2))
for n in range(sample_size):
    if n%5==0 and N_i>0:
        plot(Mi[sample_i[n/5]],n*ones(len(Mi[sample_i[n/5]])),'.')
    else:
        plot(Me[sample_e[n-1-n/5]],n*ones(len(Me[sample_e[n-1-n/5]])),'.')
xlabel('t (s)')
grid()

fig = plt.figure()
fig.canvas.manager.set_window_title('NS profile  - rand %d' % randomSeed +  ' - %02d' % t[3] +  ':%02d' % t[4] + ':%02d' % t[5])
plot(rate_e.times,rate_e.rate,'-g',label='Excitatory')
rate_sum = N_e*1.0/(N_e+N_i)*rate_e.rate
if N_i>0:
    plot(rate_e.times,rate_i.rate,'-r',label='Inhibitory')
    rate_sum += N_i*1.0/(N_e+N_i)*rate_i.rate
    plot(rate_e.times,rate_sum,'-b',label='Weighted sum')
    ylim([ -5, max(rate_e.rate.max(),rate_i.rate.max())])
else:
    ylim([ -5, rate_e.rate.max()])  
xlim(period[0]+[ rate_sum.argmax()*rate_e._bin*defaultclock.dt-.150, rate_sum.argmax()*rate_e._bin*defaultclock.dt+.150 ])
leg = legend(loc='upper right')
title('Firing rates')
xlabel('t (s)')
ylabel('Hz')
grid()

    
#raster_plot(N) 
show()






