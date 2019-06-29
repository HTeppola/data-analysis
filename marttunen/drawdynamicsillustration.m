% Illustration of HH and LIF model output. Random networks of 20 inhibitory and 80 excitatory neurons
% with binomially distributed in-degree (p=0.2) used. Spiking activity is simulated for one minute,
% and the membrane potential of first excitatory and first inhibitory neuron are considered in more
% detail.
%
% The simulation of HH model is heavy, takes a couple of hours on
% a moderately fast PC. The LIF simulation is fast, takes around 10s. The
% code should be runnable as such in UNIX/LINUX systems with MATLAB
% R2010a/R2011b, Python 2.4.3, and NEST 2.0-0-rc4. In Windows systems
% the code can be divided to two sections (first one until saving pars.mat,
% second one from loading pars.mat on), and the PyNEST code is run between
% these two sections.
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

tic;
seed = 1;

rand('state',seed);
randn('state',seed);
Ne = 80;
Ni = 20;

M = FFrandnetwork(Ne+Ni,0.2,0);

[spikes_hh,Vrec_hh] = runhh(Ne,M*0.177,61000,0.0025,0,0.9);
disp(['HH solved, toc = ' num2str(toc) ' s']);
disp(' ');
save('example.mat','spikes_hh','Vrec_hh');

syndelay=0.2;  %synaptic delay cannot be smaller than time step
syndelaystd=0; %no variance in synaptic delay
Ie=12;
Ii=12;
dt=0.2;
t_sim=61000;
wee=17.21; 
wei=-3*17.21; 
wie=4*17.21;
wii=-4*17.21;
noisestd=7.3;

save('pars.mat','Ie','Ii','wee','wei','wie','wii','Ne','Ni','dt','syndelay','syndelaystd','noisestd','t_sim','M','seed');
% unix(['python runiafts_exp.py pars']); %This has to be changed if Windows is being used
% load('pars.mat','spikefile','vmfile');
% disp(['LIF solved, toc = ' num2str(toc) ' s']);
% disp(' ');

%spikes_lif = load(spikefile,'-ascii');
%Vrec_lif = load(vmfile,'-ascii');

save('example.mat','spikes_hh','Vrec_hh','spikes_lif','Vrec_lif');

%spikes_lif = spikes_lif(:,[2 1]);
%spikes_lif(:,2) = round(spikes_lif(:,2) - 100);
%Vrec_lif = Vrec_lif(:,3)';
%Vrec_lif(round(spikes_lif(spikes_lif(:,2)==1,1))) = 31;
Vrec_hh = Vrec_hh(1,:);

figure;
set(gcf,'position',[680   558   748   420]);
maxISI = 25;
minElectrodes = 0.3*80;

SP = zeros(2,3);
inset = zeros(2,1);
for imod=1:2
  SP(imod,1) = subplot('position',[0.1, 0.1+0.45*(imod-1), 0.2, 0.35]);
  SP(imod,2) = subplot('position',[0.38, 0.1+0.45*(imod-1), 0.4, 0.35]);
  SP(imod,3) = subplot('position',[0.8, 0.1+0.45*(imod-1), 0.18, 0.35]);
%   if imod==1
%     spikes = spikes_lif;
%     Vrec = Vrec_lif;
%     resolution = 1; %1 saved value/ms
%     minBS = 0.4*80;
%   else
    spikes = spikes_hh;
    Vrec = Vrec_hh;
    resolution = 40; %40 saved values/ms
    minBS = 1.5*80;
  end
  subplot(SP(imod,1));
  Vrec = Vrec(resolution*1000:resolution*61000-1); %ignore the first second
  spikes = spikes(spikes(:,1)>= 1000,:);           %ignore the first second
  spikes(:,1) = spikes(:,1)-1000;
  sp_time = spikes(find(spikes(:,2)==1 & spikes(:,1) > 29000,1,'first'),1); % time of the considered spike
  start_Vrec_t = round((sp_time-1500)*resolution);
  stop_Vrec_t = round((sp_time+500)*resolution);
  plot((0:1/resolution:2000)/1000,Vrec(start_Vrec_t:stop_Vrec_t),'k-');
  hold on;
  yrange = [min(Vrec(start_Vrec_t:stop_Vrec_t)), max(Vrec(start_Vrec_t:stop_Vrec_t))];
  axis([0, 2, yrange(1)-(yrange(2)-yrange(1))*0.1, yrange(2)+(yrange(2)-yrange(1))*0.1]);
  plot([1.420 1.420 1.580 1.580 1.420], [yrange(1) yrange(2) yrange(2) yrange(1) yrange(1)],'r--');
  xlabel('time (s)'); ylabel('(mV)      ','rot',0);
  box off;

  inset(imod) = axes('position',[0.11, 0.28+0.45*(imod-1), 0.11, 0.15]);
  box on;
  start_Vrec_t = round((sp_time-15)*resolution);
  stop_Vrec_t = round((sp_time+25)*resolution);
  plot((0:1/resolution:40),Vrec(start_Vrec_t:stop_Vrec_t),'k-');
  axis([0, 40, yrange(1)-(yrange(2)-yrange(1))*0.1, yrange(2)+(yrange(2)-yrange(1))*0.1]);
  set(gca,'xcolor',[1 0 0], 'ycolor',[1 0 0],'xtick',[0 20 40],'ytick',[],'fontsize',7);
  xlabel('(ms)')
  
  subplot(SP(imod,2));
  plot(spikes(:,1)/1000,spikes(:,2),'k.');
  hold on;
  axis([0 60 0 101]);
  plot(sp_time/1000,1,'r.');
  [BL_rel,BS,start_t,BL] = burstanalysis(spikes(spikes(:,2)<=80,:),1,maxISI,minBS,minElectrodes);
  start_B = start_t(3);
  end_B = start_B + BL(3);
  plot([start_B-500 start_B-500 end_B+500 end_B+500 start_B-500]/1000,[0.5 100.5 100.5 0.5 0.5],'r--');
  xlabel('time (s)');
  ylabel('Neuron index');
  
  subplot(SP(imod,3));
  inds = spikes(:,1) >= start_B-500 & spikes(:,1) < end_B+500;
  plot(spikes(inds,1)-start_B,spikes(inds,2),'k.');
  axis([-200 700 0 101]);
  xlabel('time (ms)');
  set(gca,'yticklabel',[],'xtick',[-300 0 300 600]);
%  ylabel('Neuron index');
  saveas(gcf,'examplefig.pdf');
end







