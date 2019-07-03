%function [spikes,Vrec]=runhh(Ne,M,T,dt,I_mu,I_dev)
%
%  Simulates the activity in a given network using the HH-type model with
%  four ionic and three synaptic currents (Golomb et al. 2006) and synaptic
%  depression (Golomb & Amitai 1997). Euler-Maruyama method is used for
%  integration.
%
%  Input:
%    Ne    - number of excitatory cells
%    M     - connectivity matrix with synaptic weights
%    T     - simulation time (ms)
%    dt    - simulation time step (ms)
%    I_mu  - constant current injection to each cell (uA/cm2)
%    I_dev - deviation of the Brownian current injection to each cell (uA/cm2)
%    
%  Output:
%    spikes - Nx2 matrix of the spiking data, where first column shows the
%             time of spike and second column the index of the neuron that spiked
%    Vrec   - the membrane potential time course of the first excitatory and
%             first inhibitory neuron (do not use if only excitatory neurons)
%
%  Tuomo M�ki-Marttunen
%  Last modified 8.1.2013

function [spikes,Vrec]=runhh(Ne,M,T,dt,I_mu,I_dev)

if nargin < 1 || isempty(Ne)
    Ne = 80; %number of excitatory cells
end

if nargin < 2 || isempty(M)
    Ni = 20; %number of inhibitory ceclls
    N = Ne + Ni;
    M = rand(N,N) < 0.5; %the connectivity matrix is picked by random
    M(~~eye(N)) = 0;     %no autapses
else
    N = size(M,1);
    Ni = N - Ne;
end

if nargin < 3 || isempty(T)
    T = 100; %ms
end

if nargin < 4 || isempty(dt)
    dt = 0.002; %ms
end

if nargin < 5 || isempty(I_mu)
    I_mu = 0; %mS/cm^2
end

if nargin < 6 || isempty(I_dev)
    I_dev = 1; %mS/cm^2
end

spike_thr = -30; %mV; membrane potential peaks that exceed this threshold are considered spikes
N_spikes_max = 500000; %max. number of spikes stored in memory

%The parameters for membrane potential dynamics taken from Golomb et al.
%2006, Journal of Neurophysiology.
C = 1; %microF/cm^2
C_I = 1; %microF/cm^2
g_Na = 35; %mS/cm^2
g_NaP = 0.2;
V_Na = 55; % mV
theta_m = -30; %mV
sigma_m = 9.5; %mV
theta_h = -45; %mV
sigma_h = -7; %mV
theta_th = -40.5; %mV
sigma_th = -6; %mV

theta_p = -47; %mV
sigma_p = 3; %mV

g_Kdr = 6; %mS/cm^2
%Note: in Golomb et al. 2006 g_Kdr is said to be 3mS/cm^2, but in the
%corresponding modelDB simulator it is set 6mS/cm^2. 6mS/cm^2 seems to 
%give the right spike shape.
V_K = -90; %mV
theta_n = -33; %mV
sigma_n = 10; %mV
theta_tn = -27; %mV
sigma_tn = -15; %mV

g_K_slow = 1.8; %mS/cm^2
theta_z = -39; %mV
sigma_z = 5; %mV
tau_z = 75; %ms

g_L = 0.05; %mS/cm^2
V_L = -70; %mV

theta_s = -20; %mV
sigma_s = 2; %mV
k_fP = 1; %(ms)^-1
g_AMPA = 0.08;
g_EI_AMPA = 0.2;
tau_AMPA = 5; %ms
V_Glu = 0; %mV

k_xN = 1; %(ms)^-1
tau_hat_NMDA = 14.3; %ms
k_fN = 1; %(ms)^-1
g_NMDA = 0.07;
g_EI_NMDA = 0.05;
tau_NMDA = 100; %ms
sigma_NMDA = 10; %mV
MgCon = 0.7; %mM (0.7 mM used as a typical Mg-concentration in cultures)
theta_NMDA = 10.5*log(MgCon/38.3); %mV, MgCon in mM
k_t = 1; %(ms)^-1
k_v = 0.001; %(ms)^-1

g_IE_GABAA = 0.05; %mS/cm^2
g_II_GABAA = 0.05; %mS/cm^2; different from (Golomb et al. 2006), where I-I synapses were not considered
V_GABAA = -70; %mV

g_I_Na = 35; %mS/cm^2
V_I_Na = 55; %mV
g_I_Kdr = 9; %mS/cm^2
V_I_Kdr = -90; %mV

g_I_L = 0.1; %mS/cm^2
V_I_L = -65; %mV


k_fA = 1; %(ms)^-1
tau_GABAA = 10; %ms

V = -71*ones(Ne,1); %the excitatory neuron membrane potential is at rest at around -71
h = 0.98+zeros(Ne,1);
n = 0.022+zeros(Ne,1);
z = 0+zeros(Ne,1);
V_I = -64*ones(Ni,1); %the inhibitory neuron membrane potential is at rest at around -64
h_I = 0.78+zeros(Ni,1);
n_I = 0.089+zeros(Ni,1);

%The synaptic gating variables:
s_AMPA = 0+zeros(Ne,1);
s_NMDA = 0+zeros(Ne,1);
x_NMDA = 0+zeros(Ne,1);
s_GABAA = 0+zeros(Ni,1);
%The glutamatergic resource variable (see Golomb & Amitai 1997):
T_Glu = 1+zeros(Ne,1);

spikes = zeros(N_spikes_max,2);
N_spikes = 0;
spike_check_interval = 10; %the interval (#time steps) to check spikes
Vrec_save_interval = 1;    %the interval (#spike check steps) to save the membrane potential data

if nargout > 1 %save membrane potential of first cells only if asked for
    Vrec = zeros(2,ceil(T/dt/spike_check_interval/Vrec_save_interval));
end

Vold = [V; V_I]*ones(1,2); %always save the past two membrane potential values to find a spike

for ti = 1:T/dt
    I_app = I_mu + I_dev/sqrt(dt)*randn(N,1); 
    % NOTE: _divided_ by the sqrt of dt. This way, as dV = ... + dt*I_app/C, we have the coefficient of dW proportional to sqrt(dt).
    % Hence, the variance of the random increment is proportional to dt, as required by the Brownian motion framework.

    I_Na = g_Na./(1+exp(-(V-theta_m)/sigma_m)).^3.*h.*(V-V_Na);
    I_NaP = g_NaP./(1+exp(-(V-theta_p)/sigma_p)).*(V-V_Na);
    I_Kdr = g_Kdr*n.^4.*(V-V_K);
    I_K_slow = g_K_slow*z.*(V-V_K);
    I_L = g_L*(V-V_L);
    
    alpha_m = 0.5*(V_I+35)./(1-exp(-(V_I+35)/10));
    beta_m = 20*exp(-(V_I+60)/18);
    I_I_Na = g_I_Na.*(alpha_m./(alpha_m + beta_m)).^3.*h_I.*(V_I-V_I_Na);
    I_I_Kdr = g_I_Kdr*n_I.^4.*(V_I-V_I_Kdr);
    I_I_L = g_I_L*(V_I-V_I_L);

    I_AMPA = g_AMPA*(V-V_Glu).*(M(1:Ne,1:Ne)'*s_AMPA); % AMPA current to excitatory cells
    I_EI_AMPA = g_EI_AMPA*(V_I-V_Glu).*(M(1:Ne,Ne+1:N)'*s_AMPA); % AMPA current to inhibitory cells

    I_NMDA = g_NMDA./(1+exp(-(V-theta_NMDA)/sigma_NMDA)).*(V-V_Glu).*(M(1:Ne,1:Ne)'*s_NMDA); % NMDA current to excitatory cells
    I_EI_NMDA = g_EI_NMDA./(1+exp(-(V_I-theta_NMDA)/sigma_NMDA)).*(V_I-V_Glu).*(M(1:Ne,Ne+1:N)'*s_NMDA); % NMDA current to inhibitory cells

    I_IE_GABAA = g_IE_GABAA*(V-V_GABAA).*(M(Ne+1:N,1:Ne)'*s_GABAA); % GABA current to excitatory cells
    I_II_GABAA = g_II_GABAA*(V_I-V_GABAA).*(M(Ne+1:N,Ne+1:N)'*s_GABAA); % GABA current to inhibitory cells
    
    dh = (1./(1+exp(-(V-theta_h)/sigma_h)) - h)./(0.1+0.75./(1+exp(-(V-theta_th)/sigma_th)));
    dn = (1./(1+exp(-(V-theta_n)/sigma_n)) - n)./(0.1+0.5./(1+exp(-(V-theta_tn)/sigma_tn)));
    dz = (1./(1+exp(-(V-theta_z)/sigma_z)) - z)/tau_z;
    
    alpha_h =  0.35*exp(-(V_I+58)/20);
    beta_h = 5./(1+ exp(-(V_I + 28)/10));
    dh_I = (alpha_h.*(1 - h_I) - beta_h.*h_I);
    alpha_n =  0.05*(V_I+34)./(1-exp(-(V_I+34)/10));
    beta_n = 0.625*exp(-(V_I + 44)/80);
    dn_I = (alpha_n.*(1-n_I) - beta_n.*n_I);
    
    ds_AMPA = (k_fP./(1+exp(-(V-theta_s)/sigma_s)).*T_Glu.*(1-s_AMPA) - s_AMPA/tau_AMPA);
    dx_NMDA = (k_xN./(1+exp(-(V-theta_s)/sigma_s)).*(1-x_NMDA) - (1-1./(1+exp(-(V-theta_s)/sigma_s))).*x_NMDA/tau_hat_NMDA);
    ds_NMDA = (k_fN.*x_NMDA.*T_Glu.*(1-s_NMDA) - s_NMDA/tau_NMDA);
    ds_GABAA = (k_fA*1./(1+exp(-(V_I-theta_s)/sigma_s)).*(1-s_GABAA)-s_GABAA/tau_GABAA);
    dT_Glu = -k_t./(1+exp(-(V-theta_s)/sigma_s)).*T_Glu + k_v*(1-T_Glu);
    
    dV = 1/C*(-I_Na-I_NaP-I_Kdr-I_K_slow-I_L-I_AMPA-I_NMDA-I_IE_GABAA+I_app(1:Ne));
    dV_I = 1/C_I*(-I_I_Na-I_I_Kdr-I_I_L-I_EI_AMPA-I_EI_NMDA-I_II_GABAA+I_app(Ne+1:N));
    
    % Euler-Maruyama:
    V = V+dt*dV;
    h = h+dt*dh;
    n = n+dt*dn;
    z = z+dt*dz;
    V_I = V_I+dt*dV_I;
    h_I = h_I+dt*dh_I;
    n_I = n_I+dt*dn_I;
    s_AMPA = s_AMPA+dt*ds_AMPA;
    s_NMDA = s_NMDA+dt*ds_NMDA;
    x_NMDA = x_NMDA+dt*dx_NMDA;
    s_GABAA = s_GABAA+dt*ds_GABAA;
    T_Glu = T_Glu + dt*dT_Glu;
    
    %Check the spiking and save the membrane potential values if needed:
    if mod(ti,spike_check_interval) == 0
        %There was a spike iff there was a local maximum of membrane
        %potential and the maximum value is above spike_thr.
        spiking = [V; V_I] > spike_thr & Vold(:,2) < Vold(:,1) & Vold(:,1) >= [V;V_I];
        if any(spiking) && N_spikes < N_spikes_max
            spiking = find(spiking);
            spvec = [(ti-1)*dt*ones(length(spiking),1), spiking];
            spikes(N_spikes+1:min(N_spikes+length(spiking),N_spikes_max),:) = spvec(1:min(length(spiking),N_spikes_max-N_spikes),:);
            N_spikes = N_spikes + length(spiking);
        end      
        if nargout > 1 && mod(ti/spike_check_interval,Vrec_save_interval)==0
            Vrec(:,ti/spike_check_interval) = [V(1); V_I(1)];
        end
	    Vold = [[V;V_I] Vold(:,1)];
    end
    
end

end






