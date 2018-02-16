%function [BL,BS,start_t,BL_abs,mFr,Rs,Fs] = burstanalysis(firings,dt,maxISI,minSpikes,minElectrodes,ignoreBegin,sigma,divisionpars)
%
%  Analyses the bursting properties using the population spike train "firings"
%
%  Input:
%    firings - a Nx2 matrix; first column has the spike times (in ms) and the
%              second the corresponding neuron indices.
%    binwidth - the binsize for calculating the burst profile
%    maxISI - the maximum distance between spikes that can belong to the same
%             burst
%    minSpikes - the minimum number of spikes that a burst must have
%    minUnique - the minimum number of individual neurons that must spike
%                in a burst
%    ignoreBegin - ignore the spikes whose spike time is less than this value
%                  (in ms)
%    sigma - the st. deviation of Gaussian with which the profile is smoothed
%
%  Output:
%    BL - burst lengths
%    BS - burst sizes (how many spikes in a burst)
%    start_t - times of the first spike of the bursts
%    BL_abs - absolute burst lengths (duration from first to last spike)
%    mFr - maximum firing rates of the bursts
%    Rs - lengths of the rising slopes
%    Fs - lengths of the falling slopes
%
% Tuomo Mäki-Marttunen
% Last modified 8.1.2013

function [BL,BS,start_t,BL_abs,mFr,Rs,Fs] = burstanalysis(firings,dt,maxISI,minSpikes,minUnique,ignoreBegin,sigma)

if nargin < 2 || isempty(dt)
    dt = 1;
end
if nargin < 3 || isempty(maxISI)
    maxISI = 25;
end
if nargin < 4 || isempty(minSpikes)
    minSpikes = 40;
end
if nargin < 5 || isempty(minUnique)
    minUnique = 30; %the minimum number of spikers in a burst
end
if nargin < 6 || isempty(ignoreBegin)
    ignoreBegin = 0; %ignore x ms from the beginning
end
if nargin < 7 || isempty(sigma)
    sigma = 2.5; %ms
end

BS = [];
start_t = [];
BL_abs = [];
mFr = [];
Rs = [];
Fs = [];

firings(firings(:,1) < ignoreBegin,:) = [];
[~,ord] = sort(firings(:,1));
firings = firings(ord,:);
     
times = firings(:,1); %firings times (in ms)
dtimes = times(2:end)-times(1:end-1); %ISIs (in ms)
chind = [find(dtimes > maxISI); length(times)]; %the indices ahead of which lies and inter-burst period.
                                                %The last index is included so that the last burst is also included
                                                
sd = sigma/dt; %how many bins is the SD of the convoluting Gaussian
NemptySteps = round(sd*5); %consider the Gaussian up to +-5*sigma

for j=1:length(chind)-1
    chind_j = 0; %this makes it possible to study the first possible burst too. the first chind is thus 0
    if j > 1
        chind_j = chind(j-1);
    end
    SCnew = chind(j)-chind_j; %number of spikes in the considered burst
    nEl = length(unique(firings(chind_j+1:chind(j+1),2))); %number of unique spikers in the considered burst
    if SCnew < minSpikes || nEl < minUnique
        continue;
    end

    BL_absnew = times(chind(j)) - times(chind_j+1); % and burst length (in ms)

    mappi = zeros(1,floor(BL_absnew/dt)+1+2*NemptySteps); %the length of the burst profile is the number of
                                                          %bins required for the duration of the burst +-5*sigma
    mappiprofile = mappi; %mappi will have the raw global firing rate, and mappiprofile the smoothened
        
    %Save the number of spikes in bins of length dt to vector mappi:
    for i=1:length(mappi)
        mappi(i) = sum(1+floor((times(chind_j+1:chind(j)) - times(chind_j+1))/dt) == i-NemptySteps);
    end
        
    
    %Smooth the vector mappi using Gaussian of SD sigma and save to mappiprofile
    gw = 1/sqrt(2*pi*sd^2)*exp(-1/2/(sd^2)*(-length(mappi)+1:length(mappi)-1).^2);
    for i=1:length(mappi)
      mappiprofile(i) = mappi*gw((1:length(mappi))+length(mappi)-i)';
    end
        
    BS = [BS, SCnew]; %detect this as a burst, save spike count
    BL_abs = [BL_abs, BL_absnew]; % and absolute burst length (in ms)

    [a,peakind] = max(mappiprofile); %the height and index (in mappiprofile) of the peak
    mappiprofile = [0 mappiprofile 0]; %put zeros into the ends to make sure that a half-width exists
           
    % Calculate the widths of the rising slope and falling slope (both
    % using half (default) and quarter (application) mFr
    crossind1 = find(mappiprofile(1:end-1) < a/2 & mappiprofile(2:end) >= a/2);
    crossind2 = find(mappiprofile(1:end-1) >= a/2 & mappiprofile(2:end) < a/2);

    mFr = [mFr a/dt];
    Rs = [Rs dt*(peakind - (min(crossind1) + (abs(mappiprofile(min(crossind1))-a/2) > abs(mappiprofile(min(crossind1)+1)-a/2))))];
    Fs = [Fs dt*((max(crossind2) + (abs(mappiprofile(max(crossind2))-a/2) > abs(mappiprofile(max(crossind2)+1)-a/2)))-peakind)];
    %Fs = [Fs dt*((min(crossind2) + (abs(mappiprofile(min(crossind2))-a/2) > abs(mappiprofile(min(crossind2)+1)-a/2)))-peakind)];
    start_t = [start_t times(chind_j+1)]; %the start time of the burst
    if length(BL_abs)>50
    end
end
BL = Rs + Fs;

end
    







