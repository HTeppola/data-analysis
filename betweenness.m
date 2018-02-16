%function betweennesses = betweenness(M,path_length_max,Nsamp)
%
%  Calculates the node-betweenness for the nodes of the network
%
%  Input:
%    M - N x N connectivity matrix
%    path_length_max - the maximum considered path length (20 is an ad hoc
%                      limit that should be enough for connectivities p larger than 0.1)
%    Nsamp - number of samples (by default N, i.e. betweenness calculated for all nodes)
%
%  Output:
%    betweennesses - the betweenness values in a 1 x Nsamp vector
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function betweennesses = betweenness(M,path_length_max,Nsamp)

N = size(M,1);

if nargin < 2 || isempty(path_length_max)
    path_length_max = 20;
end
if nargin < 3 || isempty(Nsamp)
    Nsamp = N;
end
testedInd = randperm(N);
testedInd = sort(testedInd(1:Nsamp));

betweennesses = zeros(1,Nsamp);

% Calculate the number of paths of length 1...path_length_max from each
% node to each node
Ms = zeros(N,N,path_length_max+1);
Ms(:,:,1) = M;
accounted_for = zeros(N);
pathl = zeros(N); %Save here the shortest path lengths, zero meaning that no path was found
for i=1:path_length_max
    pathl = pathl + i*(Ms(:,:,i) > 0 & ~accounted_for);
    accounted_for = accounted_for | Ms(:,:,i) > 0;
    Ms(:,:,i+1) = Ms(:,:,i)*M;
end

for iv=1:Nsamp
  v = testedInd(iv); %the index of the node whose betweenness is being calculated
  sumvalue = 0; %the betweenness value is saved here iteratively
    
  P = M([1:v-1 v+1:N],[1:v-1 v+1:N]); % M from which v removed

  % Calculate the number of paths of length 1...path_length_max from each
  % node (except v) to each node (except v)
  Ps = zeros(N-1,N-1,path_length_max);
  Ps(:,:,1) = P;
  accounted_for = zeros(N-1);
  pathlwv = zeros(N-1);
  for i=1:path_length_max
    pathlwv = pathlwv + i*(Ps(:,:,i) > 0 & ~accounted_for);
    accounted_for = accounted_for | Ps(:,:,i) > 0;
    Ps(:,:,i+1) = Ps(:,:,i)*P;
  end
    
  for i=[1:v-1 v+1:N] %Go through all "from"-nodes except v
    
    % save to indices the "to"-nodes to be checked
    if i<v
      indices = [1:i-1 i+1:v-1 v+1:N];
    else
      indices = [1:v-1 v+1:i-1 i+1:N];
    end
    
    %iwv is an index corresponding to i, but in P
    iwv = i;
    if i > v
      iwv = iwv-1;
    end
    
    for j=indices %go through all "to"-indices
      
      %jwv an index corresponding to j, but in P
      jwv = j;
      if j > v
        jwv = jwv-1;
      end
      if pathl(i,j) > 1 %if there is a path from i to j, iv may lie on it
        if pathlwv(iwv,jwv) == 0 || pathlwv(iwv,jwv) > pathl(i,j)
          %in this case there is a path from i to j when iv is part of the
          %graph, but no path when iv is not part of the graph, hence the
          %addition of one in the betweenness value
          sumvalue = sumvalue + 1;
        else
          %otherwise, there are some paths from i to j even without iv,
          %hence, the betweenness value is added only the fraction of paths
          %that iv lies on
          sumvalue = sumvalue + 1 - Ps(iwv,jwv,pathlwv(iwv,jwv))/Ms(i,j,pathl(i,j));
        end
      end
    end
  end
  betweennesses(iv) = sumvalue;
end







