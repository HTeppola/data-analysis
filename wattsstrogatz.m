%function M=wattsstrogatz(pos,p,q)
%
%  Creates a network that has topology between locally connected and
%  random. Based on (Watts & Strogatz, 1998: Collective dynamics of
%  'small-world' networks, Nature).
%
%  Input:
%    pos - N x dim matrix of position of nodes, OR, if a number N is given
%          instead, then the N nodes are given positions in a ring
%    p - connection probability OR a vector of length N denoting the
%        probabilities of each number inputs
%    q - the rewiring probability [0,1] (note that for large q asymmetry
%        may arise near the diagonal between the above-diagonal and
%        below-diagonal elements)
%
%  Output:
%    M - the connectivity matrix such that M(i,j) denotes the existence of
%        an edge from i to j
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function M=wattsstrogatz(pos,p,q)

if nargin < 1 || isempty(pos)
    pos = 100;
end
if numel(pos)==1
    N = pos;
    pos = [cos(2*pi*(1:N)'/N) sin(2*pi*(1:N)'/N)];
else
    N = size(pos,1);
end

if nargin < 2 || isempty(p)
    p = 0.1;
end

if nargin < 3 || isempty(q)
    q = 0.1;
end

if numel(p)~=1
  if numel(p) ~= N, error('The vector p should be of size N!'); end
  %p(i) gives the probability that number of in-neighbours is i-1 (i=1,...,N)
  p = p/sum(p);
  pcs = cumsum(p);
end 

M = zeros(N,N);
D = zeros(N,N); %distance matrix (.^2)
for i=1:size(pos,2)
  D = D+(pos(:,i)*ones(1,N) - ones(N,1)*pos(:,i)').^2; 
end

% 1. The generation of locally connected network with given in-degree:
for i=1:N %go through all nodes in order
    if numel(p)==1
        nn = binornd(N-1,p);
    else
      [vain,nn] = max(rand() <= pcs);  %[~,nns(i)] not supported in all versions 
      nn = nn - 1;
    end
    if nn==0
        continue;
    end
    pind = [1:i-1 i+1:N]; %possible inputs are all but the node itself
    [dists,sorted] = sort(D(i,pind));
    inds_equallyfar = find(dists == dists(nn));
    nsofar = sum(sum(M));
    if length(inds_equallyfar) == 1 % if a unique farthest node to be chosen as input
        M(pind(sorted(1:nn)),i) = 1; %choose as inputs all from closest to the farthest-to-be-chosen
    else
        M(pind(sorted(1:min(inds_equallyfar)-1)),i) = 1; %choose each nearer than farthest-to-be-chosen
        r = randperm(length(inds_equallyfar)); % choose randomly between the ones that as as far as farthest-to-be-chosen
        M(pind(sorted(min(inds_equallyfar)-1+r(1:nn-min(inds_equallyfar)+1))),i) = 1;
    end
end

% 2. Watts-Strogatz perturbation:
for i=1:N
    A = find(M(:,i)); % find the in-neighbours of i
    for j=1:length(A)
        if rand() < q
            freeind = ~M(:,i); %possible new candidates are all the ones not yet outputting to i (excluding i itself)
            freeind(i) = 0;
            freeind(A(j)) = 1;
            B = find(freeind);
            r = ceil(rand()*length(B));
            M(A(j),i) = 0;
            M(B(r),i) = 1;
        end
    end
end






