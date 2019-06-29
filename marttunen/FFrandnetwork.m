%function M=FFrandnetwork(N,p,pow)
%
%  Creates a (recurrent) network, where feed-forward loops are promoted.
%
%  Input:
%    N - number of nodes
%    p - connection probability OR a vector of length N denoting the
%        probabilities of each number inputs
%    pow - the strength of the feed-forward loop promotion [0,inf]
%
%  Output:
%    M - the connectivity matrix such that M(i,j) denotes the existence of
%        an edge from i to j
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function M=FFrandnetwork(N,p,pow)

  M = zeros(N,N);
  if nargin < 3 || isempty(pow)
    pow = 1;
  end

  if numel(p)~=1
    if numel(p) ~= N, error('The vector p should be of size N!'); end
    %p(i) gives the probability that number of in-neighbours is i-1 (i=1,...,N)
    p = p/sum(p);
    pcs = cumsum(p);
  end

  for i=1:N %go through all nodes in order
    if numel(p)==1
      nn = binornd(N-1,p);
    else
      [vain,nn] = max(rand() <= pcs);  %[~,nns] not supported in all versions 
      nn = nn - 1;
    end
    pind = [1:i-1,i+1:N]; %possible inputs are all but the node itself
    points = ones(N,1); % initially 1 point for each node
    lastinput = 0; %this tells which node was last chosen to be an input
    for j=1:nn %go through all possible inputs in order
        
      % instead of recalculating the weight for all nodes, the
      % recalculation can be done only for those nodes that project to the
      % last chosen input and do not yet project to i
      if lastinput && any(M(:,lastinput) & ~M(:,i))
        points = points + (M(:,lastinput) & ~M(:,i)); % points added to the nodes that are now
                                                       % reaching i with two hops
      end
      if ~isinf(pow)
        weights = points(pind).^pow;
      else
        weights = points(pind) == max(points(pind));
      end
      if any(isinf(weights)) %if the weight of any node is inf, the choice will be made by random between such nodes
         weights = isinf(weights);
      end
      if all(weights==0)
          weights = weights+1;
      end
      weights = weights/sum(weights);
      r = [0; rand(1) < cumsum(weights)];
      rind = ~~(r(2:end)-r(1:end-1));
      M(pind(rind),i) = 1;
      lastinput = pind(rind);
      pind(rind) = [];
    end
  end







