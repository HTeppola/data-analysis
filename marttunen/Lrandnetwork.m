%function M=Lrandnetwork(N,p,pow,rotlen)
%
%  Creates a network, where directed loops of length rotlen are promoted.
%
%  Input:
%    N - number of nodes
%    p - connection probability OR a vector of length N denoting the
%        probabilities of each number inputs
%    pow - the strength of the loop promotion [0,inf]
%    rotlen - the promoted loop length
%
%  Output:
%    M - the connectivity matrix such that M(i,j) denotes the existence of
%        an edge from i to j
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function M=Lrandnetwork(N,p,pow,rotlen)

  M = zeros(N,N);
  if nargin < 3 || isempty(pow)
    pow = 1;
  end
  if nargin < 4 || isempty(rotlen)
    rotlen = 3;
  end
  if numel(p)==1  
    nns = binornd(N-1,p,N,1);  
  elseif numel(p)==N
    nns = zeros(N,1);
    %p(i) gives the probability that number of in-neighbours is i-1 (i=1,...,N)
    p = p/sum(p);
    pcs = cumsum(p);
    for i=1:N
      [vain,nns(i)] = max(rand() <= pcs); %[~,nns(i)] not supported in all versions 
    end
    nns = nns - 1;
  end
  
  pinds = cell(N);
  for i=1:N
      pinds{i} = [1:i-1,i+1:N];
  end
  nplaced = zeros(N,1);
  
  lastinput = 0; %this tells which node was last chosen to be an input
  for ie=1:sum(nns) %go through all the edges in a non-obvious order
    
    if lastinput %if some node has already been chosen as an input
        if nplaced(lastinput) < nns(lastinput) %if this node has to have more inputs
            i = lastinput; %choose that one
        else %otherwise pick by random
            w = nns - nplaced;
            w = w/sum(w);
            r = rand(1) < cumsum(w);
            i = find(r,1);
        end
    else %otherwise pick by random
        w = nns - nplaced;
        w = w/sum(w);
        r = rand(1) < cumsum(w);
        i = find(r,1);
    end
    
    pind = pinds{i};
    Npaths = nminpath(M,i,rotlen); %a matrix that tells the number of paths from i to other nodes

    % The points are given such that the highest points (>4) will be for
    % those nodes, to which there is a path from i of length rotlen-1, but no path
    % shorter than that. Choosing a such node as input will close the loop
    % of length rotlen. The second largest points (>3) is given to those
    % nodes, that are not reachable from i with a path of length rotlen-1,
    % but on the other hand not with any shorter path either. Choosing such
    % a node will not close a loop of length <= rotlen. The second lowest
    % points (>2) is given to nodes that are reachable from i with rotlen-1
    % edges, but which are also reachable by shorter paths. Choosing such a
    % node will close a loop of length < rotlen. The least points (1) are given
    % to the remaining nodes.
    points = (Npaths(end,:)&~any(Npaths(1:end-1,:),1)).*(4+Npaths(end,:)./((N-2)/(rotlen-2))^(rotlen-2)) +...
             (~Npaths(end,:)&~any(Npaths(1:end-1,:),1))*3 +...
             (Npaths(end,:)&any(Npaths(1:end-1,:),1)).*(2+Npaths(end,:)./((N-1)^(rotlen-2))) +...
             (~Npaths(end,:)&any(Npaths(1:end-1,:),1));
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
    r = rand(1) < cumsum(weights');
    rind = find(r,1);
    M(pind(rind),i) = 1;
    lastinput = pind(rind);
    nplaced(i) = nplaced(i) + 1;
    pind(rind) = [];
    pinds{i} = pind;
  end
end

%function nminpath returns a matrix that tells by how many paths you can
% reach nodes 1...N (horizontal) using 1...rotlen-1 (vertical) edges from
% node i
function Npaths = nminpath(M,i,rotlen) 
  N = size(M,1);
  Npaths = zeros(rotlen-1,N);
  Npaths(1,:) = M(i,:);
  
  for j=2:rotlen-1
    Npaths(j,:) = Npaths(j-1,:)*M;
  end
end







