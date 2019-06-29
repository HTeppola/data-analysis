%function cc = clusteringcoefficient(M)
%
%  Calculates the local clustering coefficients of the network
%
%  Input:
%    M - N x N connectivity matrix
%
%  Output:
%    cc - an N x 1 vector of the local clustering coefficient values
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function cc=clusteringcoefficient(M)

nn = sum(M|M')'; %the numbers of neighbours of each node (in or out or both)
cc = zeros(size(nn));
for i=1:length(nn)
  inds = find(M(i,:)'|M(:,i)); %the indices of neighbours of i
  nts = 0; %save the number of triangles here iteratively
  for j=1:length(inds) %go through the neighbours of i
    %check the links from inds(j) to other neighbours of i
    %(inds(k), where k<j). for each such link add to nts:
    %   1*1, if both the links (i-inds(j)) and (i-inds(k)) are unidirected
    %   1*2, if the link (i-inds(j)) unidirected and (i-inds(k)) bidirected
    %   2*1, if the link (i-inds(j)) bidirected and (i-inds(k)) unidirected
    %   2*2, if both the links (i-inds(j)) and (i-inds(k)) are biidirected
    nts = nts + (M(i,inds(j))+ M(inds(j),i))*sum(M(inds(j),inds([1:j-1, j+1:end])).*(M(i,inds([1:j-1,j+1:end])) + M(inds([1:j-1,j+1:end]),i)'));
  end
  cc(i) = nts/(nn(i)*(nn(i)-1)*4); %the maximum number of triangles is 8*nchoosek(nn(i),2)
end
  







