%function [path_lengths, lengths_to_self] = pathlength(M,path_length_max)
%
%  Calculates the path lengths between the nodes of the network
%
%  Input:
%    M - N x N connectivity matrix
%    path_length_max - the maximum considered path length (20 is an ad hoc
%      limit that should be enough for connectivities p larger than 0.1)
%
%  Output:
%    path_lengths - a 1 x path_length_max vector indicating the number of
%                   paths of length 1...path_length_max
%    lengths_to_self - a 1 x path_length_max vector indicating the number of
%                      paths to self of length 1...path_length_max
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function [path_lengths,lengths_to_self] = pathlength(M,path_length_max)

if nargin < 2 || isempty(path_length_max)
    path_length_max = 20;
end

N = size(M,1);

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

path_lengths = zeros(1,path_length_max);
lengths_to_self = zeros(1,path_length_max);
for i=1:path_length_max
    path_lengths(i) = sum(sum(pathl==i));
    lengths_to_self(i) = sum(diag(pathl==i));
end






