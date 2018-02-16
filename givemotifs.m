%function c = givemotifs(M)
%
%  Calculates the number of connected triple motifs in the network (see
%  Milo et al. 2002)
%
%  Input:
%    M - N x N connectivity matrix
%
%  Output:
%    c - a 1 x 13 vector of the number of motifs found
%
%  Tuomo Mäki-Marttunen
%  Last modified 8.1.2013

function c=givemotifs(M)

N = size(M,1);
c = zeros(1,13); %number of found motifs 1-13
ctot = 0; %total number of found motifs

%Check all the triples (v1,v2,v3) (order not respected), nchoosek(N,3) in total
for v=1:N 
    for v2=v+1:N
        for v3=v2+1:N
            n = M(v,v2) + M(v2,v) + M(v,v3) + M(v3,v) + M(v2,v3) + M(v3,v2); %number of connections in total in the triple
            
            %for easiness, divide into cases by the total number of edges in the motif
            if n<2
                continue %can not be connected, hence ignored
            elseif n==2
                if M(v,v2) && M(v2,v) || M(v,v3) && M(v3,v) || M(v2,v3) && M(v3,v2)
                    continue; %only one bidirected link, hence not connected, hence ignored
                elseif M(v,v2) && M(v,v3) || M(v2,v) && M(v2,v3) || M(v3,v) && M(v3,v2)
                    ctot = ctot+1;
                    c(1) = c(1) + 1; %divergent motif
                elseif M(v,v2) && M(v2,v3) || M(v,v3) && M(v3,v2) || M(v2,v3) && M(v3,v) || M(v2,v) && M(v,v3) || M(v3,v2) && M(v2,v) || M(v3,v) && M(v,v2)
                    ctot = ctot+1;
                    c(2) = c(2) + 1; %chain motif
                else %if M(v2,v) && M(v3,v) || M(v,v2) && M(v3,v2) || M(v,v3) && M(v2,v3)
                    ctot = ctot+1; 
                    c(4) = c(4) + 1; %convergent motif
                end
            elseif n==3
                ctot = ctot + 1;
                if M(v,v2) && M(v2,v) && (M(v,v3) || M(v2,v3)) || ... 
                   M(v,v3) && M(v3,v) && (M(v,v2) || M(v3,v2)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v2,v) || M(v3,v))
                    c(3) = c(3) + 1; %bidirected link and a divergent link
                elseif M(v,v2) && M(v2,v) && (M(v3,v) || M(v3,v2)) || ...
                   M(v,v3) && M(v3,v) && (M(v2,v) || M(v2,v3)) || ...
                   M(v2,v3) && M(v3,v2) && (M(v,v2) || M(v,v3))
                    c(7) = c(7) + 1; %bidirected link and a convergent link
                elseif M(v,v2) && M(v2,v3) && M(v3,v) || M(v,v3) && M(v3,v2) && M(v2,v)
                    c(9) = c(9) + 1; %directed loop
                else
                    c(5) = c(5) + 1; %feed-dorward loop
                end
            elseif n==4
                ctot = ctot + 1;
                if ~M(v,v2) && ~M(v2,v) || ~M(v,v3) && ~M(v3,v) || ~M(v2,v3) && ~M(v3,v2)
                    c(8) = c(8) + 1; %one link missing
                elseif M(v,v2) && M(v2,v) && M(v,v3) && M(v2,v3) || ...
                   M(v,v3) && M(v3,v) && M(v,v2) && M(v3,v2) || ...
                   M(v2,v3) && M(v3,v2) && M(v2,v) && M(v3,v)
                    c(6) = c(6) + 1; %bidirected link and a mutual convergence to the remaining node
                elseif M(v,v2) && M(v2,v) && M(v3,v) && M(v3,v2) || ...
                   M(v,v3) && M(v3,v) && M(v2,v) && M(v2,v3) || ...
                   M(v2,v3) && M(v3,v2) && M(v,v2) && M(v,v3)
                    c(11) = c(11) + 1; %bidirected link and a mutual divergence from the remaining node
                else
                    c(10) = c(10) + 1; %bidirected link and a chain
                end
            elseif n==5
                ctot = ctot+1;
                c(12) = c(12) + 1; %the only motif with 5 edges
            else
                ctot = ctot+1;
                c(13) = c(13) + 1; %the only motif with 6 edges
            end
        end
    end
end






