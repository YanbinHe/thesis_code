
function value = MinMaxCut(label,G,Db,Du)
%-------------------------------------------------------------------------%
% this function calculates the sum of the function
% cut(C_m)/(vol(C_m)-cut(C_m)). Here C_m is one of the clusters

% Input: label-partitioning results; A-adjacency matrix of the graph

% Output: the summation of MinMaxCut
%-------------------------------------------------------------------------%
M = length(unique(label)); % the number of clusters

% initializating cut and vol
cut = zeros(M,1);
den = zeros(M,1);
numP = size(1*full(G.W),1);
all1 = ones(1,numP);
sum_ratio = zeros(M,1);
L = diag(full(G.W)*ones(numP,1)) - full(G.W);
% calculate cut and vol for each clusters
for idxC = 1:M
    
    idxN = label; % index of each node
    idxN(idxN ~= idxC) = 0; % assign 0 to the nodes not in the cluster idxC(index of cluster)
    idxN = idxN/idxC; % cluster index is larger than 1, normalize it to 1
    idxN_inverse = all1 - idxN;% extract all the nodes not in this cluster (make 1 to 0 and 0 to 1)
    
    cut(idxC) = idxN'*L*idxN;%full(G.L)
    den(idxC) = idxN'*(Db*1*full(G.W)*Du)*idxN;
%     den(idxC) = idxN'*full(G.W)*idxN;
    sum_ratio(idxC) = cut(idxC)/den(idxC);
end

value = sum(sum_ratio);

end