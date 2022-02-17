function value = MinMaxCut_gen(label,A)

% this function calculates the sum of the function
% cut(C_m)/(vol(C_m)-cut(C_m)).

M = length(unique(label)); % the number of clusters

% initializating cut and vol
cut = zeros(M,1);
den = zeros(M,1);
numP = size(A,1);
all1 = ones(1,numP);
sum_ratio = zeros(M,1);
% calculate cut and vol for each clusters
L = diag(A*all1') - A;

for idxC = 1:M
    
    idxN = label; % index of each node
    idxN(idxN ~= idxC) = 0; % assign 0 to the nodes not in the cluster idxC(index of cluster)
    idxN = idxN/idxC; % cluster index is larger than 1, normalize it to 1
    cut(idxC) = idxN'*L*idxN;
    den(idxC) = idxN'*A*idxN;
    sum_ratio(idxC) = cut(idxC)/den(idxC);
end
value = sum(sum_ratio);

end

