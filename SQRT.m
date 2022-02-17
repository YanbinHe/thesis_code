function [sqrtd] = SQRT(D)
% D should be a positive semidefinite diagonal matrix
d = diag(D);

sqrtd = D;

posi = d > 1e-6;
sqrtd(posi,posi) = diag(d(posi).^0.5);
end

