function dx = CG(K,Iden,h,dx)
% dx is the initial setting of the iteration
% imax = 6*sqrt(size(K,1));
% Ak = sparse(K'*K + epi*speye(size(K,1),size(K,1)));
imax = size(K,1);
Kdx = K*dx;
r = h - K'*Kdx-Iden*dx;
d = r;
deltaNew = r'*r;
delta0 = deltaNew;
epi = 1e-3;
i = 0;

while i <= imax && deltaNew > epi^2*delta0
    dtemp = K*d;
    q = K'*dtemp+Iden*d;
    alpha = deltaNew/(d'*q);
    dx = dx + alpha * d;
    r = r - alpha * q;
    deltaOld = deltaNew;
    deltaNew = r'*r;
    beta = deltaNew/deltaOld;
    d = r + beta * d;
    i = i + 1;
end
end

