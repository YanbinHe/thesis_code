% clc
% clear
% load 'testGraphQ37.mat'
% solIdx = 1;
%% test 35 dimension varying update
% input description:
% W is the initial adjacency matrix
% L is the initial laplacian matrix
% solIdx: expect a numer 1, 2, 3. x->method x
% output: we can plot graph here in this function, so currently we don't
% have return anything

% what I think, we don't want 2 here since it cannot provide good results,
% based on the idea that we can give a so called static small size graph
% validation section.

% new status should be a fat matrix \in R^{#base stations \times #users}
% we keep adding users to the cellular network using users in the newStatus matrix
%% definition of some basic values
epi = 1e-4;
numClu = 4;% #tracking components
numC = 3;% the number of clusters
numN = 37;% the number of base stations
num_nodes = size(Wini, 1);

b = zeros(1,size(Wini,1));
b([1:37]) = 1;
Db = diag(b);
u = 1 - b;
Du = diag(u);
Wtilde = Db*Wini*Du;    
Wbar = 2*Wtilde;
q = (Wtilde + Wtilde')*ones(num_nodes,1);

userUpEach = 10;
numSimu = 300;
time = zeros(2,numSimu/userUpEach);
%% computing different W-related and L-related matrix
if solIdx == 1
    % filled with code generating Ws and Lbar
    W = (Wbar + Wbar');
    L = 2*Lini;
elseif solIdx == 3
    % filled with code generating Wds and Lbar
    proj_q = eye(num_nodes) - 1/(q'*q)*(q*q');
    [Upen,~,~] = svd(proj_q);% rank n-1
    U_proj = Upen(:,1:end-1);
    U_proji = U_proj;
    
    W = (U_proj'*Wtilde*U_proj)'+(U_proj'*Wtilde*U_proj);
    L = 2*U_proj'*Lini*U_proj;
end

num_tol = 1e-6;
[Ub, Eb] = eig(L);
Ub2 = Ub*pinv(SQRT(Eb),num_tol);
[Ua,Ea] = eig(0.5*((Ub2'*W*Ub2).'+(Ub2'*W*Ub2)));
V = Ub2*Ua;
Veff = V(:,abs(diag(Ea)) > num_tol);

W = sparse(W);
L = sparse(L);
% large <- small
GenEigV = Veff(:,end-numClu + 1:end); % tracked components
EigV = diag(Ea(end-numClu + 1:end,end-numClu + 1:end)); % tracked eigenvalues
%% give initial graph partitions
% labelIni = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);
labelIni = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1200
%     if labelIni(i+37) == 1
%     scatter(user_coord(1,i),user_coord(2,i),'g.');
%     hold on
%     end
% if labelIni(i+37) == 2
%     scatter(user_coord(1,i),user_coord(2,i),'b.');
%     hold on
% end
%     if labelIni(i+37) == 3
%     scatter(user_coord(1,i),user_coord(2,i),'r.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');
%%
GenEigVi = GenEigV;
EigVi = EigV;
% --------------------initialization ends-------------------------------- %
for simu = 1:numSimu/userUpEach
simu
num_nodes = size(Wini, 1) + userUpEach;
num_users = num_nodes - numN ;

newCol = newStatus(:,(simu-1)*userUpEach+1:simu*userUpEach);% new data
Worinew = [Wini,[newCol;zeros(num_users-userUpEach,userUpEach)];
           [newCol',zeros(userUpEach,num_users-userUpEach)],zeros(userUpEach)];

Db = blkdiag(Db,zeros(userUpEach,userUpEach));
Du = blkdiag(Du,eye(userUpEach,userUpEach));

tic;
Wtildenew = Db*Worinew*Du;    
Wbarnew = 2*Wtildenew;
qnew = (Wtildenew + Wtildenew')*ones(num_nodes,1);
Lorinew = diag(Worinew*ones(num_nodes,1)) - Worinew;

if solIdx == 1
    % filled with code generated Ws and Lbar
    
    Wnew = sparse(Wbarnew + Wbarnew');
    Lnew = sparse(2*Lorinew);
    
elseif solIdx == 3
    % filled with codes generating Wds and Lbar
    %% exact compute
    % proj_q = eye(num_nodes) - 1/(q'*q)*(q*q');
    % [Upen,~,~] = svd(proj_q);% rank n-1
    % U_proj = Upen(:,[1:end-1]);
    %% approximate
    proj_qnew = eye(num_nodes) - 1/(qnew'*qnew)*(qnew*qnew');
    dProj = proj_qnew - proj_q;
    U_projnew = U_proj + q/norm(q)*q'/norm(q)*dProj*U_proj;
    Wnew = (U_projnew'*Wtildenew*U_projnew)'+(U_projnew'*Wtildenew*U_projnew);% approximated Wnew
    Lnew = 2*U_projnew'*Lorinew*U_projnew;% approximated Lnew
end

dW = sparse(Wnew - blkdiag(W,zeros(userUpEach,userUpEach)));
dL = sparse(Lnew - blkdiag(L,zeros(userUpEach,userUpEach)));

% matrix update ends
% estimate lambda
% tic;
dx = sparse(zeros(size(Lnew,2),numClu));
GenEigV = [GenEigV;zeros(userUpEach,numClu)];
GenEigVi = [GenEigVi;zeros(userUpEach,numClu)];
dlambda = zeros(numClu,1);

for kkk = 1:numClu
for iii = 1:2
    
        dlambda(kkk) = ((GenEigV(:,kkk))'*(dW-(EigV(kkk))*dL)*(GenEigV(:,kkk)+dx(:,kkk)))/((GenEigV(:,kkk))'*(sparse(blkdiag(L,zeros(userUpEach,userUpEach))) + dL)*(GenEigV(:,kkk)+dx(:,kkk)));
        
        K = sparse(sparse(blkdiag(W,zeros(userUpEach,userUpEach)))+dW-(EigV(kkk)+ dlambda(kkk))*(sparse(blkdiag(L,zeros(userUpEach,userUpEach)))+dL));
        h = sparse((dlambda(kkk)*sparse(blkdiag(L,zeros(userUpEach,userUpEach)))+(EigV(kkk))*dL+dlambda(kkk)*dL-dW)*(GenEigV(:,kkk)));

%         Ak = sparse(K'*K + epi*speye(size(K,1),size(K,1)));
        qk = sparse(K'*h);
        % conjugate gradient descent
        dx(:,kkk) = CG(K,epi*speye(size(K,1),size(K,1)),qk,sparse(zeros(size(dx(:,1)))));
end
end

for kkk = 1:numClu
        dlambda(kkk) = ((GenEigV(:,kkk))'*(dW-(EigV(kkk))*dL)*(GenEigV(:,kkk)+dx(:,kkk)))/((GenEigV(:,kkk))'*(sparse(blkdiag(L,zeros(userUpEach,userUpEach))) + dL)*(GenEigV(:,kkk)+dx(:,kkk)));
end

GenEigV = GenEigV + dx;
EigV = EigV + dlambda;
[EigV,order] = sort(EigV,'ascend');
GenEigV = GenEigV(:,order);
% if solIdx == 1
    label3 = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);% using approximate
% end
time(1,simu) = toc;
% exact computation
Wnew = full(Wnew);
Lnew = full(Lnew);

tic;
Wtildenew = Db*Worinew*Du;    
Wbarnew = 2*Wtildenew;
qnew = (Wtildenew + Wtildenew')*ones(num_nodes,1);
Lorinew = diag(Worinew*ones(num_nodes,1)) - Worinew;

if solIdx == 1
    % filled with code generated Ws and Lbar
    
    Wnew = (Wbarnew + Wbarnew');
    Lnew = (2*Lorinew);
    
elseif solIdx == 3
    % filled with codes generating Wds and Lbar
    %% exact compute
    % proj_q = eye(num_nodes) - 1/(q'*q)*(q*q');
    % [Upen,~,~] = svd(proj_q);% rank n-1
    % U_proj = Upen(:,[1:end-1]);
    %% approximate
    proj_qnew = eye(num_nodes) - 1/(qnew'*qnew)*(qnew*qnew');
    dProj = proj_qnew - proj_q;
    U_projnew = U_proj + q/norm(q)*q'/norm(q)*dProj*U_proj;
    Wnew = (U_projnew'*Wtildenew*U_projnew)'+(U_projnew'*Wtildenew*U_projnew);% approximated Wnew
    Lnew = 2*U_projnew'*Lorinew*U_projnew;% approximated Lnew
end
% if solIdx == 3
%     
%     proj_qext = eye(num_nodes) - 1/(qnew'*qnew)*(qnew*qnew');
%     [Upenext,~,~] = svd(proj_qext);% rank n-1
%     U_projext = Upenext(:,[1:end-1]);
% 
%     Wnewext = (U_projext'*Wtildenew*U_projext)'+(U_projext'*Wtildenew*U_projext);
%     Lnewext = 2*U_projext'*Lorinew*U_projext;
% 
%     Projext_q = U_projext*U_projext';
% 
%     num_tol = 1e-6;
%     [Ubext, Ebext] = eig(Lnewext);
%     Ub2ext = Ubext*pinv(sqrtm(Ebext),num_tol);
%     [Uaext,Eaext] = eig(0.5*((Ub2ext'*Wnewext*Ub2ext).'+Ub2ext'*Wnewext*Ub2ext));
%     Vext = Ub2ext*Uaext;
%     Veffext = Vext(:,abs(diag(Eaext)) > num_tol);
%     Eaext = diag(Eaext);
% elseif solIdx == 1
    
    num_tol = 1e-6;
    [Ubext, Ebext] = eig(Lnew);
    Ub2ext = Ubext*pinv(SQRT(Ebext),num_tol);
    [Uaext,Eaext] = eig(0.5*((Ub2ext'*Wnew*Ub2ext).'+Ub2ext'*Wnew*Ub2ext));
    Vext = Ub2ext*Uaext;
    Veffext = Vext(:,abs(diag(Eaext)) > num_tol);
    Eaext = diag(Eaext);
% end

% if solIdx == 1
    label2 = kmeans(Veffext(:,end-numC + 1:end),numC,'Replicates',10);% base line
% end
time(2,simu) = toc;

Wnew = sparse(Wnew);
Lnew = sparse(Lnew);

dif(1,simu) = sum(abs(Eaext(end-numC + 1:end) - EigV(end-numC + 1:end)));
% ./Eaext(end-numC + 1:end)
[Projext,~] = svd(Veffext(:,end-numC + 1:end));
Projext = Projext(:,1:numC)*Projext(:,1:numC)';

[Basis1,~] = svd(GenEigV(:,end-numC + 1:end));
Basis1Proj = Basis1(:,1:numC)*Basis1(:,1:numC)';

[Basis2,~] = svd(GenEigVi(:,end-numC + 1:end));
Basis2Proj = Basis2(:,1:numC)*Basis2(:,1:numC)';
% 

dif(2,simu) = 1 - 1/3*(abs((Veffext(:,end + 1 -3)'*GenEigV(:,end + 1 -3))/(norm(Veffext(:,end + 1 -3))*norm(GenEigV(:,end + 1 -3)))) + abs((Veffext(:,end + 1 -2)'*GenEigV(:,end + 1 -2))/(norm(Veffext(:,end + 1 -2))*norm(GenEigV(:,end + 1 -2))))+ abs((Veffext(:,end + 1 -1)'*GenEigV(:,end + 1 -1))/(norm(Veffext(:,end + 1 -1))*norm(GenEigV(:,end + 1 -1)))));
dif(3,simu) = 1 - 1/3*(abs((Veffext(:,end + 1 -3)'*GenEigVi(:,end + 1 -3))/(norm(Veffext(:,end + 1 -3))*norm(GenEigVi(:,end + 1 -3))))+abs((Veffext(:,end + 1 -2)'*GenEigVi(:,end + 1 -2))/(norm(Veffext(:,end + 1 -2))*norm(GenEigVi(:,end + 1 -2))))+abs((Veffext(:,end + 1 -1)'*GenEigVi(:,end + 1 -1))/(norm(Veffext(:,end + 1 -1))*norm(GenEigVi(:,end + 1 -1)))));
dif(4,simu) = sum(abs(Eaext(end-numC + 1:end) - EigVi(end-numC + 1:end)));
% ./Eaext(end-numC + 1:end)
dif(5,simu) = norm(-Projext/trace(Projext) + Basis1Proj/trace(Basis1Proj),'fro');
dif(6,simu) = norm(-Projext/trace(Projext) + Basis2Proj/trace(Basis2Proj),'fro');
if solIdx == 3
    Proj_newq = U_projnew*U_projnew';
    dif(7,simu) = norm(-Proj_newq/trace(Proj_newq) + Projext_q/trace(Projext_q),'fro');% error of q-perpendicular subspace
end


if solIdx == 1
%     label2 = kmeans(Veffext(:,end-numC + 1:end),numC,'Replicates',10);% base line
%     label3 = kmeans(GenEigV(:,end-numC + 1:end),numC,'Replicates',10);% using approximate
    label4 = kmeans(GenEigVi(:,end-numC + 1:end),numC,'Replicates',10);% using initial
    
    % recursive bisection
    %%% Recursive bisection
%     tic;
%     b = diag(Db);
%     u = diag(Du);
%     itr = ceil(log2(numC));% #iterations to get targeted #clusters. e.g., numC = 5, we should do bisection 
%     % log2(5) = 2.3219 -> 3 times to get 8 clusters, and then merge them
%     label = RecursiveBisection(Worinew,b,u);
%     % we can use the bisection results from method 1 shown previously
%     label = label + 1;% contains 1 and 2
%     numCluster = 2; % already finish bisectioon once
%     for i = 1:itr
%         for j = 1:2^i
% 
%         % break the iteration if #cluster has been met
%         if numCluster >= numC
%             break;
%         end
% 
%         LargestClu = mode(label);
%         subdomain_v = (label == LargestClu);
%         position = find(subdomain_v == 1);
%         sub_W = Worinew(position, position);
%         sub_b = b(position);
%         sub_u = u(position);
% 
%         label_new = RecursiveBisection(sub_W,sub_b,sub_u);% only contains 0-1 labels of nodes in subgraph
%         numCluster = numCluster + 1;
%         % organize the labels since label_new only returns 0-1
% 
%         label(position(label_new == 1)) = numCluster;
% 
%         end
%     end
%     time(3,simu) = toc;
    
elseif solIdx == 3
    SpectralEmbext = U_projext*Veffext(:,end-numC + 1:end);
    label2 = kmeans(SpectralEmbext,numC,'Replicates',10);% base line
    SpectralEmbapp = U_projnew*GenEigV(:,end-numC + 1:end); % 
    label3 = kmeans(SpectralEmbapp,numC,'Replicates',10);% using approximate
    SpectralEmbini = U_proji*GenEigVi(:,end-numC + 1:end);
    label4 = kmeans(SpectralEmbini,numC,'Replicates',10);% using initial
end

G = gsp_graph(Worinew);
cost(1,simu) = MinMaxCut(label2,G,Db,Du);
cost(2,simu) = MinMaxCut(label3,G,Db,Du);
cost(3,simu) = MinMaxCut(label4,G,Db,Du);
% cost(4,simu) = MinMaxCut(label,G,Db,Du);
W = Wnew;
L = Lnew;
Wini = Worinew;
if solIdx == 3
    U_proj = U_projnew;
    proj_q = proj_qnew ;
    q = qnew;
end
end
%%
figure
plot(dif(1,:),'-.r','DisplayName','AE with update')
hold on
plot(dif(2,:),'-^r','DisplayName','DC with update')
hold on
plot(dif(5,:),'-xr','DisplayName','SD with update')
hold on
plot(dif(4,:),'-.b','DisplayName','AE without update')
hold on
plot(dif(3,:),'-Vb','DisplayName','DC without update')
hold on
plot(dif(6,:),'o-b','DisplayName','SD without update')
hold on
if solIdx == 3
plot(dif(7,:),'*-r','DisplayName','SD with (\textbf(U)_q) update')
hold on
end
xlabel('#Perturbations');
ylabel('Estimation Error');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');%
figure
plot(cost(1,:),'-^r','DisplayName','EMP')
hold on
plot(cost(2,:),'-om','DisplayName','UMP')
hold on
plot(cost(3,:),'.-b','DisplayName','WMP')
hold on
% plot(cost(4,:),'-vb','DisplayName','Recursive Bisection')
% hold on
xlabel('#Perturbations');
ylabel('Modified-MinMaxCut value');
grid on
set(gca,'fontname','DejaVuSans');%,'yscale','log'
figure
plot(time(1,:),'-.','DisplayName','Time Consumption UMP')
hold on
plot(time(2,:),'-x','DisplayName','Time Consumption EMP')
hold on
% plot(time(3,:),'-o','DisplayName','Time Consumption Recursive Bisection')
% hold on
xlabel('#Perturbations');
ylabel('Second');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');
%%
coordinate = [user_coord,Newuser_coord_status];
figure
for i = 1:num_bs_g
    X(:,i) = distance0(1,i) + unit*cos(t);
    Y(:,i) = distance0(2,i) + unit*sin(t);
end
for i = 1:num_bs_g
    plot(X(:,i),Y(:,i),'k')
    hold on
end
for i = 1:1500
    if label3(i+37) == 1
    scatter(coordinate(1,i),coordinate(2,i),'b.');
    hold on
    end
if label3(i+37) == 2
    scatter(coordinate(1,i),coordinate(2,i),'g.');
    hold on
end
    if label3(i+37) == 3
    scatter(coordinate(1,i),coordinate(2,i),'r.');
    hold on
    end
end
set(gca,'fontname','DejaVuSans');
%%
figure
for i = 1:num_bs_g
    X(:,i) = distance0(1,i) + unit*cos(t);
    Y(:,i) = distance0(2,i) + unit*sin(t);
end
for i = 1:num_bs_g
    plot(X(:,i),Y(:,i),'k')
    hold on
end
for i = 1:1500
    if label(i+37) == 1
    scatter(coordinate(1,i),coordinate(2,i),'r.');
    hold on
    end
if label(i+37) == 2
    scatter(coordinate(1,i),coordinate(2,i),'b.');
    hold on
end
    if label(i+37) == 3
    scatter(coordinate(1,i),coordinate(2,i),'g.');
    hold on
    end
end
set(gca,'fontname','DejaVuSans');
% %% 
% posi1 = label == 1;
% posi2 = label == 3;
% label(posi1) = 3;
% label(posi2) = 1;
% label ~= label3;