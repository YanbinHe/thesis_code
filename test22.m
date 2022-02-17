%% test of status-changing users using method 2-especial for bipartite graph
% clc
% clear
% load test_for_user_move8_containsAll.mat
% load test_for_user_move5.mat
% load test_for_user_Rotate1.mat
% % load test_for_user_move_new1.mat
load newnew.mat
Bini = B_eff;
%% update
numSimu = 1200;
numC = 3;
kkk = 5; % tracking the first kkk eigenvectors
epi = 1e-3;
userUpEach = 20;
epiIden = epi*speye(37);

time = zeros(2,numSimu/userUpEach);
Dx = diag(sum(Bini,2));
Dy = diag(sum(Bini,1));
Btilde = Dx^-0.5*Bini*Dy^-0.5;
[u,s,v] = svd(Btilde);

bIni = Dx^-0.5*u(:,2:kkk); % contains 2-th to kkk-th eigenvectors
muIni = diag(s(2:kkk,2:kkk)).^2;% contains 2-th to kkk-th eigenvalues

%% give initial graph partitions
W = [zeros(size(B_eff,1),size(B_eff,1)),B_eff;
    B_eff',zeros(size(B_eff,2),size(B_eff,2))];
label4Ini = kmeans(bIni(:,1:2),numC,'Replicates',10);% only cluster base station
label4_extIni = reduced_cluster(B_eff,label4Ini);
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
%     if label4_extIni(i) == 1
%     scatter(user_coord(1,i),user_coord(2,i),'g.');
%     hold on
%     end
% if label4_extIni(i) == 2
%     scatter(user_coord(1,i),user_coord(2,i),'r.');
%     hold on
% end
%     if label4_extIni(i) == 3
%     scatter(user_coord(1,i),user_coord(2,i),'b.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');

%%
bIter = bIni;
muIter = muIni;
% initial status

for ii = 1:numSimu/userUpEach
ii
newCol = B_effT(:,(ii-1)*userUpEach+1:ii*userUpEach);
Bold = Bini;
Bnew = Bini;
Bnew(:,(ii-1)*userUpEach+1:ii*userUpEach) = newCol;

% dB = Bnew - Bold;
Dx = sparse(diag(sum(Bold,2)));
Dy = sparse(diag(sum(Bold,1)));

BtildeOld = (Bold*Dy^-1*Bold');

Bnew = sparse(Bnew);

BtildeOld = (Bold*Dy^-1*Bold');

Bnew = sparse(Bnew);

% ggg = @() updateFunc(kkk,bIter,muIter,BtildeOld,Dx,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg);

Dxnew = sparse(diag(sum(Bnew,2)));
Dynew = sparse(diag(sum(Bnew,1)));

% tic
Dynew = sparse(diag(sum(Bnew,1)));
Dxnew = sparse(diag(sum(Bnew,2)));
dDx = (Dxnew - Dx);
dDy = Dynew - Dy;
BtildeNew = full(Bnew*Dynew^-1*Bnew');

dBt = (BtildeNew - BtildeOld);
% update


% ggg = @() updateFunc(kkk,dmu,db,bIter,dBt,muIter,dDx,BtildeNew,Dxnew,Dx,epi,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg);

tic;
db = zeros(size(bIter));
dmu = zeros(kkk-1,1);
dbtemp = zeros(37,1);

% dDy = Dynew - Dy;

BtildeNew = (Bnew*Dynew^-1*Bnew');

dBt = (BtildeNew - BtildeOld);
% update
db(37,4) = 0;
% zeros(size(bIter));
dmu(4,1) = 0;
% zeros(kkk-1,1);
% dbtemp = zeros(37,1);
% tic;

for EigenIdx = 1:kkk-1
    for iii = 1:2
    dmu(EigenIdx) = (((bIter(:,EigenIdx))'*(dBt-muIter(EigenIdx)*dDx))*((bIter(:,EigenIdx))+db(:,EigenIdx)))/(((bIter(:,EigenIdx))'*Dxnew)*((bIter(:,EigenIdx))+db(:,EigenIdx)));


    K = sparse(BtildeNew - ((muIter(EigenIdx) + dmu(EigenIdx))*Dxnew));
    h = sparse(dmu(EigenIdx)*Dx+(muIter(EigenIdx)+dmu(EigenIdx))*dDx-dBt)*(bIter(:,EigenIdx));
%     Atemp = (K'*K + epiIden);

    K = (BtildeNew - ((muIter(EigenIdx) + dmu(EigenIdx))*Dxnew));
    h = (dmu(EigenIdx)*Dx+(muIter(EigenIdx)+dmu(EigenIdx))*dDx-dBt)*(bIter(:,EigenIdx));
    Atemp = (K'*K + (epiIden));
    Btemp = (K'*h);
    db(:,EigenIdx) = CG(K,epiIden,Btemp,db(:,EigenIdx));
    end
end

bIter = bIter + db;
muIter = muIter + dmu;

[muIter,order] = sort(muIter,'descend');
bIter = bIter(:,order);

label3 = kmeans(bIter(:,1:2),numC,'Replicates',10);% only cluster base station
label3_ext = reduced_cluster(Bnew,label3);

time(1,ii) = toc;
% ggg = @() updateFunc(kkk,dmu,db,bIter,dBt,muIter,dDx,BtildeNew,BtildeOld,Dxnew,Dx,Dy,epi,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg);

Dxx = diag(sum(Bnew,2));
Dyy = diag(sum(Bnew,1));
Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);



Bnew = sparse(Bnew);
tic;
Dxx = sparse(diag(sum(Bnew,2)));
Dyy = sparse(diag(sum(Bnew,1)));
Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);

[cc,ss,vv] = svd(Btildee);
bb = Dxx^-0.5*cc(:,2:kkk);
uu = Dyy^-0.5*vv(:,2:kkk);
muu = diag(ss(2:kkk,2:kkk));
label2 = kmeans([[bb(:,1);uu(:,1)],[bb(:,2);uu(:,2)]],numC,'Replicates',10);% base line
time(2,ii) = toc;

BtildeNew=full(BtildeNew);
Dxnew=full(Dxnew);

Bnew = sparse(Bnew);


% tic
Dxx = sparse(diag(sum(Bnew,2)));
Dyy = sparse(diag(sum(Bnew,1)));
% Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);
Btildee = full(Dxx^-0.5*Bnew*Dyy^-1*Bnew'*Dxx^-0.5);

% hhh = @() exactFunc(kkk,numC,Bnew);
% time(3,ii) = timeit(hhh);

Dxx = sparse(diag(sum(Bnew,2)));
Dyy = sparse(diag(sum(Bnew,1)));


[ccc,eee] = svd(Btildee);
% [eee,orderExtEv] = sort(diag(eee),'descend');
% ccc = ccc(:,orderExtEv);
bbb = Dxx^-0.5*ccc(:,2:kkk);
label5 = kmeans(bbb(:,1:2),numC,'Replicates',10);% only cluster base station
label5_ext = reduced_cluster(Bnew,label5);
% time(3,ii) = toc;
% hhh = @() exactFunc(Btildee,kkk,numC,Bnew,Dxx);
% time(3,ii) = timeit(hhh);
% re-order the eigenvectors and eigenvalues

Btildee = full(Dxx^-0.5*Bnew*Dyy^-1*Bnew'*Dxx^-0.5);
% Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);

[ccc,~] = svd(Btildee);
bbb = Dxx^-0.5*ccc(:,2:kkk);
label5 = kmeans(bbb(:,1:2),numC,'Replicates',10);% only cluster base station
label5_ext = reduced_cluster(Bnew,label5);


dif(1,ii) = sum(abs(muu(1:2).^2 - muIter(1:2)));

[Proj,~] = svd(bb(:,1:2));
Proj = Proj(:,1:2)*Proj(:,1:2)';

[Basis1,~] = svd(bIter(:,1:2));
Basis1Proj = Basis1(:,1:2)*Basis1(:,1:2)';

[Basis2,~] = svd(bIni(:,1:2));
Basis2Proj = Basis2(:,1:2)*Basis2(:,1:2)';
% 

dif(2,ii) = 1 - 0.5*(abs((bb(:,1)'*bIter(:,1))/(norm(bb(:,1))*norm(bIter(:,1)))) + abs((bb(:,2)'*bIter(:,2))/(norm(bb(:,2))*norm(bIter(:,2)))));
dif(3,ii) = 1 - 0.5*(abs((bb(:,1)'*bIni(:,1))/(norm(bb(:,1))*norm(bIni(:,1))))+abs((bb(:,2)'*bIni(:,2))/(norm(bb(:,2))*norm(bIni(:,2)))));
dif(4,ii) = sum(abs(muu(1:2).^2 - muIni(1:2)));
dif(5,ii) = norm(-Proj/trace(Proj) + Basis1Proj/trace(Basis1Proj),'fro');
dif(6,ii) = norm(-Proj/trace(Proj) + Basis2Proj/trace(Basis2Proj),'fro');
% dif(5,ii) = norm(-bb(:,1:2)*bb(:,1:2)'/trace(bb(:,1:2)'*bb(:,1:2)) + bIter(:,1:2)*bIter(:,1:2)'/trace(bIter(:,1:2)'*bIter(:,1:2)),'fro');
% dif(6,ii) = norm(-bb(:,1:2)*bb(:,1:2)'/trace(bb(:,1:2)'*bb(:,1:2)) + bIni(:,1:2)*bIni(:,1:2)'/trace(bIni(:,1:2)'*bIni(:,1:2)),'fro');
% dif(5,ii) = norm(Basis1Proj*bb(:,1:2),'fro')/norm(bb(:,1:2),'fro');
% dif(6,ii) = norm(Basis2Proj*bb(:,1:2),'fro')/norm(bb(:,1:2),'fro');

% label1 = kmeans([[b1;u1],[b2;u2]],numC,'Replicates',10);% full k means
% label2 = kmeans([[bb(:,1);uu(:,1)],[bb(:,2);uu(:,2)]],numC,'Replicates',10);% base line
% label3 = kmeans(bIter(:,1:2),numC,'Replicates',10);% only cluster base station
% label3_ext = reduced_cluster(Bnew,label3);
label4 = kmeans(bIni(:,1:2),numC,'Replicates',10);% only cluster base station
label4_ext = reduced_cluster(Bnew,label4);
% label5 = kmeans(bb(:,1:2),numC,'Replicates',10);% only cluster base station
% label5_ext = reduced_cluster(Bnew,label5);

W = [zeros(size(Bnew,1),size(Bnew,1)),Bnew;
    Bnew',zeros(size(Bnew,2),size(Bnew,2))];
cost(1,ii) = MinMaxCut_gen([label4;label4_ext],W);
cost(2,ii) = MinMaxCut_gen(label2,W);
cost(3,ii) = MinMaxCut_gen([label3;label3_ext],W);
cost(4,ii) = MinMaxCut_gen([label5;label5_ext],W);
Bini = Bnew;
% iniU = iniU + 1;
end
%% plotting
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
xlabel('#Perturbations');
ylabel('Estimation Error');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');%
figure
plot(cost(1,1:60),'.-b','DisplayName','Without Update+User Assignment+K-means')%m1-lambda
hold on
plot(cost(2,1:60),'-^r','DisplayName','Baseline')%m1-all
hold on
plot(cost(3,1:60),'-xr','DisplayName','Update+User Assignment+K-means')%m1-37
hold on
plot(cost(4,1:60),'-om','DisplayName','Baseline+User Assignment+K-means')
hold on
xlabel('#Perturbations');
ylabel('MinMaxCut value');
grid on
set(gca,'fontname','DejaVuSans');%,'yscale','log'
% label3 = [label3;label3_ext];
figure
plot(time(2,:)','x-b','DisplayName','Time Consumption Baseline')%m1-lambda
hold on
plot(time(1,:)','-or','DisplayName','Time Consumption Using Update')%m1-all
hold on
% plot(time(3,:)','-.r','DisplayName','Time Consumption Baseline UserAssignment')%m1-all
% hold on
xlabel('#Perturbations');
ylabel('Second');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');
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
for i = 1:1200
    if label3_ext(i) == 3
    scatter(userRotate(1,i),userRotate(2,i),'b.');
    hold on
    end
if label3_ext(i) == 2
    scatter(userRotate(1,i),userRotate(2,i),'r.');
    hold on
end
    if label3_ext(i) == 1
    scatter(userRotate(1,i),userRotate(2,i),'g.');
    hold on
    end
end
set(gca,'fontname','DejaVuSans');