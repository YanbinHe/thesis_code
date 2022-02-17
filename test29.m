clc
clear
load test_for_user_new_cluster.mat

Bini = B_eff;
%% warmStart-versus 21 only valid-now it's not warm start, just adding columns vs 22
numSimu = 400;
numC = 3;
kkk = 5;
epi = 1e-3;
userUpEach = 5;
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
%%
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
%     scatter(user_coord(1,i),user_coord(2,i),'b.');
%     hold on
% end
%     if label4_extIni(i) == 3
%     scatter(user_coord(1,i),user_coord(2,i),'r.');
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
Bold = [Bini, zeros(size(Bini,1),userUpEach)];
Bnew = [Bini, newCol];


% dB = Bnew - Bold;
% Dx keeps constant
Dx = (diag(sum(Bold,2)));
Dy = diag(sum(Bnew,1));
% Dy = diag(sum(Bold,1));


Dynew = diag(sum(Bnew,1));

BtildeOld = (Bold*Dy^-1*Bold');

Bnew = sparse(Bnew);


Dxnew = (diag(sum(Bnew,2)));
Dynew = sparse(diag(sum(Bnew,1)));

BtildeOld = (Bold*Dy^-1*Bold');

dDx = (diag(sum(newCol,2)));% dDy = zeros(size(Dy,1),size(Dy,2));
dDy = Dynew - Dy;

% tic;
BtildeNew = (Bnew*(Dy+dDy)^-1*Bnew');
Dxnew = (diag(sum(Bnew,2)));

BtildeNew = full(Bnew*(Dynew)^-1*Bnew');
% 
dDx = (diag(sum(newCol,2)));
Bnew = sparse(Bnew);
% ggg = @() updateFunc(kkk,dmu,db,bIter,dBt,muIter,dDx,BtildeNew,BtildeOld,Dxnew,Dx,Dy,epi,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg);


dBt =(BtildeNew - BtildeOld);
% update

db = (zeros(size(bIter)));
dmu = zeros(kkk-1,1);


tic;

% tic;

for EigenIdx = 1:kkk-1
    for iii = 1:2
    dmu(EigenIdx) = (bIter(:,EigenIdx)'*(dBt-muIter(EigenIdx)*dDx)*(bIter(:,EigenIdx)+db(:,EigenIdx)))/(bIter(:,EigenIdx)'*(Dx + dDx)*(bIter(:,EigenIdx)+db(:,EigenIdx)));


    K = sparse(BtildeNew - (muIter(EigenIdx) + dmu(EigenIdx))*(Dxnew));
    h = sparse((dmu(EigenIdx) *Dx+muIter(EigenIdx)*dDx+dmu(EigenIdx)*dDx-dBt)*bIter(:,EigenIdx));

    K = (BtildeNew - (muIter(EigenIdx) + dmu(EigenIdx))*(Dxnew));
    h = ((dmu(EigenIdx)*Dx+muIter(EigenIdx)*dDx+dmu(EigenIdx)*dDx-dBt)*bIter(:,EigenIdx));

    Initial(37,1) = 0;
    db(:,EigenIdx) = CG(K,epiIden,(K'*h),Initial);
    end
end

bIter = bIter + db;
muIter = muIter + dmu;

[muIter,order] = sort(muIter,'descend');
bIter = bIter(:,order);

label3 = kmeans(bIter(:,1:2),numC,'Replicates',10);% only cluster base station
label3_ext = reduced_cluster(Bnew,label3);


time(1,ii) = toc;
% ggg = @() updateFunc(kkk,dmu,db,bIter,dBt,muIter,dDx,BtildeNew,Dxnew,Dx,epi,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg,12);

Dxx = diag(sum(Bnew,2));
Dyy = diag(sum(Bnew,1));
Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);


% ggg = @() updateFunc(kkk,dmu,db,bIter,dBt,muIter,dDx,BtildeNew,BtildeOld,Dxnew,Dx,Dy,epi,numC,Bnew,epiIden);
% time(1,ii) = timeit(ggg,12);
% time(1,ii) = toc;

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

% tic;
Dxx = diag(sum(Bnew,2));
Dyy = diag(sum(Bnew,1));
Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);

Bnew = sparse(Bnew);

Dxx = full(diag(sum(Bnew,2)));
Dyy = sparse(diag(sum(Bnew,1)));
Btildee = full(Bnew*Dyy^-1*Bnew');
% Btildee = full(Dxx^-0.5*Bnew*Dyy^-0.5);

% [cc,eigvvv] = eig(BtildeNew,Dxnew);
% [eigvvv,orderEig] = sort(eigvvv,'descend');
% ccc = ccc(:,orderEig);
% bbb = Dxx^-0.5*ccc(:,2:kkk);

% [ccc,~] = svd(Btildee);
% bbb = Dxx^-0.5*ccc(:,2:kkk);
% label5 = kmeans(bb(:,1:2),numC,'Replicates',10);% only cluster base station
% label5_ext = reduced_cluster(Bnew,label5);
% tic;
[ccc,eee] = eig(Btildee,Dxx);
[eee,orderExtEv] = sort(diag(eee),'descend');
ccc = ccc(:,orderExtEv);
bbb = ccc(:,2:kkk);%Dxx^-0.5*
label5 = kmeans(bbb(:,1:2),numC,'Replicates',10);% only cluster base station
label5_ext = reduced_cluster(Bnew,label5);

% time(3,ii) = toc;

% hhh = @() exactFunc(Btildee,kkk,numC,Bnew,Dxx);
% time(3,ii) = timeit(hhh);
% re-order the eigenvectors and eigenvalues

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
set(gca,'fontname','DejaVuSans','yscale','log');
figure
plot(cost(1,:),'.-b','DisplayName','Without Update+User Assignment+K-means')%m1-lambda
hold on
plot(cost(2,:),'-^r','DisplayName','Baseline')%m1-all
hold on
plot(cost(3,:),'-xr','DisplayName','Update+User Assignment+K-means')%m1-37
hold on
plot(cost(4,:),'-om','DisplayName','Baseline+User Assignment+K-means')
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
plot(time(3,:)','-.r','DisplayName','Time Consumption Baseline UserAssignment')%m1-all
hold on
xlabel('#Perturbations');
ylabel('Second');
grid on
set(gca,'fontname','DejaVuSans','yscale','log');%
%%
% UserCoordFinalPlot = [user_coord, userCoordTest];
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1600
%     if label3_ext(i) == 3
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'b.');
%     hold on
%     end
% if label3_ext(i) == 2
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'b.');
%     hold on
% end
%     if label3_ext(i) == 1
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'b.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');
%%
label6 = kmeans(bIter(:,1:3),numC+1,'Replicates',10);% only cluster base station
label6_ext = reduced_cluster(Bnew,label6);
%%
% figure
% for i = 1:num_bs_g
%     X(:,i) = distance0(1,i) + unit*cos(t);
%     Y(:,i) = distance0(2,i) + unit*sin(t);
% end
% for i = 1:num_bs_g
%     plot(X(:,i),Y(:,i),'k')
%     hold on
% end
% for i = 1:1600
%     if label6_ext(i) == 3
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'b.');
%     hold on
%     end
% if label6_ext(i) == 2
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'r.');
%     hold on
% end
%     if label6_ext(i) == 1
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'k.');
%     hold on
%     end
%     if label6_ext(i) == 4
%     scatter(UserCoordFinalPlot(1,i),UserCoordFinalPlot(2,i),'g.');
%     hold on
%     end
% end
% set(gca,'fontname','DejaVuSans');
%%
figure
EigneValueIdx = linspace(2,5,4);
plot(EigneValueIdx,muIni,'-vg','DisplayName','Eigenvalue without update');
hold on
plot(EigneValueIdx,muIter,'-bx','DisplayName','Eigenvalue with update');
hold on
plot(EigneValueIdx,muu.^2,'-or','DisplayName','Eigenvalue by EVD');
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
set(gca,'fontname','DejaVuSans');
xlabel('Index [k]');