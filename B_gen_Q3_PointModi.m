% run different simulation methods
% This code will first load 
clc
clear
load testAll7Modi.mat
solIdx = 1; 
Wini = Adjini; 
Lini = diag(Wini*ones(size(Wini,1),1)) - Wini;
newStatus = newStatusini;
UserStatusChanging;
%%
clc
clear
load testAll7Modi.mat
solIdx = 1; 
Wini = Adjini;
Lini = diag(Wini*ones(size(Wini,1),1)) - Wini;
newStatus = newStatusini_status;
NewUsers;