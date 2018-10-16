clear;%clc;close all;
init_H=rand*diag([0 0 0 0.1 1 0.5 0.25 0.25 0.25 0.25 1 1 1 0.5 0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1]);
 param=(1-0*0.1)*[0 0.5 1 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3...
     0 0 0 0.7 0.3 0.03 0.1 0.05 0];

%Sigma=H;
% load('simulated_dataset');
% first_obs=1;last_obs=100;
% dataset=dataset(first_obs:last_obs,:);

load('simulated_dataset');
dataset=[gap_hp,pinfobs,robs];
dataset=dataset(44:end,:);

f = @(x) likelihoodSPF(x,dataset);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',9999);
 [x,fval,exitflag,output]=fminsearch(f,param',options);
 %[fh,x,gh,H,itct,fcount,retcodeh] = csminwel(f,param',init_H,[] ,10^(-4),500);
% Sigma=nhess(@likelihood,mode);
% Sigma=inv(Sigma);
% Sigma=nearestSPD(Sigma);
% Sigma=H;mode=x;
% save MH_Candidate.mat mode Sigma;