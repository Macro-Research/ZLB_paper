clear;clc;close all;
addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
addpath('c:\users\tolga\desktop\zlb_paper\distributions');
load('init_H2.mat');
seed=round(1000*rand);
rng(800);
%  param_init=[0.04 0.8 1.2 0.027 2.7 1.26 0.46...
%      0.69 0.71 0.84 0.73 0.28 0.31 0.001 0.004 0.02 0.10 0.001];    
param_init=[0.06 0.97 1.34 0.035 3.13 1.32 0.35 0.38 0.05 0.97 0.75 0.26 0.33 0.04 0.02 0.02 0.1 0.04];

%param_init=[0.24 0.65 1 0.0073 4.27 1.39 0.46 0.87 0.87 0.79 0.16 0.03 0.30 0.05 0.01 0.02 0.11 0.01];    
%parameters: [bar_y bar_pi bar_r1 kappa tau phi_pi phi_y rho_y rho_pi rho_r 
% eta_y eta_pi eta_r1 bar_r2 eta_r2 1-p_11 1-p_22 gain]
% init_H=eye(17);

objective= @(param) likelihood(param);

[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param_init',init_H2,[] ,10^(-4),9999);
%options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',999);
%x=fminsearch('likelihood',param_init,options);

laplace_=laplace_approximator(fh,x,H)

save('estimation_results.mat');