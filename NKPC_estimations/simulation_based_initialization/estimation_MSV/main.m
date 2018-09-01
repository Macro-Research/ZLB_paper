clear;clc;close all;
addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
load('init_H8.mat');
 load('param_init.mat');
seed=round(1000*rand);
%param_init=[-0.3 0.64 0.86 0.004 2.6 1.55 0.43 0.9 0.87 0.89 0.15 0.04 0.32 0.03 0.01 0.02 0.13 0.0173];    

% rng(800);
%  param_init=[0.04 0.8 1.2 0.027 2.7 1.26 0.46...
%      0.69 0.71 0.84 0.73 0.28 0.31 0.001 0.004 0.02 0.10 0.001];    
%param_init=[0.06 0.97 1.34 0.035 3.13 1.32 0.35 0.38 0.05 0.97 0.75 0.26 0.33 0.04 0.02 0.02 0.1 0.04];

%param_init=[0 0.44 0.84 0.02 2.5 1.5 0.32 0.9 0.9 0.85 0.16 0.03 0.30 0.03 0.01 0.01 0.11 0.0013];    
%parameters: [bar_y bar_pi bar_r1 kappa tau phi_pi phi_y rho_y rho_pi rho_r 
% eta_y eta_pi eta_r1 bar_r2 eta_r2 1-p_11 1-p_22 gain]
% init_H=eye(17);
 %options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',999);
 %x=fminsearch('likelihood',param_init,options);


% objective= @(param) likelihood(param);
% H=nhess(@likelihood,param_init');

%[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param_init',init_H8,[] ,10^(-4),9999);
options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',500);
x=fminsearch('likelihood',param_init',options);
% 
 %laplace_=laplace_approximator(fh,x,H)

 save('estimation_results.mat');
 %Sigma=nhess(@likelihood,x);
 %Sigma=inv(Sigma);
% laplace2_=laplace_approximator(fh,x,Sigma);