clear;clc;close all;
addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
load('init_H7.mat');
load('param_init.mat');
seed=round(1000*rand);
% rng(800);
%param_init=[-0.2,0.4,.8,0.02,1.71,1.53,0.35,0.91,0.8,0.88,0.14,0.04,0.31,0.03,0.01,0.02,0.13,0.001];

objective= @(param) likelihood(param);

[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param_init',init_H7,[] ,10^(-4),9999);
%options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',500);
%x=fminsearch('likelihood',param_init,options);

laplace_=laplace_approximator(fh,x,H)
save('estimation_results.mat');

%   Sigma=nhess(@likelihood,x);
%   Sigma=inv(Sigma);
% laplace2_=laplace_approximator(fh,x,Sigma);