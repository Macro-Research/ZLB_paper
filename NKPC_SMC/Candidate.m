clear;%clc;close all;
load('init_H.mat');
  param=[.5474    0.6441      0.8557      0.0264     3.0583       1.4879... 
   0.4660     0.3993      0.3255      0.8718     0.7289     0.2993      0.2944  ];
%Sigma=H;
load('full_dataset.mat');
first_obs=200;last_obs=220;dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);
first_obs=200;last_obs=220;dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);

f = @(x) likelihoodKalman(x,dataset);
[fh,x,gh,H,itct,fcount,retcodeh] = csminwel(f,param',init_H,[] ,10^(-4),500);
% Sigma=nhess(@likelihood,mode);
% Sigma=inv(Sigma);
% Sigma=nearestSPD(Sigma);
Sigma=H;mode=x;
save MH_Candidate.mat mode Sigma;