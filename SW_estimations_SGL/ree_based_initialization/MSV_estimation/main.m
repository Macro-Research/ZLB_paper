   clear;clc;
   addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
   lastwarn('Success');

names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','mu_p','mu_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r1','eta_r2', 'eta_p', 'eta_w',...
           'gain','p_11','p_22','rbar_zlb'} ] ;     


load('param_init.mat');

%init_H4=nhess_diagonal(@likelihood,param_init');init_H4=inv(init_H4);
% init_H4=hessian_min_1S(@likelihood,param_init',10e-5);init_H4=inv(init_H4);
% save init_H4.mat init_H4;
load('init_H4.mat')
options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',20000);
[x fh]=fminsearch('likelihood',param_init,options);
%[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param_init',init_H4,[] ,10^(-5),9999);% laplace_=laplace_approximator(fh,x,H);
%laplace1_=laplace_approximator(fh,x,H);


