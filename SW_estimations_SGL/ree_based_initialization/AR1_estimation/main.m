   clear;clc;
   addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
   lastwarn('Success');
%rng(5) works in the baseline case with fixed beliefs 0.75
% % seed=round(1000*rand);
%  seed=417;
% rng(seed);
names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','mu_p','mu_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r1','eta_r2', 'eta_p', 'eta_w',...
           'gain','p_11','p_22','rbar_zlb'} ] ;     

% param_init=[4.2,1.1,0.75,0.7,2.36,0.7,0.5,0.5,0.46,1.5,...
%     1.43,0.89,0.07,0.05,...
%     0.57,0.18,1.81,0.42,0.13,...
%     0.97,0.36,0.99,0.54,0.35,0.65,0.55,0.5,0.5,0.5,...
%     0.41,0.55,0.4,1.24,0.11,0.05,0.19,0.85,...
%     0.03,0.9,0.7,0.05];
% [x,fh,exitflag,output,grad,H]=fminunc('likelihood',param_init);
% % init_H=eye(41);
 %load('init_H.mat');
  load('param_init.mat');
%  init_H2=.1*init_H2;
 % fmincon('likelihood',param_init);
 %init_H=0.9*init_H;
% init_H=0.1*eye(41);
%init_H=5*eye(41);
% init_H1=nhess_diagonal(@likelihood,param_init);
% init_H1=inv(init_H1);
% save init_H1.mat init_H1;

% init_H1=hessiancsd(@likelihood,param_init);
% init_H1=inv(diag(diag(init_H1)));
% save init_H1.mat init_H1;
 %load('init_H1.mat');
% init_H1=eye(41);
%init_H1=nhess_diagonal(@likelihood,param_init);
% init_H1=nearestSPD(inv(diag(diag(init_H1))));

 %load('init_H1.mat');
 %init_H_fullSample3=nhess_diagonal(@likelihood,param_init);
 %init_H_fullSample3=inv(init_H_fullSample3);init_H_fullSample3=diag(diag(init_H_fullSample3));
 %save init_H_fullSample3.mat init_H_fullSample3;
 load('init_H_fullSample3.mat')
[fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param_init,init_H_fullSample3,[] ,10^(-4),9999);% laplace_=laplace_approximator(fh,x,H);
  % options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',20000);

   %objective=@(x) likelihood(x)
   %options=optimset('Display','iter');
   %  [x fh]=fminsearch('likelihood',param_init,options);
%      [x fh]=particleswarm(objective,length(param_init));
%      x=patternsearch(objective,param_init,[],[],[],[],[],[],[],options);
%   H=nhess_diagonal(@likelihood,x);H=inv(H);
%H=diag(diag(H));
   laplace1_=laplace_approximator(fh,x,H);
 
%nhess(@likelihood,x);
% [fh,x,gh,H,itct,fcount,retcodeh] = csminwel('likelihood',param',init_H,[] ,10^(-5),1000);
% disp(x);
% results=table(names',x)
%  laplace1=laplace_approximator(fh,x,H)
% Sigma=nhess(@likelihood,x);
% Sigma=inv(Sigma);
% laplace2=laplace_approximator(fh,x,Sigma)
% %smoothed_variables(x);
% smoother;



%other solvers
%   [x fh]=particleswarm(objective,length(param_init));
%   x=patternsearch(objective,param_init,[],[],[],[],[],[],[],options);
