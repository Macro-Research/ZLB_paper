 function [prior]=NKPC_prior(param)
%clear;clc;close all;
%param=[0 0 0 0.03 2 1.5 0.5 0.5 0.5 0.9 0.3 0.3 0.3 0 0.01 0.01 0.1 0.01];  

prior1=normpdf(param(1),0,0.25);
prior2=normpdf(param(2),0,0.25);
prior3=normpdf(param(3),0,0.25);
% mu=0.03;sigma2=0.015^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior4= betapdf(param(4),a,b);
prior4=unifpdf(param(4),0,0.5);
mu=2;sigma2=0.5^2;b  = sigma2/mu;a  = mu/b;
prior5=gampdf(param(5),a,b);
mu=1.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior6=gampdf(param(6),a,b);
mu=0.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior7=gampdf(param(7),a,b);
% mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior8= betapdf(param(8),a,b);
% mu=0.5;sigma2=0.25^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior9= betapdf(param(9),a,b);
% mu=0.9;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior10= betapdf(param(10),a,b);
% [s,nu] = inverse_gamma_specification(0.1,4, 0,1, false, 'name');
% prior11= exp(lpdfig1(param(11),s,nu));
% prior12= exp(lpdfig1(param(12),s,nu));
% prior13= exp(lpdfig1(param(13),s,nu));
prior8=unifpdf(param(8),0.01,0.99);
prior9=unifpdf(param(9),0.01,0.99);
prior10=unifpdf(param(10),0.01,0.99);
 prior11= unifpdf(param(11),0.001,1);
 prior12= unifpdf(param(12),0.001,1);
 prior13= unifpdf(param(13),0.001,1);

 prior14=normpdf(param(14),0,0.25);
% 
% mu=0.05;sigma2=0.01^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior15=  betapdf(param(15),a,b);
  prior15=unifpdf(param(15),0.001,0.05);
% mu=0.02;sigma2=0.01^2;b  = sigma2/mu;a  = mu/b;
% prior15=gampdf(param(15),a,b);

 
%  mu=0.1;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
%  prior16=betapdf(param(16),a,b);
%   mu=0.3;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
%  prior17=betapdf(param(17),a,b);

prior16=unifpdf(param(16),0,0.5);
prior17=unifpdf(param(17),0,0.5);
 
 
 mu=0.035;sigma2=0.015^2;b  = sigma2/mu;a  = mu/b;
prior18=gampdf(param(18),a,b);
% prior18=unifpdf(param(18),0.0001,0.1);
%    mu=0.05;sigma2=0.01^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
%   prior18=betapdf(param(18),a,b);
%  prior19=betapdf(param(19),a,b);
%prior18=unifpdf(param(18),0,0.05);

prior=prior1*prior2*prior3*prior4*prior5*prior6*prior7*...
    prior8*prior9*prior10*prior11*prior12*prior13*prior14*prior15*...
    prior16*prior17*prior18;%*prior19;
prior=-log(prior);
 end