function [prior]=NKPC_prior(param)
prior1=normpdf(param(1),0,0.25);
mu=0.62;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior2=gampdf(param(2),a,b);
mu=1;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior3=gampdf(param(3),a,b);
mu=0.3;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior4= betapdf(param(4),a,b);
mu=2;sigma2=0.5^2;b  = sigma2/mu;a  = mu/b;
prior5=gampdf(param(5),a,b);
mu=1.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior6=gampdf(param(6),a,b);
mu=0.5;sigma2=0.25^2;b  = sigma2/mu;a  = mu/b;
prior7=gampdf(param(7),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior8= betapdf(param(8),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior9= betapdf(param(9),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior10= betapdf(param(10),a,b);
[s,nu] = inverse_gamma_specification(0.1,4, 0,1, false, 'name');
prior11= exp(lpdfig1(param(11),s,nu));
prior12= exp(lpdfig1(param(12),s,nu));
prior13= exp(lpdfig1(param(13),s,nu));


 prior14=normpdf(param(14),0.1,0.25);
% 
% mu=0.03;sigma2=0.005^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior15=  betapdf(param(15),a,b);
  prior15=unifpdf(param(15),0.005,0.05);
% mu=0.02;sigma2=0.01^2;b  = sigma2/mu;a  = mu/b;
% prior15=gampdf(param(15),a,b);

 
 mu=0.1;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
 prior16=betapdf(param(16),a,b);
  mu=0.3;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
 prior17=betapdf(param(17),a,b);
 
 
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