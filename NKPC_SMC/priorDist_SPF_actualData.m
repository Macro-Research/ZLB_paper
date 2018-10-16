function [prior]=priorDist_SPF_actualData(param)

%normal regime parameters
prior1=normpdf(param(1),0.5,0.2);
% mu=0.62;sigma2=0.15^2;b  = sigma2/mu;a  = mu/b;
prior2=normpdf(param(2),0.5,0.2);
% mu=0.5;sigma2=0.15^2;b  = sigma2/mu;a  = mu/b;
prior3=normpdf(param(3),1,0.2);
mu=0.2;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
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
prior11=exp(lnpdfig(param(11),0.7,2));
prior12=exp(lnpdfig(param(12),0.3,2));
prior13=exp(lnpdfig(param(13),0.3,2));

%zlb regime parameters: monetary policy: phi_pinf phi_y rho_r, and 
%shock st. dev: eta_y eta_pinf eta_r in that order

prior14=normpdf(param(14),0,0.2);
prior15=normpdf(param(15),0,0.2);
prior16=normpdf(param(16),0,0.2);
prior17=exp(lnpdfig(param(17),0.7,2));
prior18=exp(lnpdfig(param(18),0.3,2));
prior19=exp(lnpdfig(param(19),0.03,2));
mu=0.1;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior20=betapdf(param(20),a,b);
mu=0.05;sigma2=0.025^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior21=betapdf(param(21),a,b);
prior22=normpdf(param(22),0,0.2);

prior=prior1*prior2*prior3*prior4*prior5*prior6*prior7*...
    prior8*prior9*prior10*prior11*prior12*prior13*...
    prior14*prior15*prior16*prior17*prior18*prior19*prior20*prior21*prior22;
 prior=-log(prior);
end