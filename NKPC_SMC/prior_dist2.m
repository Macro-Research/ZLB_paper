function [prior]=prior_dist2(param)
prior1=normpdf(param(1),0,0.1);
mu=0;sigma2=0.15^2;b  = sigma2/mu;a  = mu/b;
prior2=gampdf(param(2),a,b);
mu=0;sigma2=0.15^2;b  = sigma2/mu;a  = mu/b;
prior3=gampdf(param(3),a,b);
mu=0.2;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior4= betapdf(param(4),a,b);
mu=3;sigma2=0.5^2;b  = sigma2/mu;a  = mu/b;
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

%prior14=unifpdf(param(14),0,0.5);


prior=prior1*prior2*prior3*prior4*prior5*prior6*prior7*...
    prior8*prior9*prior10*prior11*prior12*prior13;
prior=-log(prior);
end