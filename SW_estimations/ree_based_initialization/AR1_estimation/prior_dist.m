function [prior]=prior_dist(param)
% param=param_init;
% names={'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
%     'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
%            'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
%            'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','mu_p','mu_w','rho_ga',...
%            'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r1','eta_r2', 'eta_p', 'eta_w',...
%            'gain','p_11','p_22','rbar_zlb'}  ;    

prior1=normpdf(param(1),4,1.5);%phi
prior2=normpdf(param(2),1.5,0.375);%sigma_c
% prior2=normpdf(param(2),1.5,0.15);%sigma_c

mu=0.7;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior3= betapdf(param(3),a,b);%lambda
% mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior3= betapdf(param(3),a,b);%lambda

mu=0.5;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior4= betapdf(param(4),a,b);%xi_w
% mu=0.75;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior4= betapdf(param(4),a,b);%xi_w

prior5=normpdf(param(5),2,0.75);%sigma_l
% 
mu=0.5;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior6= betapdf(param(6),a,b);%xi_p
% mu=0.75;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior6= betapdf(param(6),a,b);%xi_p

mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior7= betapdf(param(7),a,b);%iota_w

mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior8= betapdf(param(8),a,b);%iota_p

%prior7=1;prior8=1;

mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior9= betapdf(param(9),a,b);%psi

prior10=normpdf(param(10),1.25,0.125);%phi_p

prior11=normpdf(param(11),1.5,0.25);%r_pi

mu=0.75;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior12= betapdf(param(12),a,b);

prior13=normpdf(param(13),0.125,0.05);
prior14=normpdf(param(14),0.125,0.05);

mu=0.625;sigma2=0.1^2;b  = sigma2/mu;a  = mu/b;
prior15=gampdf(param(15),a,b);

mu=0.25;sigma2=0.1^2;b  = sigma2/mu;a  = mu/b;
prior16=gampdf(param(16),a,b);

prior17=normpdf(param(17),0,2);
prior18=normpdf(param(18),0.4,0.1);
prior19=normpdf(param(19),0.3,0.05);

mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior20= betapdf(param(20),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior21= betapdf(param(21),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior22= betapdf(param(22),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior23= betapdf(param(23),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior24= betapdf(param(24),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior25= betapdf(param(25),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior26= betapdf(param(26),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior27= betapdf(param(27),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior28= betapdf(param(28),a,b);
%prior27=1;prior28=1;
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior29= betapdf(param(29),a,b);

[s,nu] = inverse_gamma_specification(0.1,4, 0,1, false, 'name');
% [s,nu] = inverse_gamma_specification(0.5,16, 0,1, false, 'name');
prior30= exp(lpdfig1(param(30),s,nu));
prior31= exp(lpdfig1(param(31),s,nu));
prior32= exp(lpdfig1(param(32),s,nu));
prior33= exp(lpdfig1(param(33),s,nu));
prior34= exp(lpdfig1(param(34),s,nu));

prior36= exp(lpdfig1(param(36),s,nu));
prior37= exp(lpdfig1(param(37),s,nu));

mu=0.03;sigma2=0.01^2;b  = sigma2/mu;a  = mu/b;
prior35=gampdf(param(35),a,b);%policy shock at the zlb regime


 mu=0.035;sigma2=0.015^2;b  = sigma2/mu;a  = mu/b;
prior38=gampdf(param(38),a,b);%gain

 mu=0.1;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
 prior39=betapdf(param(39),a,b);%1-p_11
  mu=0.3;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
 prior40=betapdf(param(40),a,b);%1-p_22
  prior41=normpdf(param(41),0.05,0.25);%rbar_zlb



prior=prior1*prior2*prior3*prior4*prior5*prior6*prior7*prior8*...
    prior9*prior10*prior11*prior12*...
    prior13*prior14*prior15*prior16*prior17*...
    prior18*prior19*prior20*prior21*prior22*prior23*prior24*...
    prior25*prior26*prior27*prior28*prior29*prior30*prior31*...
    prior32*prior33*prior34*prior35*prior36*prior37*prior38*...
    prior39*prior40*prior41;
prior=-log(prior);
end