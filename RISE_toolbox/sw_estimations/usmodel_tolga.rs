var   labobs robs pinfobs dy dc dinve dw 
 mc zcap rk k q    c inve y lab pinf w r eps_a  eps_b eps_g eps_i  eps_r  eps_w eps_p kp

observables labobs robs pinfobs dy dc dinve dw

exogenous eta_a eta_b eta_g eta_i eta_r eta_p eta_w
 

parameters curvw rho_ga curvp l_bar pi_bar beta_const mu_w mu_p alpha 
psi beta phi delta sigma_c lambda   phi_p
iota_w xi_w iota_p xi_p sigma_l phi_w 
r_pi r_dy r_y rho 
rho_a  rho_b rho_g  rho_i rho_r rho_p rho_w  
gamma_bar 
r_bar G
sig_a sig_b sig_g sig_i sig_r sig_p sig_w



parameterization


delta,.025;
phi_w,1.5;
G,0.18;
curvp,10;
curvw,10;
    xi_w,0.8,0.5,0.1,beta_pdf;
	xi_p,0.8,0.5,0.1,beta_pdf;
	iota_w,0.65,0.5,0.15,beta_pdf;
	iota_p,0.3,0.5,0.15,beta_pdf;
	phi_p,1.5,1.25,0.125,normal_pdf;
	alpha,0.17,0.3,0.05,normal_pdf;
	psi,0.75,0.5,0.15,beta_pdf;
    lambda,0.7,0.7,0.1,beta_pdf;
    sigma_c,1.3,1.5,0.375,normal_pdf;
    phi,4.5,4,1.5,normal_pdf;
    sigma_l,1.96,2,0.75,normal_pdf;


	rho_a,0.95 ,0.5,0.2,beta_pdf;
    rho_b,0.5,0.5,0.2,beta_pdf;
	rho_g,0.95,0.5,0.2,beta_pdf;
    rho_i,0.8,0.5,0.2,beta_pdf;
    rho_r,0.5,0.5,0.2,beta_pdf;
    rho_p,0.6,0.5,0.2,beta_pdf;
    rho_w,0.6,0.5,0.2,beta_pdf;
    rho_ga,0.5,0.5,0.2,beta_pdf;
    mu_p,0.5,0.5,0.2,beta_pdf;
    mu_w,0.5,0.5,0.2,beta_pdf;
    
    r_pi,1.8,1.5,0.25,normal_pdf;
	rho,0.8258,0.75,0.1,beta_pdf;
	r_y,0.08,0.125,0.05,normal_pdf;
	r_dy,0.23,0.125,0.05,normal_pdf;

	pi_bar,0.7,0.625,0.1,gamma_pdf;
	l_bar,0,0,2,normal_pdf;
	beta_const,0.19,0.25,0.1,gamma_pdf;
    gamma_bar,0.4,0.4,0.1,normal_pdf;

    sig_a,0.5,0.1,2,inv_gamma_pdf;
	sig_b,1,0.1,2,inv_gamma_pdf;
	sig_g,0.55,0.1,2,inv_gamma_pdf;
	sig_i,0.3,0.1,2,inv_gamma_pdf;
	sig_r,0.2,0.1,2,inv_gamma_pdf;
    sig_p,0.1,0.1,2,inv_gamma_pdf;
	sig_w,0.4,0.1,2,inv_gamma_pdf;




model
# gamma=gamma_bar/100 + 1;
# cpie = pi_bar / 100 + 1;
# beta_bar=beta*gamma^(-sigma_c);
# cr=cpie/(beta*gamma^(-sigma_c));
# crk=(beta^(-1))*(gamma^sigma_c) - (1-delta);
# cw = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*crk^alpha))^(1/(1-alpha));
# cikbar=(1-(1-delta)/gamma);
# cik=(1-(1-delta)/gamma)*gamma;
# clk=((1-alpha)/alpha)*(crk/cw);
# cky=phi_p*(clk)^(alpha-1);
# ciy=cik*cky;
# ccy=1-G-cik*cky;
# crkky=crk*cky;
# cwhlc=(1/phi_w)*(1-alpha)/alpha*crk*cky/ccy;
# cwly=1-crk*cky;


	
	      mc =  alpha*rk+(1-alpha)*(w) - 1*eps_a - 0*(1-alpha)*eps_a ;
	      zcap =  (1/(psi/(1-psi)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+beta_bar*gamma))* (  inve(-1) + beta_bar*gamma*inve(1)+(1/(gamma^2*phi))*q ) +eps_i ;
          q = -r+pinf(1)-0*eps_b +(1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma))))*eps_b + (crk/(crk+(1-delta)))*rk(1) +  ((1-delta)/(crk+(1-delta)))*q(1) ;
	      c = (lambda/gamma)/(1+lambda/gamma)*c(-1) + (1/(1+lambda/gamma))*c(+1) +((sigma_c-1)*cwhlc/(sigma_c*(1+lambda/gamma)))*(lab-lab(+1)) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma))*(r-pinf(+1) + 0*eps_b) +eps_b ;
	      y = ccy*c+ciy*inve+eps_g  +  1*crkky*zcap ;
	      y = phi_p*( alpha*k+(1-alpha)*lab +eps_a );
	      pinf =  (1/(1+beta_bar*gamma*iota_p)) * ( beta_bar*gamma*pinf(1) +iota_p*pinf(-1) 
               +((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curvp+1)*(mc)  )  + eps_p; 
	      w =  (1/(1+beta_bar*gamma))*w(-1)
               +(beta_bar*gamma/(1+beta_bar*gamma))*w(1)
               +(iota_w/(1+beta_bar*gamma))*pinf(-1)
               -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf
               +(beta_bar*gamma)/(1+beta_bar*gamma)*pinf(1)
               +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curvw+1))*
               (sigma_l*lab + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*c(-1) -w) 
               + 1*eps_w ;
	      r =  r_pi*(1-rho)*pinf
               +r_y*(1-rho)*(y-phi_p*eps_a)     
               +r_dy*(y-phi_p*eps_a-y(-1)+phi_p*eps_a(-1))
               +rho*r(-1)
               +eps_r  ;
	      eps_a = rho_a*eps_a(-1)  + sig_a*eta_a;
	      eps_b = rho_b*eps_b(-1) + sig_b*eta_b;
	      eps_g = rho_g*(eps_g(-1)) + sig_g*eta_g + rho_ga*sig_a*eta_a;
	      eps_i = rho_i*eps_i(-1) + sig_i*eta_i;
	      eps_r = rho_r*eps_r(-1) + sig_r*eta_r;
	      eps_p = rho_p*eps_p(-1) + sig_p*eta_p - mu_p*sig_p*eta_p(-1);
	      eps_w = rho_w*eps_w(-1) + sig_w*eta_w - mu_w*sig_w*eta_w(-1) ;
	         
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*gamma^2*phi*eps_i ;

// measurement equations

dy=y-y(-1)+gamma_bar;
dc=c-c(-1)+gamma_bar;
dinve=inve-inve(-1)+gamma_bar;
dw=w-w(-1)+gamma_bar;
pinfobs = 1*(pinf) + pi_bar;
robs =    1*(r) + r_bar;
labobs = lab + l_bar;



