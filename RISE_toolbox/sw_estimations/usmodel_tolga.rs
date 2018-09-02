endogenous labobs robs   pinfobs  dy      dc      dinve  dw      mc    zcap  rk    k_s   q   c  i y     l   pinf   w       
 r      eps_a   eps_b eps_g eps_i  eps_r eps_p eps_w k 

observables labobs robs pinfobs dy dc dinve dw

exogenous eta_a eta_b eta_g eta_i eta_r eta_p eta_w
 

parameters curv_w curv_p G  delta  phi_w 
l_bar  pi_bar beta_const alpha gamma_bar 
psi  phi  lambda phi_p iota_w xi_w iota_p xi_p 
sigma_c sigma_l r_pi r_dy r_y rho 
rho_a rho_b rho_p rho_w  rho_i  rho_r  rho_ga  rho_g Mu_w Mu_p 
sig_a sig_b sig_g sig_i sig_r sig_p sig_w


parameterization
	phi,2.5,4,1.5,normal_pdf;
	sigma_c,1.5,1.5,0.375,normal_pdf;
	lambda,0.5,0.7,0.1,beta_pdf;
	xi_w,0.57,0.75,0.1,beta_pdf;
	sigma_l,0.96,2,0.75,normal_pdf;
	xi_p,0.9,0.75,0.1,beta_pdf;
	iota_w,0.3,0.5,0.1,beta_pdf;
	iota_p,0.65,0.5,0.1,beta_pdf;
	psi,0.3,0.5,0.15,beta_pdf;
	phi_p,1.3,1.25,0.125,normal_pdf;
	r_pi,1.2,1.5,0.25,normal_pdf;
	rho,0.8258,0.75,0.1,beta_pdf;
	r_y,0.13,0.125,0.05,normal_pdf;
	r_dy,0.12,0.125,0.05,normal_pdf;
	pi_bar,0.7,0.625,0.1,gamma_pdf;
	beta_const,0.19,0.25,0.1,gamma_pdf;
	l_bar,-2.25,0,2,normal_pdf;
	gamma_bar,0.4,0.4,0.1,normal_pdf;
	alpha,0.16,0.3,0.1,normal_pdf;

	rho_a,0.5 ,0.5,0.2,beta_pdf;
    rho_b,0.5,0.5,0.2,beta_pdf;
	rho_g,0.5,0.5,0.2,beta_pdf;
    rho_i,0.5,0.5,0.2,beta_pdf;
    rho_r,0.5,0.5,0.2,beta_pdf;
    rho_p,0.5,0.5,0.2,beta_pdf;
    rho_w,0.5,0.5,0.2,beta_pdf;
Mu_p,0.5,0.5,0.2,beta_pdf;
Mu_w,0.5,0.5,0.2,beta_pdf;
rho_ga,0.5,0.5,0.2,beta_pdf;

	sig_a,0.1,0.1,2,inv_gamma_pdf;
	sig_b,0.77,0.1,2,inv_gamma_pdf;
	sig_g,0.55,0.1,2,inv_gamma_pdf;
	sig_i,1.3,0.1,2,inv_gamma_pdf;
	sig_r,0.2,0.1,2,inv_gamma_pdf;
    sig_p,0.4,0.1,2,inv_gamma_pdf;
	sig_w,0.3,0.1,2,inv_gamma_pdf;


model;
	// Transformation of estimated parameters to model parameters
	# PI_star = 1 + pi_bar/100;
	# gamma = 1 + gamma_bar/100 ;
	# beta = 1/(1 + beta_const/100);
	// Convenience variable
	# beta_bar = beta*gamma^(-sigma_c);
	// Steady state values
	# Rk = (beta^(-1)) * (gamma^sigma_c) - (1-delta);
	# W = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*Rk^alpha))^(1/(1-alpha));
	# I_K_bar = (1-(1-delta)/gamma);
	# I_K = (1-(1-delta)/gamma)*gamma;
	# L_K = ((1-alpha)/alpha)*(Rk/W);
	# K_Y = phi_p*(L_K)^(alpha-1);
	# I_Y = I_K * K_Y;
	# C_Y = 1 - G - I_K*K_Y;
	# Z_Y = Rk*K_Y;
	# WL_C = (1/phi_w)*(1-alpha)/alpha*Rk*K_Y/C_Y;
	# r_bar=((PI_star/(beta*gamma^(-sigma_c)))-1)*100;
	
	
	mc = alpha*rk + (1-alpha)*w - eps_a;
	zcap =  ((1 - psi)/psi) * rk;
	rk =  w + l - k_s;
	k_s =  k(-1) + zcap;
	i = (1/(1 + beta_bar*gamma)) * (i(-1) + (beta_bar * gamma) * i(1) + (1/(gamma^2*phi)) * q) + eps_i;
	q = ((1-delta)/(Rk+(1-delta)))*q(1) + (Rk/(Rk+(1-delta))) * rk(1) - r + pinf(+1) + (1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma)))) * eps_b ;
	c = (lambda/gamma)/(1+lambda/gamma) * c(-1) + (1/(1+lambda/gamma)) * c(+1) + ((sigma_c-1)*WL_C/(sigma_c*(1+lambda/gamma))) * (l - l(+1)) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma)) * (r - pinf(+1)) + eps_b;
	y = C_Y * c + I_Y * i + eps_g + Z_Y * zcap;
	y = phi_p * (alpha * k_s + (1-alpha) * l + eps_a);
	pinf = (1/(1+beta_bar*gamma*iota_p)) * (beta_bar*gamma*pinf(+1) + iota_p * pinf(-1) + ((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curv_p+1) * mc) + eps_p ; 
	w =  (1/(1+beta_bar*gamma))*w(-1)
	   +(beta_bar*gamma/(1+beta_bar*gamma))*w(1)
	   +(iota_w/(1+beta_bar*gamma))*pinf(-1)
	   -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf
	   +(beta_bar*gamma)/(1+beta_bar*gamma)*pinf(1)
	   +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curv_w+1))*
	   (sigma_l*l + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*c(-1) -w) 
	   + 1*eps_w ;
	r =  r_pi * (1-rho) * pinf + r_y * (1-rho) * (y-phi_p*eps_a) + r_dy * ( y - phi_p*eps_a - (y(-1) - phi_p*eps_a(-1))) + rho * r(-1) + eps_r;
	eps_a = rho_a * eps_a(-1) + sig_a*eta_a;
	eps_g = rho_g * eps_g(-1) + sig_g*eta_g + rho_ga * sig_a*eta_a;
	eps_i = rho_i * eps_i(-1) + sig_i*eta_i;
    eps_b = rho_b * eps_b(-1) + sig_b*eta_b; 
    eps_p = rho_p * eps_p(-1)  + sig_p*eta_p - Mu_p*sig_p*eta_p(-1);
    eps_w = rho_w * eps_w(-1) +sig_w*eta_w -Mu_w*sig_w*eta_w(-1);
   eps_r = rho_r*eps_r(-1) + sig_r*eta_r;

    k = (1-I_K_bar) * k(-1) + I_K_bar * i + I_K_bar*gamma^2*phi*eps_i;

  dy = y - y(-1) + gamma_bar;
	dc = c - c(-1) + gamma_bar;
	dinve = i - i(-1) + gamma_bar;
	dw = w - w(-1) + gamma_bar;
	pinfobs = pinf + pi_bar;
	robs = r + r_bar;
	labobs = l + l_bar;




