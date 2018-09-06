var

	labobs $log(l^{obs})$
    robs   $d(log(r^{obs})) $
    pinfobs $d(log(\pi^{obs}))$ 
    dy      $ d( log(y^{obs})) $
    dc      $ d( log(c^{obs})) $
    dinve   $ d( log(i^{obs})) $
    dw      $ d( log(w^{obs})) $
	mc      $mc$
   zcap     $z$
   rk       $rk$
   k_s      $ks$
   q       $q$
   c       $c$
   i        $i$
   y        $y$
   l        $l$
  pinf      $\pi$
  w         $w$
  r          $r$
 outputGap   $Gap$
	eps_a    $\eps_a$
eps_b $\eps_b$
eps_g $\eps_g$

eps_i $\eps_i$
eps_r $\eps_r$
eps_p  $\eps_p$
eps_w $\eps_w$
k     $k$
;
 

parameters

curv_w $\eps_w$
curv_p $\eps_p$
G $G$
delta $\delta$
phi_w $\phi_w$

l_bar  $\bar{l}$
pi_bar $\bar{\pi}$
beta_const $\bar{\beta}$
alpha $\alpha$
gamma_bar $\bar{\gamma}$

psi $\psi$
phi $\phi$
lambda $\lambda$
phi_p $\phi_p$
iota_w $\iota_w$
xi_w $\xi_w$
iota_p $\iota_p$
xi_p $\xi_p$

sigma_c $\sigma_c$
sigma_l $\sigma_l$

r_pi $\r_{\pi}$
r_dy $r_{\delta y}$
r_y $r_y$
rho $\rho$

rho_a $\rho_a$
rho_b $\rho_b$
rho_p $\rho_p$
rho_w  $\rho_w$
rho_i  $\rho_i$
rho_r  $\rho_r$
rho_ga  $\rho_{\ga}$
rho_g   $\rho_g$
Mu_w   $\mu_w$
Mu_p  $\mu_p$
; 

varexo
eta_a $\eta_a$
eta_b $\eta_b$
eta_g  $\eta_g$
eta_i $\eta_i$
eta_r $\eta_r$
eta_p $\eta_p$
eta_w  $\eta_w$
;

delta = 0.025;
phi_w = 1.5;
G = 0.18;
curv_p = 10;
curv_w = 10;

phi=7.5493;
sigma_c=1.5607;
lambda=0.7803 ;
xi_w= 0.8991 ;
sigma_l= 1.4122 ;
xi_p= 0.9166;
iota_w= 0;
iota_p=0 ;
psi=0.6147;
phi_p= 1.4584 ;



r_dy=0.1 ;
rho=0.5;

pi_bar=0.6281 ;
beta_const=0.1201 ;
l_bar= 1.0740 ;
gamma_bar=0.3005  ;
alpha= 0.1671  ;
rho_a= 0.9872 ;
rho_b=0.2971 ;
rho_g=0.9552;
rho_i= 0.8086;
rho_r=0.2708;
rho_p=0.4653 ;
rho_w= 0.0450;
Mu_p= 0;
Mu_w=0;
rho_ga= 0;


shocks;
	var eta_a; stderr  0.4084;
	var eta_b; stderr 0.1900;
	var eta_g; stderr 0.3971;
	var eta_i; stderr 0.2868;
	var eta_r; stderr 0.1141 ;
	var eta_p; stderr 0.1082;
	var eta_w; stderr 0.4205;
end;


//original priors
estimated_params;
	phi,2.5,1.1,15,NORMAL_PDF,4,1.5;
	sigma_c,1.5,0.25,3,NORMAL_PDF,1.50,0.375;
	lambda,0.5,0.001,0.99,BETA_PDF,0.7,0.1;
	xi_w,0.57,0.3,0.99,BETA_PDF,0.5,0.1;
	sigma_l,0.96,0.25,10,NORMAL_PDF,2,0.75;
	xi_p,0.9,0.5,0.99,BETA_PDF,0.5,0.1;
	iota_w,0.3,0.01,0.99,BETA_PDF,0.5,0.15;
	iota_p,0.65,0.01,0.99,BETA_PDF,0.5  ,0.15;
	psi,0.3,0.01,1,BETA_PDF,0.5,0.15;
	phi_p,1.3,1.0,3,NORMAL_PDF,1.25,0.125;
	r_pi,1.2,1.0,3,NORMAL_PDF,1.5,0.25;
	rho,0.8258,0.5,0.975,BETA_PDF,0.75,0.10;
	r_y,0.13,0.001,0.5,NORMAL_PDF,0.125,0.05;
	r_dy,0.12,0.001,0.5,NORMAL_PDF,0.125,0.05;
	pi_bar,0.7,0.1,2.0,GAMMA_PDF,0.625,0.1;
	beta_const,0.19,0.01,2.0,GAMMA_PDF,0.25,0.1;
	l_bar,-2.25,-10.0,10.0,NORMAL_PDF,0.0,2.0;
	gamma_bar,0.4,0.1,0.8,NORMAL_PDF,0.4,0.10;	
	alpha,0.16,0.01,1.0,NORMAL_PDF,0.3,0.05;

	rho_a,0.5 ,0,.99,BETA_PDF,0.5,0.2;
    rho_b,0.5,0,.99,BETA_PDF,0.5,0.2;
	rho_g,0.5,0,.999999,BETA_PDF,0.5,0.2;
    rho_i,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_r,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_p,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_w,0.5,0,.99,BETA_PDF,0.5,0.2;
Mu_p,0.5,0,.99,BETA_PDF,0.5,0.2;
Mu_w,0.5,0,.99,BETA_PDF,0.5,0.2;
rho_ga,0.5,0,.99,BETA_PDF,0.5,0.25;

	stderr eta_a,0.55,0.01,3,INV_GAMMA_PDF,0.1,2;
	stderr eta_b,0.77,0.025,5,INV_GAMMA_PDF,0.1,2;
	stderr eta_g,0.55,0.01,3,INV_GAMMA_PDF,0.1,2;
	stderr eta_i,1.3,0.01,3,INV_GAMMA_PDF,0.1,2;
	stderr eta_r,0.2,0.01,3,INV_GAMMA_PDF,0.1,2;
	stderr eta_p,0.4,0.01,3,INV_GAMMA_PDF,0.1,2;
	stderr eta_w,0.3,0.01,3,INV_GAMMA_PDF,0.1,2;

end;



varobs dy dc dinve labobs pinfobs dw robs;

model(linear); 
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
	
	
	mc = alpha*rk + (1-alpha)*w - eps_a;//ok
	zcap =  ((1 - psi)/psi) * rk;//ok
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
	eps_a = rho_a * eps_a(-1) + eta_a;
	eps_g = rho_g * eps_g(-1) + eta_g + rho_ga * eta_a;
	eps_i = rho_i * eps_i(-1) + eta_i;
    eps_b = rho_b * eps_b(-1) + eta_b; 
    eps_p = rho_p * eps_p(-1)  + eta_p - Mu_p*eta_p(-1);
    eps_w = rho_w * eps_w(-1) + eta_w -Mu_w*eta_w(-1);
   eps_r = rho_r*eps_r(-1) + eta_r;

    k = (1-I_K_bar) * k(-1) + I_K_bar * i + I_K_bar*gamma^2*phi*eps_i;
outputGap = y - phi_p*eps_a;

  dy = y - y(-1) + gamma_bar;
	dc = c - c(-1) + gamma_bar;
	dinve = i - i(-1) + gamma_bar;
	dw = w - w(-1) + gamma_bar;
	pinfobs = pinf + pi_bar;
	robs = r + r_bar;
	labobs = l + l_bar;
end;

  estimation(optim=('MaxIter',200),nograph,datafile=usmodel_data,mode_compute=3,first_obs=71,presample=4,lik_init=2,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.20,mh_drop=0.2);


options_.noprint=1;
//clean_current_folder;