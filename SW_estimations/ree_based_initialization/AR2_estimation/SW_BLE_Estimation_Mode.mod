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
dx
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
beta1_i beta2_i betar_i betac_i betainv_i betay_i betapinf_i betaw_i
beta1_q beta2_q betar_q betac_q betainv_q betay_q betapinf_q betaw_q
beta1_rk beta2_rk betar_rk betac_rk betainv_rk betay_rk betapinf_rk betaw_rk
beta1_pinf beta2_pinf betar_pinf betac_pinf betainv_pinf betay_pinf betapinf_pinf betaw_pinf
beta1_c beta2_c betar_c betac_c betainv_c betay_c betapinf_c betaw_c
beta1_l beta2_l betar_l betac_l betainv_l betay_l betapinf_l betaw_l
beta1_w beta2_w betar_w betac_w betainv_w betay_w betapinf_w betaw_w
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



beta1_rk=0.9;
beta1_q=0.9;
beta1_c=0.9;
beta1_i=0.9;
beta1_l=0;
beta1_pinf=0.9;
beta1_w=0.9;

beta2_rk=0;
beta2_q=0;
beta2_c=0;
beta2_i=0;
beta2_l=0;
beta2_pinf=0;
beta2_w=0;     

betar_rk=0;
betar_q=0;
betar_c=0;
betar_i=0;
betar_l=0;
betar_pinf=0;
betar_w=0;     

betac_rk=0;
betac_q=-0;
betac_c=0;
betac_i=0;
betac_l=0;
betac_pinf=0;
betac_w=0; 

betainv_rk=0;
betainv_q=0;
betainv_c=0;
betainv_i=0;
betainv_l=0;
betainv_pinf=0;
betainv_w=0; 

betay_rk=0;
betay_q=0;
betay_c=0;
betay_i=0;
betay_l=0;
betay_pinf=0;
betay_w=0; 

betapinf_rk=0;
betapinf_q=0;
betapinf_c=0;
betapinf_i=0;
betapinf_l=0;
betapinf_pinf=0;
betapinf_w=0; 

betaw_rk=0;
betaw_q=0;
betaw_c=0;
betaw_i=0;
betaw_l=0;
betaw_pinf=0;
betaw_w=0; 

delta = 0.025;
phi_w = 1.5;
G = 0.18;
curv_p = 10;
curv_w = 10;

phi=4.3869; 
sigma_c=1.0249;
lambda=0.7922 ;
xi_w= 0.6673 ;
sigma_l= 2.3583 ;
xi_p= 0.7389;
iota_w= 0;
iota_p=0 ;
psi=0.4552;
phi_p= 1.5003 ;

r_dy=0.1;
rho=0.5;

/*load('policyParametersBLE.mat');
r_pi =policyParametersBLE(1);
r_y=policyParametersBLE(2);*/


pi_bar=0.5745 ;
beta_const=0.1756 ;
l_bar= 1.8140 ;
gamma_bar=0.4168  ;
alpha= 0.1325  ;
rho_a= 0.9743 ;
rho_b=0.3339 ;
rho_g=0.9860;
rho_i= 0.4576;
rho_r=0.3301;
rho_p=0.0981 ;
rho_w= 0.0438;
Mu_p= 0;
Mu_w=0;
rho_ga= 0;



shocks;
	var eta_a; stderr  0.4084;
	var eta_b; stderr 0.5541;
	var eta_g; stderr 0.3977;
	var eta_i; stderr 1.2428;
	var eta_r; stderr 0.1075 ;
	var eta_p; stderr 0.1871;
	var eta_w; stderr 0.8479;
end;


iota_w=0;iota_p=0;lambda=0;
Mu_p= 0;Mu_w=0;rho_ga= 0;




//original priors
estimated_params;
	phi,2.5,1.1,15,NORMAL_PDF,4,1.5;
	sigma_c,1.5,0.25,3,NORMAL_PDF,1.50,0.375;
	lambda,0.5,0.001,0.99,BETA_PDF,0.7,0.1;
	xi_w,0.57,0.3,0.99,BETA_PDF,0.75,0.05;
	sigma_l,0.96,0.25,10,NORMAL_PDF,2,0.5;
	xi_p,0.9,0.01,0.99,BETA_PDF,0.75,0.05;
	//iota_w,0.3,0.01,0.99,BETA_PDF,0.5,0.15;
	//iota_p,0.65,0.01,0.99,BETA_PDF,0.5  ,0.15;
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
	rho_g,0.5,0,.9999,BETA_PDF,0.5,0.2;
    rho_i,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_r,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_p,0.5,0,.99,BETA_PDF,0.5,0.2;
    rho_w,0.5,0,.99,BETA_PDF,0.5,0.2;
//Mu_p,0.5,0,.99,BETA_PDF,0.5,0.2;
//Mu_w,0.5,0,.99,BETA_PDF,0.5,0.2;
//rho_ga,0.5,0,.99,BETA_PDF,0.5,0.25;

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
	


#forecast_i=beta1_i^2*i(-1)+betar_i*r(-1)+betac_i*c(-1)+betainv_i*i(-1)+betay_i*y(-1)+betapinf_i*pinf(-1)+betaw_i*w(-1);
#forecast_q=beta1_q^2*q(-1)+betar_q*r(-1)+betac_q*c(-1)+betainv_q*i(-1)+betay_q*y(-1)+betapinf_q*pinf(-1)+betaw_q*w(-1);
#forecast_rk=beta1_rk^2*rk(-1)+betar_rk*r(-1)+betac_rk*c(-1)+betainv_rk*i(-1)+betay_rk*y(-1)+betapinf_rk*pinf(-1)+betaw_rk*w(-1);
#forecast_pinf=beta1_pinf^2*pinf(-1)+betar_pinf*r(-1)+betac_pinf*c(-1)+betainv_pinf*i(-1)+betay_pinf*y(-1)+betapinf_pinf*pinf(-1)+betaw_pinf*w(-1);
#forecast_c=beta1_c^2*c(-1)+betar_c*r(-1)+betac_c*c(-1)+betainv_c*i(-1)+betay_c*y(-1)+betapinf_c*pinf(-1)+betaw_c*w(-1);
#forecast_l=beta1_l^2*l(-1)+betar_l*r(-1)+betac_l*c(-1)+betainv_l*i(-1)+betay_l*y(-1)+betapinf_l*pinf(-1)+betaw_l*w(-1);
#forecast_w=beta1_w^2*w(-1)+betar_w*r(-1)+betac_w*c(-1)+betainv_w*i(-1)+betay_w*y(-1)+betapinf_w*pinf(-1)+betaw_w*w(-1);

	mc = alpha*rk + (1-alpha)*w - eps_a;
	zcap =  ((1 - psi)/psi) * rk;
	rk =  w + l - k_s;
	k_s =  k(-1) + zcap;
	i = (1/(1 + beta_bar*gamma)) * (i(-1) + (beta_bar * gamma) * forecast_i + (1/(gamma^2*phi)) * q) + eps_i;
	q = ((1-delta)/(Rk+(1-delta)))*forecast_q + (Rk/(Rk+(1-delta))) * forecast_rk - r + forecast_pinf + (1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma)))) * eps_b ;
	c = (lambda/gamma)/(1+lambda/gamma) * c(-1) + (1/(1+lambda/gamma)) * forecast_c + ((sigma_c-1)*WL_C/(sigma_c*(1+lambda/gamma))) * (l - forecast_l) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma)) * (r - forecast_pinf) + eps_b;
	y = C_Y * c + I_Y * i + eps_g + Z_Y * zcap;
	y = phi_p * (alpha * k_s + (1-alpha) * l + eps_a);
	pinf = (1/(1+beta_bar*gamma*iota_p)) * (beta_bar*gamma*forecast_pinf + iota_p * pinf(-1) + ((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curv_p+1) * mc) + eps_p ; 
	w =  (1/(1+beta_bar*gamma))*w(-1)
	   +(beta_bar*gamma/(1+beta_bar*gamma))*forecast_w
	   +(iota_w/(1+beta_bar*gamma))*pinf(-1)
	   -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf
	   +(beta_bar*gamma)/(1+beta_bar*gamma)*forecast_pinf
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
dx=outputGap-outputGap(-1);
end;

 estimation(optim=('MaxIter',200),
//datafile=usmodel_data,
datafile=raf_dataset,
mode_compute=1,
nograph,
	nodiagnostic,
	optim=('Algorithm','active-set'),
//first_obs=120,//estimation smpl in ble paper
//nobs=88,
first_obs=147,
nobs=60,
//nobs=76,
//nobs=75,
//nobs=172,
//nobs=156,
//first_obs=208,
presample=4,
kalman_algo=1,
lik_init=1,
prefilter=0,
mh_replic=00000,
mh_nblocks=1,
mh_jscale=0.35,
mh_drop=0.2);
//shock_decomposition(parameter_set=posterior_mode)pinf,y;
//forecast=1);  //labobs robs pinfobs dy dc dinve dw; 
//stoch_simul(ar=10,conditional_variance_decomposition=[1:100]);
stoch_simul(ar=10,irf=100);

//stoch_simul(ar=10);
//options_.noprint=1;
//clean_current_folder;