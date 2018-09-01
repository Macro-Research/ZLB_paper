var y pi r u_y u_pi   gap_hp pinfobs robs; 
varexo eps_y eps_pi eps_r;
parameters y_bar pi_bar r_bar kappa tau phi_pi phi_y rho_y rho_pi rho_r 
a_1 a_2 a_3
b_11 b_12 b_13 b_21 b_22 b_23 b_31 b_32 b_33;

y_bar=0;pi_bar=0;r_bar=0; kappa= 0.01 ; tau=3; phi_pi=1.5; phi_y = 0.5; rho_r =0; rho_y =0.9 ; rho_pi =0.9; 

    
shocks;
var eps_y ; stderr 0.7;
var eps_pi; stderr 0.3 ;
var eps_r ; stderr 0;
end;

estimated_params;
y_bar,0.8,0,2,NORMAL_PDF,0.4,0.25;
pi_bar,0.7,-2,2,GAMMA_PDF,0.62,0.25;
r_bar,0.9,-2,2,GAMMA_PDF,0.5,0.25;
kappa,0.1,0,1,BETA_PDF,0.3,0.15;
tau,5,0,10,GAMMA_PDF,2,0.5;
phi_pi,1.5,0,10,GAMMA_PDF,1.5,0.25;
phi_y,0.5,0,10,GAMMA_PDF,0.5,0.25;
rho_y,0.5,0,1,BETA_PDF,0.5,0.2;
rho_pi,0.5,0,1,BETA_PDF,0.5,0.2;
rho_r,0.5,0,1,BETA_PDF,0.5,0.2;
stderr eps_y,0.07,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eps_pi,1,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eps_r,0.4,0.01,3,INV_GAMMA_PDF,0.1,2;

a_1,0,-5,5,NORMAL_PDF,0,1;
a_2,0,-5,5,NORMAL_PDF,0,1;
a_3,0,-5,5,NORMAL_PDF,0,1;

b_11,0.95,0,0.9999,NORMAL_PDF,0.65,1;
b_12,0,-0.9999,0.9999,NORMAL_PDF,0.51,1;
b_13,0,-0.9999,0.9999,NORMAL_PDF,0.31,1;
b_21,0,-0.9999,0.9999,NORMAL_PDF,-0.05,1;
b_22,0.95,,0.9999,NORMAL_PDF,0.89,1;
b_23,0,-0.9999,0.9999,NORMAL_PDF,0.08,1;
b_31,0,-0.9999,0.9999,NORMAL_PDF,0.01,1;
b_32,0,-0.9999,0.9999,NORMAL_PDF,0.25,1;
b_33,0.95,0,0.9999,NORMAL_PDF,0.51,1;

end;
varobs gap_hp pinfobs robs;

model(linear);
#beta=0.99;
#y_forecast=(a_1 + a_1*b_11 + a_2*b_12 + a_3*b_13)+
 (b_11^2 + b_12*b_21 + b_13*b_31)*y(-1)+
(b_11*b_12 + b_12*b_22 + b_13*b_32)*pi*(-1)+
(b_11*b_13 + b_12*b_23 + b_13*b_33)*r(-1);

#pinf_forecast=(a_2 + a_1*b_21 + a_2*b_22 + a_3*b_23)+
(b_11*b_21 + b_21*b_22 + b_23*b_31)*y(-1)+
(b_22^2 + b_12*b_21 + b_23*b_32)*pi(-1)+
(b_13*b_21 + b_22*b_23 + b_23*b_33)*r(-1);

y=y_forecast-(1/tau)*(r-pinf_forecast)+u_y;
pi=beta*pinf_forecast+kappa*y+u_pi;
r=rho_r*r(-1)+(1-rho_r)*(phi_pi*pi+phi_y*y)+eps_r;
u_y=rho_y*u_y(-1)+eps_y;
u_pi=rho_pi*u_pi(-1)+eps_pi;
gap_hp = y_bar+y;
pinfobs=pi_bar+pi;
robs=r_bar+r;
end;
//resid(1);
//steady;


estimation(datafile=us_dataset,
mode_compute=1,
nograph,nodiagnostic,
optim=('Algorithm','active-set'),
first_obs=44,
//nobs=79,
mh_replic=0,mh_jscale=0.51,mh_drop=0.2,mh_nblocks=1);
stoch_simul(ar=10,nograph,periods=10000);
//clean_current_folder;

