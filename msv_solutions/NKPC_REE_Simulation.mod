var y pi r u_y u_pi ;
//   gap_hp pinfobs robs; 
varexo eps_y eps_pi eps_r;
parameters y_bar pi_bar r_bar kappa tau phi_pi phi_y rho_y rho_pi rho_r ;

y_bar=0;pi_bar=0;r_bar=0; kappa= 0.01 ; tau=3; phi_pi=1.5; phi_y = 0.5; rho_r =0.9; rho_y =0.5 ; rho_pi =0.5; 

    
shocks;
var eps_y ; stderr 0.7;
var eps_pi; stderr 0.3 ;
var eps_r ; stderr 0.3;
end;

/*estimated_params;
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
end;*/
//varobs gap_hp pinfobs robs;

model(linear);
#beta=0.99;

y=y(+1)-(1/tau)*(r-pi(+1))+u_y;
pi=beta*pi(+1)+kappa*y+u_pi;
r=rho_r*r(-1)+(1-rho_r)*(phi_pi*pi+phi_y*y)+eps_r;
u_y=rho_y*u_y(-1)+eps_y;
u_pi=rho_pi*u_pi(-1)+eps_pi;
//gap_hp = y_bar+y;
//pinfobs=pi_bar+pi;
//robs=r_bar+r;
end;
//resid(1);
//steady;


/*estimation(datafile=full_dataset,
mode_compute=1,
nograph,nodiagnostic,
optim=('Algorithm','active-set'),
first_obs=44,
//nobs=79,
mh_replic=0,mh_jscale=0.51,mh_drop=0.2,mh_nblocks=1);*/
stoch_simul(ar=10,nograph);
//clean_current_folder;

