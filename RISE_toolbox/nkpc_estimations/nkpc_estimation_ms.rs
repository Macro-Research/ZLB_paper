endogenous y pinf r u_y u_pi gap_cbo pinfobs robs
observables gap_cbo pinfobs robs 
exogenous eps_y eps_pi eps_r
parameters beta y_bar pi_bar kappa tau rho_y rho_pi sig_pi sig_y 
parameters(coef,2) r_bar phi_pi phi_y rho_r sig_r
parameters coef_tp_1_2, coef_tp_2_1 

model


y=y(+1)-(1/tau)*(r-pinf(+1))+u_y;
pinf=beta*pinf(+1)+kappa*y+u_pi;
r=rho_r*r(-1)+(1-rho_r)*(phi_pi*pinf+phi_y*y)+sig_r*eps_r;
u_y=rho_y*u_y(-1)+sig_y*eps_y;
u_pi=rho_pi*u_pi(-1)+sig_pi*eps_pi;
gap_cbo = y_bar+y;
pinfobs=pi_bar+pinf;
robs=r_bar+r;

parameterization;
coef_tp_1_2,0.1,0.1,0.05,beta_pdf;
coef_tp_2_1,0.3,0.3,0.1,beta_pdf;
y_bar,0.24,0,0.25,normal_pdf;
pi_bar,0.65,0.62,0.25,gamma_pdf;
r_bar(coef,1),1.05,1,0.25,gamma_pdf;
r_bar(coef,2),0.01,0.1,0.25,normal_pdf;
kappa,0.0073,0.3,0.15,beta_pdf;
tau,4,2,0.5,gamma_pdf;
phi_pi(coef,1),1.4,1.5,0.25,gamma_pdf;
phi_pi(coef,2),0;
phi_y(coef,1),0.46,0.5,0.25,gamma_pdf;
phi_y(coef,2),0;
rho_y,0.87,0.5,0.2,beta_pdf;
rho_pi,0.87,0.5,0.2,beta_pdf;
rho_r(coef,1),0.79,0.5,0.2,beta_pdf;
rho_r(coef,2),0;
sig_y,0.16,0.1,2,inv_gamma_pdf;
sig_pi,0.03,0.1,2,inv_gamma_pdf;
sig_r(coef,1),0.30,0.1,2,inv_gamma_pdf;
sig_r(coef,2),0.02,0.005,0.05,uniform_pdf;
beta,0.99;

