function[AA BB CC DD EE FF rho E F G] = NKPC_sysmat_MSV(param);
y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
 rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
 eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 

AA= [1,0,1/tau;
 -gamma,1,0;
 phi_y*(rho_r-1),phi_pinf*(rho_r-1),1];
BB= [0,0,0;
 0,0,0;
 0,0,rho_r];

CC=[1,1/tau,0;
 0,99/100,0;
 0,0,0];
EE= [0,0,0;
 0,0,0;
 0,0,1];

DD=[1,0;0,1;0,0];

FF=[1,0,0;0,1,0];

rho=[1,0;0,1];

 E=[y_bar;pi_bar;r_bar];
 F=[1,0,0,0,0;
          0,1,0,0,0;
          0,0,1,0,0];
 G=[0,0,0,0,0;
 0,0,0,0,0;
 0,0,0,0,0];
end