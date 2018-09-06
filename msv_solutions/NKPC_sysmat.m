function[ A B C D E F G] = NKPC_sysmat(param)

y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
 rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
 eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 

A= [1,0,1/tau,-1,0;
 -gamma,1,0,0,-1;
 phi_y*(rho_r-1),phi_pinf*(rho_r-1),1,0,0;
 0,0,0,1,0;
 0,0,0,0,1];
B= [0,0,0,0,0;
 0,0,0,0,0;
 0,0,rho_r,0,0;
 0,0,0,rho_y,0;
 0,0,0,0,rho_pinf];
C=[1,1/tau,0,0,0;
 0,99/100,0,0,0;
 0,0,0,0,0;
 0,0,0,0,0;
 0,0,0,0,0];
D= [0,0,0;
 0,0,0;
 0,0,1;
 1,0,0;
 0,1,0];
E=[y_bar;pi_bar;r_bar];
F=[1,0,0,0,0;
         0,1,0,0,0;
         0,0,1,0,0];
G=[0,0,0,0,0;
0,0,0,0,0;
0,0,0,0,0];

end