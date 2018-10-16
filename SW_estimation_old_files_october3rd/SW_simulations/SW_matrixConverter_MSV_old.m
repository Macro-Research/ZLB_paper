function [ AA BB CC DD EE RHO FF GG] = SW_matrixConverter_MSV_old( parameters )




syms mc zcap rk k_s q c i1 y l pinf w r k ...%sticky-contemporaneous
    mcp zcapp rkp k_sp qp cp ip yp lp pinfp wp rp kp ...%sticky-expectations
    mcm zcapm rkm k_sm qm cm im ym lm pinfm wm rm km...%sticky-lagged
    eps_a eps_b eps_g eps_i eps_r eps_p eps_w...%structural shocks
    epsm_a epsm_b epsm_g epsm_i epsm_r epsm_p epsm_w...
    epsp_a epsp_b epsp_g epsp_i epsp_r epsp_p epsp_w...
    eta1_a eta1_b eta1_g eta1_i eta1_r eta1_p eta1_w...
    etam1_a etam1_b etam1_g etam1_i etam1_r etam1_p etam1_w;
    


NumVar= 20;
NumShock=7;

 delta   = parameters(1);
 G     = parameters(2);
 phi_w   =parameters(3); 
 curv_p  = parameters(4);
 curv_w  = parameters(5);
 
phi     =parameters(6) ;    
sigma_c =  parameters(7) ;
lambda   =  parameters(8)   ;
xi_w     =  parameters(9)    ;
sigma_l  = parameters(10);
xi_p     = parameters(11) ;
iota_w   =    parameters(12);
iota_p   =  parameters(13);
psi      =   parameters(14) ;
phi_p    =  parameters(15)  ;

r_pi     =     parameters(16);
rho      =  parameters(17) ;
r_y      = parameters(18) ;
r_dy=     parameters(19);




pi_bar   = parameters(20);
beta_const=   parameters(21);
l_bar    =parameters(22) ;
gamma_bar =parameters(23);
alpha    =   parameters(24);

rho_a= parameters(25);
rho_b =  parameters(26);
rho_g=    parameters(27);
rho_i=  parameters(28); 
rho_r =parameters(29) ;
rho_p=  parameters(30) ;
rho_w= parameters(31);
Mu_p= parameters(32) ;
Mu_w =parameters(33) ;
rho_ga=     parameters(34)  ;








PI_star = 1 + pi_bar/100;
gamma = 1 + gamma_bar/100 ;
beta = 1/(1 + beta_const/100);
beta_bar = beta*gamma^(-sigma_c);
Rk = (beta^(-1)) * (gamma^sigma_c) - (1-delta);
W = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*Rk^alpha))^(1/(1-alpha));
I_K_bar = (1-(1-delta)/gamma);
I_K = (1-(1-delta)/gamma)*gamma;
L_K = ((1-alpha)/alpha)*(Rk/W);
K_Y = phi_p*(L_K)^(alpha-1);
I_Y = I_K * K_Y;
C_Y = 1 - G - I_K*K_Y;
Z_Y = Rk*K_Y;
WL_C = (1/phi_w)*(1-alpha)/alpha*Rk*K_Y/C_Y;
r_bar=((PI_star/(beta*gamma^(-sigma_c)))-1)*100;

fs1=-mc + alpha*rk + (1-alpha)*w - eps_a==0;
	fs2=-zcap +  ((1 - psi)/psi) * rk==0;
	fs3=-rk +  w + l - k_s==0;
	fs4=-k_s +  km + zcap==0;
	fs5=-i1 + (1/(1 + beta_bar*gamma)) * (im + (beta_bar * gamma) * ip + (1/(gamma^2*phi)) * q) + eps_i==0;
	fs6=-q + ((1-delta)/(Rk+(1-delta)))*qp + (Rk/(Rk+(1-delta))) * rkp - ...
        r + pinfp + (1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma)))) * eps_b==0 ;
	fs7=-c + (lambda/gamma)/(1+lambda/gamma) * cm + (1/(1+lambda/gamma)) * cp +...
        ((sigma_c-1)*WL_C/(sigma_c*(1+lambda/gamma))) * (l - lp) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma)) * (r - pinfp) + eps_b==0;
	fs8=-y + C_Y * c + I_Y * i1 + eps_g + Z_Y * zcap==0;
	fs9=-y + phi_p * (alpha * k_s + (1-alpha) * l + eps_a)==0;
	fs10=-pinf + (1/(1+beta_bar*gamma*iota_p)) * (beta_bar*gamma*pinfp +...
        iota_p * pinfm + ((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curv_p+1) * mc) + eps_p ==0; 
	fs11=-w +  (1/(1+beta_bar*gamma))*wm...
	   +(beta_bar*gamma/(1+beta_bar*gamma))*wp...
	   +(iota_w/(1+beta_bar*gamma))*pinfm...
	   -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf...
	   +(beta_bar*gamma)/(1+beta_bar*gamma)*pinfp...
	   +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curv_w+1))*...
	   (sigma_l*l + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*cm -w)... 
	   + 1*eps_w ==0;
	fs12=-r +  r_pi * (1-rho) * pinf + r_y * (1-rho) * (y-phi_p*eps_a) + r_dy * ( y - phi_p*eps_a - (ym - phi_p*epsm_a)) + rho * rm + eps_r==0;
    fs13=-k + (1-I_K_bar) * km + I_K_bar * i1 + I_K_bar*gamma^2*phi*eps_i==0;
    

    
    
    err1=-eps_a + rho_a * epsm_a + eta1_a==0;
	err2=-eps_b + rho_b * epsm_b + eta1_b==0;
	err3=-eps_g + rho_g * epsm_g + eta1_g + rho_ga * eta1_a==0;
	err4=-eps_i + rho_i * epsm_i + eta1_i==0;
	err5=-eps_r + rho_r * epsm_r + eta1_r==0;
	err6=-eps_p + rho_p * epsm_p + eta1_p - Mu_p * etam1_p==0;
	err7=-eps_w + rho_w * epsm_w + eta1_w - Mu_w * etam1_w==0;
	





    equations_endo = [fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8,fs9,fs10,fs11,fs12,fs13...
       ];
   
   equations_exo=[ err1,err2,err3,err4,err5,err6,err7];
    
    
     Contemp_endo = [mc,zcap,rk,k_s,q,c,i1,y,l,pinf,w,r,k]; ...
     Contemp_exo= [eps_a,eps_b,eps_g,eps_i,eps_r,eps_p,eps_w ];
           
     Expectation_endo=[mcp zcapp rkp k_sp   qp cp ip yp lp pinfp wp rp kp];...
    
     
     Lagged_endo=[mcm zcapm rkm k_sm   qm cm im ym lm pinfm wm rm km];...
    Lagged_exo=[   epsm_a  epsm_b epsm_g epsm_i  epsm_r  epsm_p epsm_w];
     
     Shocks = [eta1_a ,eta1_b ,eta1_g ,eta1_i ,eta1_r ,eta1_p ,eta1_w];
     Shocks_lagged=[etam1_a etam1_b etam1_g etam1_i etam1_r etam1_p etam1_w];
     
AA = equationsToMatrix(equations_endo, Contemp_endo);   
BB = equationsToMatrix(equations_endo, Lagged_endo);
CC = equationsToMatrix(equations_endo, Expectation_endo);
DD = equationsToMatrix(equations_endo, Contemp_exo);
EE = equationsToMatrix(equations_endo,Shocks);
RHO=equationsToMatrix(equations_exo,Lagged_exo);
FF=equationsToMatrix(equations_exo,Shocks);
GG=equationsToMatrix(equations_exo,Shocks_lagged);

AA=double(AA);
BB=double(BB);
CC=double(CC);
DD=double(DD);
EE=double(EE);
RHO=double(RHO);
FF=double(FF);
GG=double(GG);

    
AA=-AA;    
end