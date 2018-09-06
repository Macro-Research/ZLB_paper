clear;clc;close all;

syms delta G phi_w curvp curvw phi sigma_c lambda xi_w sigma_l xi_p iota_w iota_p...
    psi phi_p r_pi rho r_y r_dy pi_bar beta1_const l_bar    gamma_bar alpha...
    rho_a rho_b  rho_g rho_i rho_r  rho_p rho_w mu_p mu_w rho_ga;


syms mc zcap rk k1   q c inve y lab pinf w r kp eps_a  eps_b eps_g eps_i  eps_r  eps_p eps_w dy dc dw dinve ...
    mcm zcapm rkm k1m   qm cm invem ym labm pinfm wm rm kpm eps_am  eps_bm eps_gm eps_im  eps_rm  eps_pm eps_wm dym dcm dwm dinvem...
 mcp zcapp rkp k1p   qp cp invep yp labp pinfp wp rp kpp eps_ap  eps_bp eps_gp eps_ip  eps_rp  eps_pp eps_wp dyp dcp dwp dinvep...
    eta_a eta_b eta_g eta_i eta_r eta_p eta_w...
    eta_am eta_bm eta_gm eta_im eta_rm eta_pm eta_wm;

NumEndo=17;
NumExo=7;
NumShock=7;


cpie=1+pi_bar/100;
gamma = 1 + gamma_bar/100 ;
beta1 = 1/(1 + beta1_const/100);
beta1_bar=beta1*gamma^(-sigma_c);
cr=cpie/(beta1*gamma^(-sigma_c));
r_bar=(cr-1)*100;
crk=(beta1^(-1))*(gamma^sigma_c) - (1-delta);
cw = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*crk^alpha))^(1/(1-alpha));
cikbar=(1-(1-delta)/gamma);
cik=(1-(1-delta)/gamma)*gamma;
clk=((1-alpha)/alpha)*(crk/cw);
cky=phi_p*(clk)^(alpha-1);
ciy=cik*cky;
ccy=1-G-cik*cky;
crkky=crk*cky;
cwhlc=(1/phi_w)*(1-alpha)/alpha*crk*cky/ccy;
cwly=1-crk*cky;




fs1=-mc +  alpha*rk+(1-alpha)*(w) - 1*eps_a - 0*(1-alpha)*eps_a==0 ;%ok...
	 fs2=-zcap +  (1/(psi/(1-psi)))* rk==0 ;%ok
	 fs3=-rk +  w+lab-k1 ==0;
	 fs4=-k1 +  kpm+zcap ==0;%ok
	fs5=-inve + (1/(1+beta1_bar*gamma))* (  invem + beta1_bar*gamma*invep+(1/(gamma^2*phi))*q ) +eps_i ==0;%ok
     fs6=-q + -r+pinfp-0*eps_b +(1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma))))*eps_b + (crk/(crk+(1-delta)))*rkp +  ((1-delta)/(crk+(1-delta)))*qp==0 ;%ok
fs7=-c + (lambda/gamma)/(1+lambda/gamma)*cm + (1/(1+lambda/gamma))*cp +((sigma_c-1)*cwhlc/(sigma_c*(1+lambda/gamma)))*(lab-labp) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma))*(r-pinfp + 0*eps_b) +eps_b ==0;%ok
fs8=-y + ccy*c+ciy*inve+eps_g  +  1*crkky*zcap ==0;%ok
fs9=-y + phi_p*( alpha*k1+(1-alpha)*lab +eps_a )==0;%ok
fs10=-pinf +  (1/(1+beta1_bar*gamma*iota_p)) * ( beta1_bar*gamma*pinfp +iota_p*pinfm...
               +((1-xi_p)*(1-beta1_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curvp+1)*(mc)  )  + eps_p==0;%ok 
fs11=-w +  (1/(1+beta1_bar*gamma))*wm...
               +(beta1_bar*gamma/(1+beta1_bar*gamma))*wp...
               +(iota_w/(1+beta1_bar*gamma))*pinfm...
               -(1+beta1_bar*gamma*iota_w)/(1+beta1_bar*gamma)*pinf...
               +(beta1_bar*gamma)/(1+beta1_bar*gamma)*pinfp...
               +(1-xi_w)*(1-beta1_bar*gamma*xi_w)/((1+beta1_bar*gamma)*xi_w)*(1/((phi_w-1)*curvw+1))*...
               (sigma_l*lab + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*cm -w) ...
               + 1*eps_w==0 ;
fs12=-r +  r_pi*(1-rho)*pinf...
               +r_y*(1-rho)*(y-phi_p*eps_a)  ...   
               +r_dy*(y-phi_p*eps_a-ym+phi_p*eps_am)...
               +rho*rm...
               +eps_r ==0 ;
fs13=-kp +  (1-cikbar)*kpm+cikbar*inve + cikbar*gamma^2*phi*eps_i==0 ;  %ok             
fs14= -dy + y - ym ==0;
fs15= -dc + c - cm ==0;
fs17= -dinve+inve-invem==0;
fs16= -dw + w - wm ==0;

err1=-eps_a + rho_a*eps_am  + eta_a==0;
err2=-eps_b + rho_b*eps_bm + eta_b==0;
err3=-eps_g + rho_g*(eps_gm) + eta_g + rho_ga*eta_a==0;
err4=-eps_i + rho_i*eps_im + eta_i==0;
err5=-eps_r + rho_r*eps_rm + eta_r==0;
err6=-eps_p + rho_p*eps_pm + eta_p - mu_p*eta_pm==0;
err7=-eps_w + rho_w*eps_wm + eta_w - mu_w*eta_wm==0 ;
	         


  
    F = [fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8,fs9,fs10,fs11,fs12,fs13 fs14 fs15 fs16 fs17...
        err1,err2,err3,err4,err5,err6,err7];
    
    
    Contemp = [mc,zcap,rk,k1,q,c,inve,y,lab,pinf,w,r,kp dy dc dinve dw...
               eps_a,eps_b,eps_g,eps_i,eps_r,eps_p,eps_w ];
           
     Expectation=[mcp zcapp rkp k1p   qp cp invep yp labp pinfp wp rp kpp dyp dcp dinvep dwp...
         eps_ap  eps_bp eps_gp eps_ip  eps_rp  eps_pp eps_wp];
     
     Lagged=[mcm zcapm rkm k1m   qm cm invem ym labm pinfm wm rm kpm dym dcm dinvem dwm...
         eps_am  eps_bm eps_gm eps_im  eps_rm  eps_pm eps_wm];
     
     Shocks = [eta_a ,eta_b ,eta_g ,eta_i ,eta_r ,eta_p ,eta_w];
     Shocks_lagged=[eta_am eta_bm eta_gm eta_im eta_rm eta_pm eta_wm];
         
    
AA = equationsToMatrix(F, Contemp);   
BB = equationsToMatrix(F, Lagged);
CC = equationsToMatrix(F, Expectation);
DD = equationsToMatrix(F, Shocks);
EE = equationsToMatrix(F,Shocks_lagged);
E =[gamma_bar;gamma_bar;gamma_bar;gamma_bar;pi_bar;r_bar;l_bar];
 F=[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
clearvars -except AA BB CC DD EE E F;
% AA=double(AA);
% BB=double(BB);
% CC=double(CC);
% DD=double(DD);
% EE=double(EE);
% RHO=double(RHO);
% FF=double(FF);
% GG=double(GG);


AA=-AA;    
