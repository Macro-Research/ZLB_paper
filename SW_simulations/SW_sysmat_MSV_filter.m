function [ AA BB CC DD EE RHO FF GG E F] = SW_sysmat_MSV_filter( parameters )

 delta   = parameters(1);
 G     = parameters(2);
 phi_w   =parameters(3); 
 curvp  = parameters(4);
 curvw  = parameters(5);
 
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
beta1_const=   parameters(21);
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
rho_ga=     parameters(32)  ;
mu_p=0;mu_w=0;


% cpie=1+(pi_bar/100);%ok
% gamma=1+(gamma_bar/100);%ok
% beta1=1/(1 + beta1_const/100);
% beta1_bar=beta1*gamma^(-sigma_c);
% crk = (beta1^(-1)) * (gamma^sigma_c) - (1-delta);
% cw = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*crk^alpha))^(1/(1-alpha));
% cikbar = (1-(1-delta)/gamma);
% cik = (1-(1-delta)/gamma)*gamma;
% clk = ((1-alpha)/alpha)*(crk/cw);
% cky = phi_p*(clk)^(alpha-1);
% ciy = cik * cky;
% ccy = 1 - G - cik*cky;
% crkky = crk*cky;
% cwhlc = (1/phi_w)*(1-alpha)/alpha*crk*cky/ccy;
% cwly=1-crk*cky;
% cr=((cpie/(beta1*gamma^(-sigma_c)))-1)*100;
% r_bar=cr;

cpie=1.005;
gamma = 1 + gamma_bar/100 ;
beta1 = 1/(1 + beta1_const/100);
beta1_bar=beta1*gamma^(-sigma_c);
cr=cpie/(beta1*gamma^(-sigma_c));
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

%gamma_bar=(gamma-1)*100;
r_bar=(cr-1)*100;
%pi_bar=(cpie-1)*100;
%l_bar=0;


%-----------------------------------------
AA=[1,0,-alpha,0,0,0,0,0,0,0,alpha-1,0,0,0,0,0,0;
0,1,(psi-1)/psi,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,1,1,0,0,0,0,-1,0,-1,0,0,0,0,0,0;
0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,-1/(phi*(gamma_bar/100+1)^2*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)),0,1,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0;
0,0,0,0,0,1,0,0,-(phi_p*(alpha-1)*(sigma_c-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1))/(alpha*phi_w*sigma_c*(lambda/(gamma_bar/100+1)+1)*(G+phi_p*((delta-1)/(gamma_bar/100+1)+1)*(gamma_bar/100+1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1)-1)),0,0,-(lambda/(gamma_bar/100+1)-1)/(sigma_c*(lambda/(gamma_bar/100+1)+1)),0,0,0,0,0;
0,-phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1),0,0,0,G+phi_p*((delta-1)/(gamma_bar/100+1)+1)*(gamma_bar/100+1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1)-1,-phi_p*((delta-1)/(gamma_bar/100+1)+1)*(gamma_bar/100+1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1),1,0,0,0,0,0,0,0,0,0;
0,0,0,-alpha*phi_p,0,0,0,1,phi_p*(alpha-1),0,0,0,0,0,0,0,0;
-(((xi_p*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)-1)*(xi_p-1))/(xi_p*((iota_p*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)*(curvp*(phi_p-1)+1)),0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
0,0,0,0,0,(((xi_w*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)-1)*(xi_w-1))/(xi_w*(lambda/(gamma_bar/100+1)-1)*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)*(curvw*(phi_w-1)+1)),0,0,-(sigma_l*((xi_w*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)-1)*(xi_w-1))/(xi_w*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)*(curvw*(phi_w-1)+1)),((iota_w*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)/((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1),(((xi_w*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)-1)*(xi_w-1))/(xi_w*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)*(curvw*(phi_w-1)+1))+1,0,0,0,0,0,0;
0,0,0,0,0,0,0,r_y*(rho-1)-r_dy,0,r_pi*(rho-1),0,1,0,0,0,0,0;
0,0,0,0,0,0,-(delta-1)/(gamma_bar/100+1)-1,0,0,0,0,0,1,0,0,0,0;
0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0;
0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0;
0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0;
0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,1];

BB=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;
0,0,0,0,0,0,1/((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1),0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,lambda/((gamma_bar/100+1)*(lambda/(gamma_bar/100+1)+1)),0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,iota_p/((iota_p*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1),0,0,0,0,0,0,0;
0,0,0,0,0,(lambda*((xi_w*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)-1)*(xi_w-1))/(xi_w*(gamma_bar/100+1)*(lambda/(gamma_bar/100+1)-1)*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)*(curvw*(phi_w-1)+1)),0,0,0,iota_w/((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1),1/((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1),0,0,0,0,0,0;
0,0,0,0,0,0,0,-r_dy,0,0,0,rho,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,-(delta-1)/(gamma_bar/100+1),0,0,0,0;
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0;
0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0];

CC=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,(gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)),0,0,0,0,0,0,0,0,0,0;
0,0,(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c),0,-(delta-1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c),0,0,0,0,1,0,0,0,0,0,0,0;
0,0,0,0,0,1/(lambda/(gamma_bar/100+1)+1),0,0,-(phi_p*(alpha-1)*(sigma_c-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1))/(alpha*phi_w*sigma_c*(lambda/(gamma_bar/100+1)+1)*(G+phi_p*((delta-1)/(gamma_bar/100+1)+1)*(gamma_bar/100+1)*(-((alpha-1)*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)*((alpha^alpha*(1-alpha)^(1-alpha))/(phi_p*(delta+(beta1_const/100+1)*(gamma_bar/100+1)^sigma_c-1)^alpha))^(1/(alpha-1)))/alpha)^(alpha-1)-1)),-(lambda/(gamma_bar/100+1)-1)/(sigma_c*(lambda/(gamma_bar/100+1)+1)),0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,(gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c*((iota_p*(gamma_bar/100+1))/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)),0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,(gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)),(gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c*((gamma_bar/100+1)/((beta1_const/100+1)*(gamma_bar/100+1)^sigma_c)+1)),0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

DD=[-1,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,1,0,0,0;
0,-(sigma_c*(lambda/(gamma_bar/100+1)+1))/(lambda/(gamma_bar/100+1)-1),0,0,0,0,0;
0,1,0,0,0,0,0;
0,0,1,0,0,0,0;
phi_p,0,0,0,0,0,0;
0,0,0,0,0,1,0;
0,0,0,0,0,0,1;
phi_p*r_y*(rho-1)-phi_p*r_dy,0,0,0,1,0,0;
0,0,0,phi*((delta-1)/(gamma_bar/100+1)+1)*(gamma_bar/100+1)^2,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0];

EE=[0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0];

FF=[1,0,0,0,0,0,0;
0,1,0,0,0,0,0;
rho_ga,0,1,0,0,0,0;
0,0,0,1,0,0,0;
0,0,0,0,1,0,0;
0,0,0,0,0,1,0;
0,0,0,0,0,0,1];

GG=[0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,0,0;
0,0,0,0,0,-mu_p,0;
0,0,0,0,0,0,-mu_w];

RHO=[rho_a,0,0,0,0,0,0;
0,rho_b,0,0,0,0,0;
0,0,rho_g,0,0,0,0;
0,0,0,rho_i,0,0,0;
0,0,0,0,rho_r,0,0;
0,0,0,0,0,rho_p,0;
0,0,0,0,0,0,rho_w];

E =[gamma_bar;gamma_bar;gamma_bar;gamma_bar;pi_bar;r_bar;l_bar];
 F=[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

end