// copy van usmodel_hist_dsge_f19_7_71

var rm rp ams ebma  labobs robs pinfobs dy dc dinve dw  ewma epinfma  zcapf rkf kf pkf    cf invef yf labf wf rrf mc zcap rk k pk    c inve y lab pinf w r a  b g qs  ms  spinf sw kpf kp GS2 ;    
 
varexo ea eb eg  eqs  em  epinf ew  eam ;  
 
parameters conster2 cmab curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap  csadjcost ctou csigma chabb ccs cinvs cfc 
cindw cprobw cindp cprobp csigl clandaw 
crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
sig_a sig_b sig_g      sig_qs sig_pinf sig_w sig_am 
ctrend cg  
cdummy crp crhoams ;

parameters(zlb,2)  crpi crdy cry crr sig_m crams crb czlb 
parameters zlb_tp_1_2, zlb_tp_2_1

//parameters(volrs,2,"normal times","zlb constraint") flagrs 

model; 

// ? robs>=0.0;

# cpie=1+constepinf/100;
# cgamma=1+ctrend/100 ;
# cbeta=1/(1+constebeta/100);
# clandap=cfc;
# cbetabar=cbeta*cgamma^(-csigma);
# cr=cpie/(cbeta*cgamma^(-csigma));
# crk=(cbeta^(-1))*(cgamma^csigma) - (1-ctou);
# cw = (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
# cikbar=(1-(1-ctou)/cgamma);
# cik=(1-(1-ctou)/cgamma)*cgamma;
# clk=((1-calfa)/calfa)*(crk/cw);
# cky=cfc*(clk)^(calfa-1);
# ciy=cik*cky;
# ccy=1-cg-cik*cky;
# crkky=crk*cky;
# cwhlc=(1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
# cwly=1-crk*cky;
# conster=(cr-1)*100;

// flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-1*b+0*(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+1*b) + 0*b ;
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
//	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;
          kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(1+cbetabar*cgamma)*(cgamma^2*csadjcost)*qs ;

// sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -rp+pinf(1)-b  +0*(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rp-pinf(+1)+b ) +0*b ;
	      y = ccy*c+ciy*inve+g  +  1*crkky*zcap ;
	      y = cfc*( calfa*k+(1-calfa)*lab +a );
	      pinf =  (1/(1+cbetabar*cgamma*cindp)) * ( cbetabar*cgamma*pinf(1) +cindp*pinf(-1) 
               +((1-cprobp)*(1-cbetabar*cgamma*cprobp)/cprobp)/((cfc-1)*curvp+1)*(mc)  )  + spinf ; 
	      w =  (1/(1+cbetabar*cgamma))*w(-1)
               +(cbetabar*cgamma/(1+cbetabar*cgamma))*w(1)
               +(cindw/(1+cbetabar*cgamma))*pinf(-1)
               -(1+cbetabar*cgamma*cindw)/(1+cbetabar*cgamma)*pinf
               +(cbetabar*cgamma)/(1+cbetabar*cgamma)*pinf(1)
               +(1-cprobw)*(1-cbetabar*cgamma*cprobw)/((1+cbetabar*cgamma)*cprobw)*(1/((clandaw-1)*curvw+1))*
               (csigl*lab + (1/(1-chabb/cgamma))*c - ((chabb/cgamma)/(1-chabb/cgamma))*c(-1) -w) 
               + 1*sw ;
//	      r =  crpi*(1-crr)*pinf
//               +cry*(1-crr)*(y-yf)     
//               +crdy*(y-yf-y(-1)+yf(-1))
//               +crr*r(-1)
//               +ms + crams*ams + crb*sig_b*eb ;
	      r =  crpi*(1-crr)*pinf
               +cry*(1-crr)*(y-yf)     
               +crdy*(y-yf-y(-1)+yf(-1))
               +crr*r(-1)
               +ms + crams*ams - crb*sig_b*eb 
               -czlb*(1-crr)*(conster-0.25/4/2+0*conster) 
            ;

rm=r+ams;
rp=r+crp*ams;
	      a = crhoa*a(-1)  + sig_a*ea  ;
	      b = crhob*b(-1) + ebma - cmab*ebma(-1) ;
              ebma=sig_b*eb;
	      g = crhog*(g(-1)) + sig_g*eg + cgy*sig_a*ea;
	      qs = crhoqs*qs(-1) + sig_qs*eqs;
	      ms = crhoms*ms(-1) + sig_m*em;
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=sig_pinf*epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=sig_w*ew; 
          ams=crhoams*ams(-1)+sig_am*eam;

//	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*(1+cbetabar*cgamma)*cgamma^2*csadjcost*qs ;

// measurment equations

dy=y-y(-1)+ctrend;
dc=c-c(-1)+ctrend;
dinve=inve-inve(-1)+ctrend;
dw=w-w(-1)+ctrend;
pinfobs = 1*(pinf) + constepinf;
robs =    1*(r) + conster;
labobs = lab + constelab;
GS2 =4*(  1/8*(rm+rm(1)+rm(2)+rm(3)+rm(4)+rm(5)+rm(6)+rm(7))+conster+conster2 );

steady_state_model
    dy=ctrend;
    dc=ctrend;
    dinve=ctrend;
    dw=ctrend;
    pinfobs =constepinf;
    //robs =conster-czlb*(1-crr)*(conster-0.25/4/2+0*conster) ;
    labobs =constelab;
    GS2 = 4*(conster+conster2);
   //b=czlb*(1-crr)*(conster-0.25/4/2+0*conster) ;
   //r=-czlb*(1-crr)*(conster-0.25/4/2+0*conster) ;


observables dy dc dinve labobs pinfobs dw robs GS2;

parameterization
	ctou      ,  .025;
	clandaw   ,  1.5;
	cg        ,  0.18;
	curvp     ,  10;
	curvw     ,  10;

%flagrs(volrs,1), 0;
%flagrs(volrs,2), 1; 
%volrs_tp_1_2, 0;
%volrs_tp_2_1, 0;

crhoa     ,  .9842,  .85   ,  .1    ,  beta_pdf,0.01,0.9999 ;
crhob     ,  .8973,  .5    ,  .15   ,  beta_pdf,0.01,0.9999;
crhog     ,  .9727,  .85   ,  .1    ,  beta_pdf,0.01,0.9999;
crhoqs    ,  .6878,  .5    ,  .15   ,  beta_pdf,0.01,0.9999;
crhoms    ,  .2547,  .15   ,  .1    ,  beta_pdf,0.01,0.9999;
crhopinf  ,  .8267,  .85   ,  .1    ,  beta_pdf,0.01,0.9999;
crhow     ,  .9936,  .85   ,  .1    ,  beta_pdf,0.01,0.9999;  %old ones: .3,.7 
crhoams   ,  .8323,  .5    ,  .15   ,  beta_pdf,0.01,.999;
cmab      ,  .6770,  .5    ,  .15   ,  beta_pdf,0.01,0.9999;
cmap      ,  .7347,  .5    ,  .15   ,  beta_pdf,0.01,0.9999;
cmaw      ,  .9824,  .5    ,  .15   ,  beta_pdf,0.01,0.9999;

csadjcost ,  4.6731,  4    ,  1.0   ,  normal_pdf,1,15;
csigma    ,  1.1791,  1.5  ,  .25   ,  normal, 0.25,10  ;
chabb     ,  0.6222,  0.7  ,  0.1   ,  beta_pdf,0.001,0.99;
cprobw    ,  0.8383,  0.6  ,  0.1   ,  beta_pdf,0.1,0.95;%0.8202
csigl     ,  2.5038,  3.0  ,  0.75  ,  normal_pdf,0.25,10;
cprobp    ,  0.8836,  0.6  ,  0.1   ,  beta_pdf,0.1,0.95;
cindw     ,  0.4637,  0.5  ,  0.15  ,  beta_pdf,0.01,0.99;
cindp     ,  0.2191,  0.5  ,  0.15  ,  beta_pdf,0.01,0.99;
czcap     ,  0.6604,  0.5  ,  0.15  ,  beta_pdf,0.01,0.99;
cfc       ,  1.5222,  1.5  ,  0.25  ,  normal_pdf,1,3;

	crpi(zlb,1)    ,  2.0791,  1.5   ,  0.25,  normal_pdf,1,3;
	crr(zlb,1)     ,  0.8144,  0.75  ,  0.10,  beta_pdf,0.5,0.975;%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 0.83
	cry(zlb,1)     ,  0.1146,  0.125 ,  0.05,  normal_pdf,0.001,0.5;
    crdy(zlb,1)    ,  0.3249,  0.25  ,  0.10,  normal_pdf,0.001,0.5;
	crpi(zlb,2)    ,  .0;%1.01
	crr(zlb,2)     ,  .0;%0.999 or 1
	cry(zlb,2)     ,  .0;
    crdy(zlb,2)    ,  .0;
	
constepinf,  0.6454,  0.5  ,  0.05  ,  gamma_pdf,0.1,2;
constebeta,  0.1206,  0.25 ,  0.1   ,  gamma_pdf,0.01,2;
conster2  ,  0.1445,  0.2  ,  0.05  ,  gamma_pdf,0.01,2;
constelab ,  3.1744,  0.0  ,  2.5   ,  normal_pdf,-10,10;
ctrend    ,  0.4059,  0.4  ,  0.05  ,  normal_pdf,0.1,0.8;
cgy       ,  0.4846,  0.5  ,  0.25  ,  normal_pdf,0.01,2;
calfa     ,  0.1592,  0.3  ,  0.10  ,  normal_pdf,0.01,1;
//crams,0.0435      ,0.250,0.10 ,normal_pdf,-2.0,2.0;
//crb,0.5        ,0.250,0.10 ,normal_pdf,-2.0,2.0;
        crams(zlb,1),0.1874        ,0.250,0.1 ,normal_pdf,-2.0,2.0;
        crb(zlb,1),0.0395      ,0.250,0.1 ,normal_pdf,-2.0,2.0;
        crams(zlb,2),0;
        crb(zlb,2),0;
crp,0.3695        ,0.5,0.15 ,beta_pdf,0.001,0.999;

sig_a,   0.4615  ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_b,   0.8780  ,0.1,10,inv_gamma_pdf,0.01,3;
sig_g,   0.4792  ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_qs,  0.3814  ,0.1,10,inv_gamma_pdf,0.01,3;	
//	sig_m,   0.1648  ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_m(zlb,1),   0.2456 ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_m(zlb,2),   0.0077,0.025,.0125,uniform_pdf,0.00,0.05;
sig_pinf,0.1409  ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_w,   0.3607  ,0.1,10,inv_gamma_pdf,0.01,3;	
sig_am,  0.1857  ,0.10,10,inv_gamma_pdf,0.01,3;

    zlb_tp_1_2  	 ,0.0045   ,.10      ,.05,beta_pdf,0.001,0.999;	
	zlb_tp_2_1  	 ,0.3303   ,.40      ,.15,beta_pdf,0.001,0.999; 

czlb(zlb,1),0;
czlb(zlb,2),1;

