endogenous ewma epinfma  mc zcap rk k pk c inve y lab pinf w r output_gap
 a b g qs ms  spinf sw kp

endogenous labobs robs pinfobs dy dc dinve dw
 
exogenous ea eb eg  eqs  em  epinf ew   
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap csadjcost ctou csigma chabb cfc 
cindw cprobw cindp cprobp csigl clandaw  
crhoa  crhob crhog crhoqs  crhopinf crhow  
ctrend cg, sig_a, sig_b,sig_w,sig_pinf,sig_qs,sig_g   crhoms

parameters(coef,2) crpi crr cry crdy sig_m czlb
parameters coef_tp_1_2, coef_tp_2_1 

steady_state_model

b=0;



model 

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





% sticky price - wage economy

	      mc =  calfa*rk+(1-calfa)*(w) - 1*a - 0*(1-calfa)*a ;
	      zcap =  (1/(czcap/(1-czcap)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+cbetabar*cgamma))* (  inve(-1) + cbetabar*cgamma*inve(1)+(1/(cgamma^2*csadjcost))*pk ) +qs ;
          pk = -r+pinf(1)-0*b +(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b + (crk/(crk+(1-ctou)))*rk(1) +  ((1-ctou)/(crk+(1-ctou)))*pk(1) ;
	      c = (chabb/cgamma)/(1+chabb/cgamma)*c(-1) + (1/(1+chabb/cgamma))*c(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(lab-lab(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(r-pinf(+1) + 0*b) +b ;
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
	      r =  crpi*(1-crr)*pinf +
                cry*(1-crr)*(output_gap)+
                crdy*(output_gap-output_gap(-1))+crr*r(-1)
                -czlb*(1-crr)*(conster-0.02)    
                + ms;
          output_gap=y-cfc*a;

	      a = crhoa*a(-1)  + sig_a*ea;	%
	      b = crhob*b(-1) + sig_b*eb -(1-crhob)*czlb*(1-crr)*(conster-0.02)   ;
	      g = crhog*(g(-1)) + sig_g*eg + cgy*sig_a*ea;
	      qs = crhoqs*qs(-1) + sig_qs*eqs;
	      ms = crhoms*ms(-1) + sig_m*em;
	      spinf = crhopinf*spinf(-1) + epinfma - cmap*epinfma(-1);
	          epinfma=sig_pinf*epinf;
	      sw = crhow*sw(-1) + ewma - cmaw*ewma(-1) ;
	          ewma=sig_w*ew; 
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

% measurment equations

		dy=y-y(-1)+ctrend;
		dc=c-c(-1)+ctrend;
		dinve=inve-inve(-1)+ctrend;
		dw=w-w(-1)+ctrend;
		pinfobs = 1*(pinf) + constepinf;
		robs =    1*(r) + conster;
		labobs = lab + constelab;

observables dy dc dinve dw pinfobs robs labobs


parameterization
	% fixed parameters

	ctou,.025;
	clandaw,1.5;
	cg,0.18;
	curvp,10;
	curvw,10;
	
	% estimated parameters initialisation
    %coef_tp_1_2,0.005;
    %coef_tp_2_1,0.3;
    coef_tp_1_2,0.03,0.1,0.05,beta_pdf;
    coef_tp_2_1,0.3,0.3,0.1,beta_pdf;
	sig_a,   0.4,0.1,2,inv_gamma_pdf;
	sig_b,   0.06,0.1,2,inv_gamma_pdf;
	sig_g,   0.37,0.1,2,inv_gamma_pdf;
	sig_qs,  0.28,0.1,2,inv_gamma_pdf;

	sig_pinf,0.10,0.1,2,inv_gamma_pdf;
	sig_w,   0.42,0.1,2,inv_gamma_pdf;
	
	calfa,   0.16, 0.3,0.05,normal_pdf;
	csigma,   1.1,1.5,0.375,normal_pdf;
	cfc,      1.5, 1.25,0.125,normal_pdf;//phi_p
	cgy,      0.66,0.5,0.25,normal_pdf;//rho_ga
	
	csadjcost, 4.5,4,1.5,normal_pdf;//phi
	chabb,     0.61,0.7,0.1,beta_pdf;//lambda    
	cprobw,    0.85,0.5,0.1,beta_pdf;//xi_w
	csigl,     2.5,2,0.75,normal_pdf;//sigma_l
	cprobp,    0.88,0.5,0.1,beta_pdf;//xi_p
	cindw,     0.48,0.5,0.15,beta_pdf;//iota_w
	cindp,     0.22,0.5,0.15,beta_pdf;//iota_p
	czcap,     0.83,0.5,0.15,beta_pdf;//psi



	crhoa,     0.98,0.5,0.2,beta_pdf;
	crhob,     0.2,0.5,0.2,beta_pdf;
	crhog,     0.98,0.5,0.2,beta_pdf;
	crhoqs,    0.68,0.5,0.2,beta_pdf;

	crhopinf,  0.15,0.5,0.2,beta_pdf;
	crhow,     0.05,0.5,0.2,beta_pdf;
	cmap ,     0;
	cmaw  ,    0;
	% derived from steady state
	constebeta, 0.25,.25,.01,gamma_pdf;
	
	ctrend,     0.4,0.4,0.1,normal_pdf;
	constepinf, 0.67,0.625,0.1,gamma_pdf;
	constelab,  1.5,0,2,normal_pdf;

	crpi(coef,1),      2,1.5,0.25,normal_pdf;
	crr(coef,1),       0.81,0.75,0.1,beta_pdf;
	cry(coef,1),       0.11,0.125,0.05,normal_pdf;
	crdy(coef,1),      0.3,0.125,0.05,normal_pdf;



	crpi(coef,2),      0;
	crr(coef,2),       0;
	cry(coef,2),       0;
	crdy(coef,2),      0;

	sig_m(coef,1),   0.1,0.1,2,inv_gamma_pdf;
    sig_m(coef,2),0.027,0.03,0.01,gamma_pdf;

    %conster,1.22,1,0.5,normal_pdf;
    %conster(coef,2),0.05,0.05,0.25,normal_pdf;

	crhoms,    .25,0.5,0.2,beta_pdf;
    %crhoms(coef,2),0;
    czlb(coef,1),0;
    czlb(coef,2),1;




steady_state_model
dy=ctrend;
dc=ctrend;
dinve=ctrend;
dw=ctrend;
pinfobs=constepinf;
robs=conster;
labobs=constelab;
