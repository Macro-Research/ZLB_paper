endogenous ewma epinfma  zcapf rkf kf pkf cf invef yf labf wf rrf mc zcap rk k pk c inve y lab pinf w r a
b g qs ms  spinf sw kpf kp

endogenous labobs robs pinfobs dy dc dinve dw
 
exogenous ea eb eg  eqs  em  epinf ew   
 
parameters curvw cgy curvp constelab constepinf constebeta cmaw cmap calfa 
czcap csadjcost ctou csigma chabb cfc 
cindw cprobw cindp cprobp csigl clandaw  crpi crdy cry crr 
crhoa crhoas crhob crhog crhols crhoqs crhoms crhopinf crhow  
ctrend cg, sig_a, sig_b,sig_w,sig_pinf,sig_m,sig_qs,sig_g


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

% flexible economy

	      0*(1-calfa)*a + 1*a =  calfa*rkf+(1-calfa)*(wf)  ;
	      zcapf =  (1/(czcap/(1-czcap)))* rkf  ;
	      rkf =  (wf)+labf-kf ;
	      kf =  kpf(-1)+zcapf ;
	      invef = (1/(1+cbetabar*cgamma))* (  invef(-1) + cbetabar*cgamma*invef(1)+(1/(cgamma^2*csadjcost))*pkf ) +qs ;
          pkf = -rrf-0*b+(1/((1-chabb/cgamma)/(csigma*(1+chabb/cgamma))))*b +(crk/(crk+(1-ctou)))*rkf(1) +  ((1-ctou)/(crk+(1-ctou)))*pkf(1) ;
	      cf = (chabb/cgamma)/(1+chabb/cgamma)*cf(-1) + (1/(1+chabb/cgamma))*cf(+1) +((csigma-1)*cwhlc/(csigma*(1+chabb/cgamma)))*(labf-labf(+1)) - (1-chabb/cgamma)/(csigma*(1+chabb/cgamma))*(rrf+0*b) + b ;
	      yf = ccy*cf+ciy*invef+g  +  crkky*zcapf ;
	      yf = cfc*( calfa*kf+(1-calfa)*labf +a );
	      wf = csigl*labf 	+(1/(1-chabb/cgamma))*cf - (chabb/cgamma)/(1-chabb/cgamma)*cf(-1) ;
	      kpf =  (1-cikbar)*kpf(-1)+(cikbar)*invef + (cikbar)*(cgamma^2*csadjcost)*qs ;

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
	      r =  crpi*(1-crr)*pinf +cry*(1-crr)*(y-yf)+crdy*(y-yf-y(-1)+yf(-1))+crr*r(-1)+ms  ;
	      a = crhoa*a(-1)  + sig_a*ea;	%
	      b = crhob*b(-1) + sig_b*eb;
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
	crhoas,1; 	 % this parameter does not enter the model
	crhols,    0.9928;  % this parameter does not enter the model
	ctou,.025;
	clandaw,1.5;
	cg,0.18;
	curvp,10;
	curvw,10;
	
	% estimated parameters initialisation
	sig_a,   0.4618,0.1,2,inv_gamma_pdf;
	sig_b,   0.18185,0.1,2,inv_gamma_pdf;
	sig_g,   0.6090,0.1,2,inv_gamma_pdf;
	sig_qs,  0.46017,0.1,2,inv_gamma_pdf;
	sig_m,   0.2397,0.1,2,inv_gamma_pdf;
	sig_pinf,0.1455,0.1,2,inv_gamma_pdf;
	sig_w,   0.2089,0.1,2,inv_gamma_pdf;
	
	calfa,   0.24, 0.3,0.05,normal_pdf;
	csigma,   1.5,1.5,0.375,normal_pdf;
	cfc,      1.5, 1.5,0.125,normal_pdf;//phi_p
	cgy,      0.51,0.5,0.25,normal_pdf;//rho_ga
	
	csadjcost, 6.0144,4,1.5,normal_pdf;//phi
	chabb,     0.6361,0.7,0.1,beta_pdf;//lambda    
	cprobw,    0.8087,0.5,0.1,beta_pdf;//xi_w
	csigl,     1.9423,2,0.75,normal_pdf;//sigma_l
	cprobp,    0.6,0.5,0.1,beta_pdf;//xi_p
	cindw,     0.3243,0.5,0.15,beta_pdf;//iota_w
	cindp,     0.47,0.5,0.15,beta_pdf;//iota_p
	czcap,     0.2696,0.5,0.15,beta_pdf;//psi
	crpi,      1.488,1.5,0.25,normal_pdf;
	crr,       0.8762,0.75,0.1,beta_pdf;
	cry,       0.0593,0.125,0.05,normal_pdf;
	crdy,      0.2347,0.125,0.05,normal_pdf;
	crhoa,     0.9977,0.5,0.2,beta_pdf;
	crhob,     0.5799,0.5,0.2,beta_pdf;
	crhog,     0.9957,0.5,0.2,beta_pdf;
	crhoqs,    0.7165,0.5,0.2,beta_pdf;
	crhoms,    .3,0.5,0.2,beta_pdf;
	crhopinf,  0.8,0.5,0.2,beta_pdf;
	crhow,     0.8,0.5,0.2,beta_pdf;
	cmap ,     0.7,0.5,0.2,beta_pdf;
	cmaw  ,    0.7,0.5,0.2,beta_pdf;
	% derived from steady state
	constebeta, 0.25,.25,.01,gamma_pdf;
	
	ctrend,     0.3982,0.4,0.1,normal_pdf;
	constepinf, 0.7,0.625,0.1,gamma_pdf;
	constelab,  1.2918,0,2,normal_pdf;



