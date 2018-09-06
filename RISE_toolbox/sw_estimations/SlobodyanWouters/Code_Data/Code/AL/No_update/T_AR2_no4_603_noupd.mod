var labobs robs pinfobs dy dc dinve dw  ewma epinfma  mc zcap rk k pk    c inve y lab pinf w r a  b g qs  ms  spinf sw kp 
                cl invel labl pinfl pkl rkl wl;    
varexo  ea eb eg  eqs  em  epinf ew  ;  
parameters  curvw cgy curvp constelab constepinf constebeta   calfa 
            czcap cbeta csadjcost ctou csigma chabb ccs cinvs cfc 
            cindw cprobw cindp cprobp csigl clandaw 
            crdpi crpi crdy cry crr 
            crhoa crhoas crhob crhog crhols crhoqs crhoms   
            ctrend conster cg cgamma clandap cbetabar cr cpie crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly;
 
// fixed parameters

ctou=       0.025;
clandaw=    1.5;
cg=         0.18;
curvp=     10.0;
curvw=     10.0;

// estimated parameters

ctrend =    0.4;
constebeta= 0.25;
constepinf= 0.625;
constelab=  0.0;
calfa=      0.24;
cgy=        0.5;
csadjcost=  6.0144;
csigma=     1.5;
chabb=      0.6361;    
cprobw=     0.8087;
csigl=      1.9423;
cprobp=     0.6;
cindw=      0.3243;
cindp=      0.47;
czcap=      0.2696;
cfc=        1.5;
crpi=       1.488;
crr=        0.8762;
cry=        0.0593;
crdy=       0.2347;

crhoa=      0.9977;
crhob=      0.5799;
crhog=      0.9957;
crhols=     0.9928;
crhoqs=     0.7165;
crhoms=     0.0;
crhopinf=   0.5;
crhow=      0.5;
cmap =      0.5;
cmaw  =     0.5;

// derived parameters from steady state : see stst_f19.m

cpie=     1+constepinf/100;
cgamma=   1+ctrend/100 ;
cbeta=    1/(1+constebeta/100);
conster=  (cr-1)*100;
clandap=  cfc;
cbetabar= cbeta*cgamma^(-csigma);
cr=       cpie/(cbeta*cgamma^(-csigma));
crk=      (cbeta^(-1))*(cgamma^csigma) - (1-ctou);
cw =      (calfa^calfa*(1-calfa)^(1-calfa)/(clandap*crk^calfa))^(1/(1-calfa));
cikbar=   (1-(1-ctou)/cgamma);
cik=      (1-(1-ctou)/cgamma)*cgamma;
clk=      ((1-calfa)/calfa)*(crk/cw);
cky=      cfc*(clk)^(calfa-1);
ciy=      cik*cky;
ccy=      1-cg-cik*cky;
crkky=    crk*cky;
cwhlc=    (1/clandaw)*(1-calfa)/calfa*crk*cky/ccy;
cwly=     1-crk*cky;


model(linear); 
#stst_f19;


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
               +1*sw ;
            r = crpi*(1-crr)*pinf +cry*(1-crr)*(y-cfc*a) +crdy*(y-y(-1)-cfc*(a-a(-1))) +crr*r(-1) +ms  ;
            kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*cgamma^2*csadjcost*qs ;

            a     = crhoa*a(-1)   + ea;
            b     = crhob*b(-1)   + eb;
            g     = crhog*(g(-1)) + eg + cgy*ea;
            qs    = crhoqs*qs(-1)+ eqs;
            ms    = crhoms*ms(-1)+ em;
            spinf = epinfma;
                  epinfma=epinf;
            sw    = ewma;
                  ewma=ew; 

            dy      = y-y(-1)+ctrend;
            dc      = c-c(-1)+ctrend;
            dinve   = inve-inve(-1)+ctrend;
            dw      = w-w(-1)+ctrend;
            pinfobs = 1*(pinf) + constepinf;
            robs    = 1*(r) + conster;
            labobs  = lab + constelab;

            cl = c(-1);
            invel = inve(-1);
            labl = lab(-1);
            pinfl = pinf(-1);
            pkl = pk(-1);
            rkl = rk(-1);
            wl = w(-1);

            end; 

shocks;
var ea;
stderr 0.4618;
var eb;
stderr 1.8513;
var eg;
stderr 0.6090;
var eqs;
stderr 0.6017;
var em;
stderr 0.2397;
var epinf;
stderr 0.1455;
var ew;
stderr 0.2089;
end;


//stoch_simul(irf=20) y c inve pinf r w lab zcap df ;
//datatomfile('ddd',[]);

estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
stderr ea,0.4618,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eb,0.1818513,0.025,5,INV_GAMMA_PDF,0.1,2;
//stderr eb,1.818513,0.01,3,INV_GAMMA_PDF,.1,2;
stderr eg,0.6090,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eqs,0.46017,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr em,0.2397,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr epinf,0.1455,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr ew,0.2089,0.01,3,INV_GAMMA_PDF,0.1,2;
crhoa,.9676 ,.01,.9999,BETA_PDF,0.5,0.20;
crhob,.2703,.01,.9999,BETA_PDF,0.5,0.20;
crhog,.9930,.01,.9999,BETA_PDF,0.5,0.20;
crhoqs,.5724,.01,.9999,BETA_PDF,0.5,0.20;
crhoms,.3,.01,.9999,BETA_PDF,0.5,0.20;
//crhopinf,.8692,.01,.9999,BETA_PDF,0.5,0.20;
//crhow,.9546,.001,.9999,BETA_PDF,0.5,0.20;
//cmap,.7652,0.01,.9999,BETA_PDF,0.5,0.2;
//cmaw,.8936,0.01,.9999,BETA_PDF,0.5,0.2;
csadjcost,6.3325,2,15,NORMAL_PDF,4,1.5;
csigma,1.2312,0.25,3,NORMAL_PDF,1.50,0.375;
chabb,0.7205,0.001,0.99,BETA_PDF,0.7,0.1;
cprobw,0.7937,0.3,0.95,BETA_PDF,0.5,0.1;
//csigl,2.8401,0.25,10,NORMAL_PDF,2,0.75;
csigl,2.8401,0.25,10,NORMAL_PDF,2,0.5;
cprobp,0.7813,0.25,0.95,BETA_PDF,0.5,0.10;
cindw,0.4425,0.01,0.99,BETA_PDF,0.5,0.15;
cindp,0.3291,0.01,0.99,BETA_PDF,0.5,0.15;
czcap,0.2648,0.01,1,BETA_PDF,0.5,0.15;
cfc,1.4672,1.0,3,NORMAL_PDF,1.25,0.125;
crpi,1.7985,1.0,3,NORMAL_PDF,1.5,0.25;
crr,0.8258,0.5,0.975,BETA_PDF,0.75,0.10;
cry,0.0893,0.001,0.5,NORMAL_PDF,0.125,0.05;
crdy,0.2239,0.001,0.5,NORMAL_PDF,0.125,0.05;
constepinf,0.7,0.1,2.0,GAMMA_PDF,0.625,0.1;//20;
constebeta,0.7420,0.01,2.0,GAMMA_PDF,0.25,0.1;//0.20;
constelab,1.2918,-10.0,10.0,NORMAL_PDF,0.0,2.0;
ctrend,0.3982,0.1,0.8,NORMAL_PDF,0.4,0.10;
cgy,0.05,0.01,2.0,NORMAL_PDF,0.5,0.25;
calfa,0.24,0.01,1.0,NORMAL_PDF,0.3,0.05;
//gain,0.01,0.0001,0.07,GAMMA_PDF,0.04,0.03;
//ro,.9676,.01,1.,BETA_PDF,0.5,0.288675134594813;
//sigma,0.01,0,1,BETA_PDF,0.5,0.288675134594813;

end; 

//varobs y c inve lab pinf w r;

//observation_trends;
//y (ctrend);
//c (ctrend);
//inve (ctrend);
//w (ctrend);
//end; 
 
varobs dy dc dinve labobs pinfobs dw robs;

estimation(optim=('MaxIter',200),datafile=usmodel_data_SW_data_april2009,mode_compute=0,mode_file=T_AR2_no4_603_noupd_mode,first_obs=71,presample=4,nobs=176,kalman_algo=603,lik_init=2,prefilter=0,mh_replic=0,mh_nblocks=2,mh_jscale=0.2,mh_drop=0.2);

 