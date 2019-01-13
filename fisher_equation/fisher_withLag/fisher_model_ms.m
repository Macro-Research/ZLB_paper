clear;clc;close all;
%Fisher equation with two (Markov) regimes, msv learning. 
%rng(4);
N=50000;burn_in=100;
alpha1=5;
alpha2=2;
rho=0.9;
iota_p=0.25;

eta_sigma=0.1;uu_sigma=0.001;
eta=normrnd(0,eta_sigma,[N 1]);
uu=normrnd(0,uu_sigma,[N 1]);
r=zeros(N,1);
pinf=zeros(N,1);

p_11=0.99;p_22=.99;
p_12=1-p_11;p_21=1-p_22;
Q=[p_11,1-p_11;1-p_22,p_22];
regime=zeros(N,1);
regime(1)=1;

learning=zeros(N,2);
learningCovariance=nan(N,2,2);

ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
alpha_erg=alpha1*ergodic_states(1)+alpha2*ergodic_states(2);

b_1=(alpha_erg-sqrt(alpha_erg^2-4*iota_p))/2
b_2=(alpha_erg+sqrt(alpha_erg^2-4*iota_p))/2%non-fundamental,E-unstable

a_1=2/(alpha_erg-2*rho+sqrt(alpha_erg^2-4*iota_p))
a_2=2/(alpha_erg-2*rho-sqrt(alpha_erg^2-4*iota_p))%non-fundamental,E-unstable

kappa=(b_1^2+iota_p)/(alpha1*alpha2);
PP=-(p_11*alpha2+p_22*alpha1)*kappa;
QQ=(p_11*p_22-p_12*p_21)*kappa;

MSS_eig1=-PP/2 + sqrt ( PP^2/4 - QQ);
MSS_eig2=-PP/2 - sqrt ( PP^2/4 - QQ);
MSS_spectral=max(abs(MSS_eig1),abs(MSS_eig2));

E_stabilityMarkov_1=(rho+b_1)/(alpha_erg);
E_stabilityMarkov_2=(2*b_1)/alpha_erg;
E_stabilityMarkov=max(E_stabilityMarkov_1,E_stabilityMarkov_2);

E_stabilityRegime1_1=(rho+b_1)/(alpha1);
E_stabilityRegime1_2=(2*b_1)/alpha1;
E_stabilityRegime1=max(E_stabilityRegime1_1,E_stabilityRegime1_2);

E_stabilityRegime2_1=(rho+b_1)/(alpha2);
E_stabilityRegime2_2=(2*b_1)/alpha2;
E_stabilityRegime2=max(E_stabilityRegime2_1,E_stabilityRegime2_2);

thr_ES=max(2*b_1,rho+b_1);
thr_MSS=MSS_spectral;

if E_stabilityMarkov<1
    disp('MARKOV CHAIN E-STABILITY SATISFIED');
else disp('MARKOV CHAIN E-STABILITY NOT SATISFIED');
end

if E_stabilityRegime1<1
    disp('REGIME 1 E-STABILITY SATISFIED');
else disp('REGIME 1 E-STABILITY NOT SATISFIED');
end

if E_stabilityRegime2<1
    disp('REGIME 2 E-STABILITY SATISFIED');
else disp('REGIME 2 E-STABILITY NOT SATISFIED');
end

if MSS_spectral<1
    disp('MEAN SQUARE STABILITY SATISFIED');
else disp('MEAN SQUARE STABILITY NOT SATISFIED');
end

aa=-1*rand(N,1);aa(1)=a_1;
bb=-1*rand(N,1);bb(1)=b_1;
rr=0.001*eye(2);

for jj=2:N
%     gain=1/jj;
   gain=0.01;
    regime(jj)=findRegime(regime(jj-1),p_11,p_22);
    r(jj)=rho*r(jj-1)+eta(jj);

    
    alpha_jj=regime(jj)*alpha1+(1-regime(jj))*alpha2;
    
pinf(jj)=(1/alpha_jj)*(aa(jj-1)*rho+bb(jj-1)*aa(jj-1)+1)*r(jj)+...
         (1/alpha_jj)*(bb(jj-1)^2+iota_p)*pinf(jj-1)+...
         uu(jj);
    

[aa(jj) bb(jj) rr largestEig(jj),projectionFac_flag(jj)] =...
    l_LS(pinf(jj),[r(jj) pinf(jj-1)]',aa(jj-1),bb(jj-1),rr,gain);
   
end

[mode(aa(burn_in:end)) mode(bb(burn_in:end))]
figure('Name','Fisher eqn-learning coef','units','normalized','outerposition',[0 0 1 1]);

subplot(3,1,1);
plot(aa,'lineWidth',3);
hold on;
plot(ones(N,1)*a_1,'--');
title('aa-cofficient on r_t');
subplot(3,1,2);
plot(bb,'lineWidth',3);
hold on;
plot(ones(N,1)*b_1,'--');
title('bb-coefficient on \pi_t');
subplot(3,1,3);
plot(pinf,'lineWidth',3);
title('Inflation');

figure('Name','Projection Facility','units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(largestEig,'lineWidth',3);
title('Largest Eigenvalue');
subplot(2,1,2);
plot(projectionFac_flag,'lineWidth',3);
title('Activity of Projection Facility');

figure('Name','distribution of learning parameters','units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1);
histogram(aa,20);
title('frequency dist of a_t: learning coef on r_t');
subplot(2,1,2);
histogram(bb,20);
title('frequency dist of b_t: learning coef on \pi_t');