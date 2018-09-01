%function[aa_init cc_init dd_init rr_init]=simulation_MSV(parameters)
clear;clc;close all;
rng(1);
parameters=[-0.32 0.63 0.93 0.014 2.66 1.53 0.41 0.91 0.89 0.89 0.22 0.03 .32 0.03 0.01 0.02 0.13 0.008];    
gain=parameters(18);

numEndo=3;numExo=2;
param1=parameters(1:13);
param2=parameters(1:13);
param2(3)=parameters(14);
param2(13)=parameters(15);
param2(6)=0;param2(7)=0;param2(10)=0;
q_11=1-parameters(16);q_22=1-parameters(17);
N=5000;

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);


[AA1 BB1 CC1 DD1 EE1 FF1 rho1]= NKPC_sysmat_MSV(param1);
[AA2 BB2 CC2 DD2 EE2 FF2 rho2]= NKPC_sysmat_MSV(param2);

AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

X=zeros(5,N);X(:,1)=zeros(5,1);

eps_y(:,1) = normrnd(0,sigma_y1,[N,1]);
eps_y(:,2) = normrnd(0,sigma_y2,[N,1]);
eps_pinf(:,1) = normrnd(0,sigma_pinf1,[N,1]);    
eps_pinf(:,2) = normrnd(0,sigma_pinf2,[N,1]); 
eps_r(:,1) = normrnd(0,sigma_r1,[N,1]);
eps_r(:,2) = normrnd(0,sigma_r2,[N,1]);
errors1=[eps_y(:,1) eps_pinf(:,1) eps_r(:,1)]' ; 
errors2=[eps_y(:,2) eps_pinf(:,2) eps_r(:,2)]' ; 
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=0;

Q=[q_11,1-q_11;1-q_22,q_22];
ergodic_states=[(1-q_22)/(2-q_11-q_22);(1-q_11)/(2-q_11-q_22)];

alpha1=zeros(3,1);%constant. coef for learning
beta1=zeros(3,3);%coef. on lagged endo variables
cc1=rand(3,2);%coef on shocks
rr=10*eye(4);%auxiliary learning matrix



for tt=2:N
    
%      gain=1/tt;
gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(3)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(3)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;

gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho1];
gamma2_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_1_tilde;zeros(numExo,1)];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[EE1;FF1];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho2];
gamma2_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_2_tilde;zeros(numExo,1)];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[EE2;FF2];

regime(tt)=findRegime(regime(tt-1),q_11,q_22);


X(:,tt)=regime(tt)*(gamma2_1 + gamma1_1*X(:,tt-1)+gamma3_1*errors1(:,tt))+...
    (1-regime(tt))*(gamma2_2 + gamma1_2*X(:,tt-1)+gamma3_2*errors2(:,tt));


thetaOld=[alpha1 beta1(:,3) cc1];
[theta rr largestEig(tt) pr_flag(tt)] =msv_learning2(X(1:3,tt),[1;X(3,tt-1);X(4:5,tt)],thetaOld,rr,gain);
alpha1=theta(1,:)';beta1(:,3)=theta(2,:)';cc1=theta(3:4,:)';

alpha_tt(:,tt)=alpha1(1:3);
beta_tt(tt,:)=beta1(:,3);
cc_tt(tt,:,:)=cc1;
learningCovariance(tt,:,:)=rr;

end

% ll=round(N/10);
% aa_init=zeros(3,1);
% cc_init=mean(beta_tt(end-ll:end,:))';
% dd_init=squeeze(mean(cc_tt(end-ll:end,:,:)));
% rr_init=squeeze(mean(learningCovariance(end-ll:end,:,:)));

  %end

figure;
subplot(3,1,1);
plot(alpha_tt(1,:));
subplot(3,1,2);
plot(alpha_tt(2,:));
subplot(3,1,3);
plot(alpha_tt(3,:));

figure;
subplot(3,1,1);
plot(beta_tt(:,1));
subplot(3,1,2);
plot(beta_tt(:,2));
subplot(3,1,3);
plot(beta_tt(:,3));


figure;
subplot(3,2,1);
plot(squeeze(cc_tt(:,1,1)));
subplot(3,2,2);
plot(squeeze(cc_tt(:,1,2)));
subplot(3,2,3);
plot(squeeze(cc_tt(:,2,1)));
subplot(3,2,4);
plot(squeeze(cc_tt(:,2,2)));
subplot(3,2,5);
plot(squeeze(cc_tt(:,3,1)));
subplot(3,2,6);
plot(squeeze(cc_tt(:,3,2)));
