%Monte Carlo Experiment with ZLB regime: switching unscented kalman filter
clear;clc;close all;%tic
addpath('optimization_routines');
%------------------SIMULATION
seed=round(1000*rand);%rng(99)
param1=[0.2 0.35 0.72 0.06 2.48 1.32 0.45 0.7 0.24 0.84 0.78 0.25 0.31 ];
param2=[0.2 0.35 0.034 0.06 2.48 0    0 0.7 0.24 0 0.78 0.25 0.039 ];
r_bar=param1(3);

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);

[A1 B1 C1 D1, E1 F1 G1]=NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]=NKPC_sysmat_regime1(param2);
N=10000;
numVar=5;

A1_inv=A1^(-1);A2_inv=A2^(-1);

X=zeros(numVar,N);
X(:,1)=zeros(numVar,1);
eps_y(:,1) = normrnd(0,sigma_y1,[N,1]);
eps_y(:,2) = normrnd(0,sigma_y2,[N,1]);
eps_pinf(:,1) = normrnd(0,sigma_pinf1,[N,1]);    
eps_pinf(:,2) = normrnd(0,sigma_pinf2,[N,1]); 
eps_r(:,1) = normrnd(0,sigma_r1,[N,1]);
eps_r(:,2) = normrnd(0,sigma_r2,[N,1]);
errors1=[eps_y(:,1) eps_pinf(:,1) eps_r(:,1)]' ; 
errors2=[eps_y(:,2) eps_pinf(:,2) eps_r(:,2)]' ; 
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
p_11=0.98;p_22=0.5;
Q=[p_11,1-p_11;1-p_22,p_22];
alpha1=zeros(numVar,1);
% beta1=diag([0.8642;0.7897;0;0;0]);
beta1=0*eye(numVar);
% r_t=eye(numVar+1);
rr=zeros(numVar,1);
shadowRate=zeros(N,1);
learningMatrix=zeros(numVar,N,2);
for t=2:N
    gain=0.26;
    disp(t)
  shadow(:,t)=A1_inv*( B1*X(:,t-1)+C1*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D1*errors1(:,t));   
  shadowRate(t)=shadow(3,t);


  regime(t)=findRegime(regime(t-1),p_11,p_22);

    
X(:,t) = (regime(t))* A1_inv * ( B1*X(:,t-1)+C1*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D1*errors1(:,t))+...
(1-regime(t))* A2_inv * ( B2*X(:,t-1)+C2*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D2*errors2(:,t));


    for jj=1:numVar
   [alpha1(jj),beta1(jj,jj), rr(jj)]=l_SAC_CGL...
       (X(jj,t),X(jj,t-1),alpha1(jj),beta1(jj,jj),rr(jj),gain);
   alphaTotal(jj,t)=alpha1(jj);betaTotal(jj,jj,t)=beta1(jj,jj);
% [alpha1(jj) beta1(jj,jj)]=l_SAC_DGL(S_filtered(1:tt,jj),tt);
%     learning_filtered(jj,tt-1,:)=[alpha1(jj),beta1(jj,jj)];
    end
    
end

figure;
title('state variables');
for jj=1:numVar
    subplot(numVar,1,jj)
    plot(X(jj,:),'lineWidth',3);
end

figure;
title('alphas');
for jj=1:numVar
    subplot(numVar,1,jj)
   
   plot(alphaTotal(jj,:));
end

figure;
 title('betas')
 for jj=1:numVar
     for ii=1:numVar
 subplot(numVar,numVar,(jj-1)*numVar+ii);
 plot(squeeze(betaTotal(jj,ii,:)),'lineWidth',3);
 ylim([0 1]);
    end
 end