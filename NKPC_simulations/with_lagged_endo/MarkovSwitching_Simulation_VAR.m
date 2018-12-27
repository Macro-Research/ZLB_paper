%Simulation of 3-equation NKPC with Markov-switching, MSV-learning with
%least squares.
clear;clc;close all;%tic
addpath('c:/users/tolga/desktop/zlb_paper/optimization_routines');
%------------------SIMULATION
seed=round(1000*rand);%rng(99)
param1=[0.58 0.66 0.80 0.04 2.6414 1.5574 0.2889 0.3742 0.3256 0.9073 0.7396 0.2986 0.2991 ];
param2=[0.58 0.66 0.80 0.04 2.6414 0 0 0.3742 0.3256 0 0.7396 0.2986 0.02991 ];
r_bar=param1(3);

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat(param2);
N=2000;
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
p_11=.98;p_22=0.9; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
alpha1=zeros(numVar,1);
% beta1=diag([0.8642;0.7897;0;0;0]);
beta1=0.1*eye(numVar);
% r_t=eye(numVar+1);
 r1=10*eye(numVar+1);
% r1=zeros(2,2,2);
% r1(:,:,1)=eye(2);
% r1(:,:,2)=eye(2);
learningMatrix=zeros(numVar,N,2);

for t=2:N
    gain=0.04;
disp(t)

  regime(t)=findRegime(regime(t-1),p_11,p_22);

    
X(:,t) = (regime(t))* A1_inv * (( B1+C1*beta1^2)*X(:,t-1)+C1*(alpha1+beta1*alpha1)+D1*errors1(:,t))+...
    (1-regime(t))* A2_inv * (( B2+C2*beta1^2)*X(:,t-1)+C2*(alpha1+beta1*alpha1)+D2*errors2(:,t));


[alpha1(1:3),beta1(1:3,1:3),r1(1:4,1:4),largestEig(t),projectionFlag(t)]=...
    l_LS(X(1:3,t),[1;X(1:3,t-1)],alpha1(1:3),beta1(1:3,1:3),r1(1:4,1:4),gain);
   

alphaTotal(:,t)=alpha1;
betaTotal(:,:,t)=beta1;
rrTotal(:,:,t)=r1;

end

figure('Name','state variables');
title('state variables');
for jj=1:numVar
    subplot(numVar,1,jj)
    plot(X(jj,:),'lineWidth',3);
end

figure('Name','alphas');
title('alphas');
for jj=1:numVar
    subplot(numVar,1,jj)
   
   plot(alphaTotal(jj,:));
end

figure('Name','betas');
 title('betas')
 for jj=1:numVar
     for ii=1:numVar
 subplot(numVar,numVar,(jj-1)*numVar+ii);
 plot(squeeze(betaTotal(jj,ii,:)),'lineWidth',3);
%  ylim([0 1]);
    end
 end
 
 figure('Name','PLM-Covar Matrix');
  for jj=1:numVar+1
     for ii=1:numVar+1
 subplot(numVar+1,numVar+1,(jj-1)*(numVar+1)+ii);
 plot(squeeze(rrTotal(jj,ii,:)),'lineWidth',3);
%  ylim([0 1]);
    end
 end
 
 figure('Name','Projection Facility');
 plot(largestEig,'lineWidth',3);
 hold on;
 plot(projectionFlag,'lineWidth',3);
 legend('largest eigenvalue','projection facility');
 
 
 figure('Name','Regime');
 area(regime);
 title('regime');