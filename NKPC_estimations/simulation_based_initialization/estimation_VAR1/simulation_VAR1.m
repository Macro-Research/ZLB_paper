 function[alpha_init,beta_init,rr_init]=simulation_VAR1(parameters)

%clear;clc;close all;
% load('estimation_results.mat');
% parameters=x;
%parameters=[0.1 0.1 1 0.03 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 0.03 0.03 0.05 0.15 0.01];    

rng(1);
gain=parameters(18);

param1=parameters(1:13);
param2=parameters(1:13);
param2(3)=parameters(14);
param2(13)=parameters(15);
param2(6)=0;param2(7)=0;param2(8)=0;
q_11=1-parameters(16);q_22=1-parameters(17);

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime1(param2);
N=5000;
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
regime(1)=0;
ergodic_states=[(1-q_22)/(2-q_11-q_22);(1-q_11)/(2-q_11-q_22)];
alpha1=zeros(numVar,1);
beta1=0*eye(numVar);
r1=eye(4);

for t=2:N
%gain=1/t;

  regime(t)=findRegime(regime(t-1),q_11,q_22);

    
X(:,t) = (regime(t))* A1_inv * (( B1+C1*beta1^2)*X(:,t-1)+C1*(alpha1+beta1*alpha1)+D1*errors1(:,t))+...
    (1-regime(t))* A2_inv * (( B2+C2*beta1^2)*X(:,t-1)+C2*(alpha1+beta1*alpha1)+D2*errors2(:,t));

      
    [alpha1(1:3),beta1(1:3,1:3),r1,largestEig(t),pr_flag(t)] =...
         msv_learning(X(1:3,t),[1;X(1:3,t-1)],...
      alpha1(1:3),beta1(1:3,1:3),r1,gain);

    

alphaTotal(:,t)=alpha1;
betaTotal(t,:,:)=beta1;
rrTotal(t,:,:,:)=r1;

end

ll=round(N/5);
alpha_init=zeros(5,1);
beta_init=squeeze(mean(betaTotal(end-ll:end,:,:)));
rr_init=squeeze(mean(rrTotal(end-ll:end,:,:,:)));
end
