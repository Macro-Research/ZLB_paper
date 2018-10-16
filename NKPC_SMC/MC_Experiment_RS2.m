clear;clc;%close all;tic
%------------------SIMULATION
% seed=round(1000*rand);rng(560)
param1=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
param2=[0 0 0 0.01 3 0 0 0.5 0.5 0 0.7 0.3 0.03 ];
r_bar=param1(3);

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime2(param2);
N=100;
numVar=5;

A1_inv=A1^(-1);A2_inv=A2^(-1);

X=zeros(numVar,N);
X(:,1)=rand(numVar,1);
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
p_LN=0.05*ones(N,1);p_NL=0.05*ones(N,1); 
Q=[1-p_LN,p_LN;p_NL,1-p_NL];
alpha1=zeros(numVar,1);
% beta1=diag([0.8642;0.7897;0;0;0]);
beta1=0.5*eye(numVar);
phi=2000;
shadow=X;
shadowRate=zeros(N,1);
learningMatrix=zeros(numVar,N,2);
for t=2:N
    disp(t)
  shadow(:,t)=A1_inv*( B1*X(:,t-1)+C1*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D1*errors1(:,t));   
  shadowRate(t)=shadow(3,t);

 diff=shadowRate(t)-(-r_bar);
 if diff<0
     p_LN(t)=0;
     p_NL(t)=1;
 else 
     p_LN(t)=1;
     p_NL(t)=0;
 end
%  
 regime(t)=p_NL(t);

 
% regime(t)=findRegime(regime(t-1),p_LN(t),p_NL(t));
    
X(:,t) = (1-regime(t))* A1_inv * ( B1*X(:,t-1)+C1*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D1*errors1(:,t))+...
(regime(t))* A2_inv * ( B2*X(:,t-1)+C2*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D2*errors2(:,t));
    

end


M=5000;
dataset=X(1:3,:)';
p_LN=0.15*ones(M,N);p_NL=0.05*ones(M,N);
p_LN_upd=0.15*ones(M,N);p_NL_upd=0.05*ones(M,N);
r_bar=param1(3);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime2(param2);
A1_inv=A1^(-1);A2_inv=A2^(-1);


Sigma1=diag([sigma_y1^2;sigma_pinf1^2;sigma_r1^2]);


eps_y2=param2(end-2);eps_pinf2=param2(end-1);eps_r2=param2(end);
Sigma2=diag([sigma_y2^2;sigma_pinf2^2;sigma_r2^2]);


gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A1_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A1_inv*D2;
S0=zeros(numVar,1);P0=0.1*eye(numVar);

H=0.01*diag(var(dataset));


ne        = size(Sigma1,1);
[~, ns] = size(F1);
T         = size(dataset,1);
sqrtSigma    = gamma3_1*chol(Sigma1)';

% matrix for store
all_s_up  = zeros(T, ns, M);   % resampled 
lik       = zeros(T,1);
Neff      = zeros(T,1);


temp_s = S0;
temp_P = P0;
s_up   = repmat(temp_s, 1, M) + chol(temp_P)'*randn(ns, M);
weights=nan(T,M);
markov_state=nan(T,M); markov_state_upd=nan(T,M);
PHI=2000;
shadowRate1=nan(T,M);
resample=0;
for tt=1:1:T

    
    
    
  
end
