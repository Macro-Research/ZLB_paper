clear;clc;close all;tic
path('LRE',path);
seed=round(1000*rand);rng(560)
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
N=1000;
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
% p_LN(1)=rand;p_NL(1)=rand;
Q=[1-p_LN,p_LN;p_NL,1-p_NL];
alphaTotal=zeros(numVar,1);
betaTotal=diag([0.8642;0.7897;0;0;0]);
%betaTotal=-0.1*diag(ones(numVar,1));
rTotal=zeros(numVar,1);
phi=2000;
shadow=X;
shadowRate=zeros(N,1);
learningMatrix=zeros(numVar,N,2);
for t=2:N
    disp(t)
  shadow(:,t)=A1_inv*( B1*X(:,t-1)+C1*(alphaTotal+betaTotal^2*(X(:,t-1)-alphaTotal))+D1*errors1(:,t));   
  shadowRate(t)=shadow(3,t);
%  if shadowRate(t)<-r_bar
%      regime(t)=1;
%  else regime(t)=0;
%  end
%  X(:,t)=shadow(:,t);
% p_LN(t)=1/(1+exp(-phi*(shadowRate(t)+r_bar))); 
% p_NL(t)=1/(1+exp(+phi*(shadowRate(t)+r_bar))); 
%   p_LN(mm,tt)=1/(1+exp(-PHI*(shadowRate(tt,mm)+r_bar))); 
%  p_NL(mm,tt)=1/(1+exp(+PHI*(shadowRate(tt,mm))+r_bar)); 
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
    
X(:,t) = (1-regime(t))* A1_inv * ( B1*X(:,t-1)+C1*(alphaTotal+betaTotal^2*(X(:,t-1)-alphaTotal))+D1*errors1(:,t))+...
(regime(t))* A2_inv * ( B2*X(:,t-1)+C2*(alphaTotal+betaTotal^2*(X(:,t-1)-alphaTotal))+D2*errors2(:,t));
    

for i=1:numVar
%  [alphaTotal(i),betaTotal(i,i),rTotal(i)]=...
%  recursive_update(X(1,:),t,alphaTotal(i),betaTotal(i,i),rTotal(i));
% %  [alphaTotal(i),betaTotal(i,i)]=...
% %  learning_update(X(i,1:t),t);
% % betaTotal(i,i)=ywl_function(X(i,1:t));
% % alphaTotal(i)=ywl_function_alpha(X(i,1:t),alphaTotal(i));
% 
learningMatrix(i,t,1)=alphaTotal(i);
learningMatrix(i,t,2)=betaTotal(i,i);
end

 



end
plotLength=500;
figure;

for i=1:numVar
    subplot(numVar,1,i)
    plot(X(i,end-plotLength:end),'lineWidth',2);
    title('Endogenous Variables');
    if i==3
      hold on;
      plot(shadowRate(end-plotLength:end),'--','lineWidth',2);
          legend('Realized Int.','Shadow Rate');
    end

end


figure;
plot(regime(end-plotLength+1:end),'lineWidth',3);
title('Realized Regimes','lineWidth',3);

figure;

subplot(2,1,1);
plot(p_LN(end-plotLength+1:end),'lineWidth',3);
subplot(2,1,2);
plot(p_NL(end-plotLength+1:end),'lineWidth',3);
title('Regime Probabilities','lineWidth',3);
learningMatrix(i,t,1)=alphaTotal(i);
learningMatrix(i,t,2)=betaTotal(i,i);
figure;
for i=1:2;
    for j=1:2
        subplot(2,2,(i-1)*2+j);
        plot(learningMatrix(i,end-plotLength+1:end,j),'lineWidth',3);
        title('learning parameters');
    end
end
disp('Average Regime Prob.')
[mean(regime) 1-mean(regime)]
seed