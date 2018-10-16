% function [lik] = COPF_regimeSwitching2(param,dataset)
% y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
% gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
%  rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
%  eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 
param1=[.18 0.91 0.5 0.038 2.2373 1.5189 0.67017 0.73465 0.77657 0.90779 0.49265  0.12 0.16 ];
param2=[.18 0.91 0 0.038 2.2373 0 0 0.73465 0.77657 0.90779 0.49265  0.12 0.16 ];
M=1000;
load('full_dataset.mat');first_obs=195;initLearn=1;last_obs=215;
dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);
l=3;N=length(dataset);numVar=5;
p_LN=0.15*ones(M,N);p_NL=0.01*ones(M,N);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime2(param2);
A1_inv=A1^(-1);A2_inv=A2^(-1);
eps_y1=param1(end-2);eps_pinf1=param1(end-1);eps_r1=param1(end);
Sigma1=diag([eps_y1^2;eps_pinf1^2;eps_r1^2]);

eps_y2=param2(end-2);eps_pinf2=param2(end-1);eps_r2=param1(end);
Sigma2=diag([eps_y2^2;eps_pinf2^2;eps_r2^2]);

beta1=diag([0.8642;0.7897;0;0;0]);
alpha1=0*ones(numVar,1);Ainv=A1^(-1);
gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A1_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A1_inv*D2;
S0=zeros(numVar,1);P0=1*eye(numVar);
% H=diag([0.2943;0.1139;0.1805]);
H=0.1*diag(var(dataset));


ne        = size(Sigma,1);
[~, ns] = size(F1);
T         = size(dataset,1);
sqrtSigma    = gamma3_1*chol(Sigma)';

% matrix for store
all_s_up  = zeros(T, ns, M);   % resampled 
lik       = zeros(T,1);
Neff      = zeros(T,1);


temp_s = S0;
temp_P = P0;
s_up   = repmat(temp_s, 1, M) + chol(temp_P)'*randn(ns, M);
weights=nan(T,M);
markov_state=nan(T,M);
for tt=1:1:T
 shadowMatrix=gamma1_1*s_up;
 shadowRate(tt,:)=shadowMatrix(3,:)';
    for mm=1:M
  p_LN(mm,tt)=1/(1+exp(-2000*(shadowRate(tt,mm)-param(3)))); 
 p_NL(mm,tt)=1/(1+exp(+2000*(shadowRate(tt,mm)-param(3)))); 
 
        if tt>1
    markov_state(tt,mm)=findRegime(markov_state(tt-1,mm),p_LN(mm,tt),p_NL(mm,tt));
        else markov_state(tt,mm)=0;
        end
    end
    
    yy = dataset(tt,:);

    markov_matrix1=repmat(markov_state(tt,:),[numVar 1]);
    markov_matrix2=repmat(markov_state(tt,:),[3 1]);
    s_pred           =markov_matrix1.*(gamma1_2*s_up)+(ones(numVar,M)-markov_matrix1).*(gamma1_1*s_up);
    P_pred           =gamma3_1*Sigma*gamma3_1';
    P_pred           = nearestSPD(P_pred);
    
    v1  = repmat(yy'-E1, 1, M) - F1*s_pred;
    v2  = repmat(yy'-E2, 1, M) - F2*s_pred;
    v=markov_matrix2.*v2+(ones(3,M)-markov_matrix2).*v1;
    Fe  = F1*P_pred*F1' + H;
    s_upd = s_pred + P_pred * F1' * (Fe\v);
  
    P_upd = P_pred - P_pred * F1' * (Fe\F1) * P_pred;
    P_upd = nearestSPD(P_upd);
    
    s_fore       = s_upd + chol(P_upd)' * randn(ns,M); 

    % Weights Calculation 
    perror  = repmat(yy'-E1, 1, M) - F1*s_fore;
    adjustment        = mvnpdf(s_fore',s_pred',P_pred) ./ mvnpdf(s_fore',s_upd',P_upd);
    density        = mvnpdf(perror',zeros(1,size(yy,2)), H).* adjustment; 
    
    % Sum of density (normalized weights)
    weights(tt,:) = density/mean(density); 
    
    
    
   % Effective sample size
   Neff(tt,1) = (M^2)/sum(weights(tt,:).^2);
   
   
%    if Neff(tt,1)>=M/2
%        s_up = s_fore; 
%        
%    else
        if norm(weights(tt,:))>0
       s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');  % Resampling if ESS falls below a threshold
        else ind=1;
        end    
%    end
    % Store results
    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;

    
end
   lik=sum(lik(5:end));
% end
%   
%     for i=1:T for j=1:5
%           particle_states(i,j)=mean(all_s_up(i,j,:));
%           particle_var(i,j)=var(all_s_up(i,j,:));
%       end;end
%   
%   figure;
% subplot(5,1,1);
% plot(particle_states(:,1),'lineWidth',3);
% subplot(5,1,2);
% plot(particle_states(:,2),'lineWidth',3);
% subplot(5,1,3);
% plot(particle_states(:,3),'lineWidth',3);
% subplot(5,1,4);
% plot(particle_states(:,4),'lineWidth',3);
% subplot(5,1,5);
% plot(particle_states(:,5),'lineWidth',3);
% 
% figure;
% markov_state1=mean(markov_state');
% p_LN1=mean(p_LN);
% p_NL1=mean(p_NL);
% subplot(3,1,1);
% plot(markov_state1,'lineWidth',3);
% subplot(3,1,2);
% plot(p_LN1,'lineWidth',3);
% subplot(3,1,3);
% plot(p_NL1,'lineWidth',3);

