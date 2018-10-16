function [lik] = COPF_lik(param)
% clear;clc;close all;
% param=[.03 0.82 0.4 0.0196 2.59 1.43 0.31 0.42 0.31 0.89 0.75 0.5 0.5 ];
%param=[.18 0.91 .37 0.038 2.2373 1.5189 0.67017 0.73465 0.77657 0.90779 0.49265  0.12 0.16 ];
M=400; ind=0;
load('full_dataset.mat');first_obs=120;initLearn=1;last_obs=length(gap_hp);
dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);l=3;N=length(dataset);numVar=5;

[A B C D E F G]= NKPC_sysmat_regime1(param);%regime 1 is not binding
eps_y=param(end-2);eps_pinf=param(end-1);eps_r=param(end);
r_auxiliary=0*ones(numVar,1);
Sigma=diag([eps_y^2;eps_pinf^2;eps_r^2]);

 beta1=diag([0.8642;0.7897;0;0;0]);
%beta1=diag([0 0 0 0 0]);
alpha1=0*ones(numVar,1);Ainv=A^(-1);

S0=zeros(numVar,1);P0=1*eye(numVar);
% H=diag([0.2943;0.1139;0.1805]);
H=0.1*diag(var(dataset));

gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;x0=zeros(5,1);P0=1*eye(5);


ne        = size(Sigma,1);
[~, ns] = size(F);
T         = size(dataset,1);
sqrtSigma    = gamma3*chol(Sigma)';

% matrix for store
all_s_up  = zeros(T, ns, M);   % resampled 
lik       = zeros(T,1);
Neff      = zeros(T,1);


temp_s = S0;
temp_P = P0;
s_up   = repmat(temp_s, 1, M) + chol(temp_P)'*randn(ns, M);
P_upd=eye(numVar);
weights=nan(T,M);

for tt=1:1:T
    gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;

    yy = dataset(tt,:);
%     
  
    s_pred           = gamma1*s_up;
    P_pred           = gamma3*Sigma*gamma3';
    P_pred           = nearestSPD(P_pred);
    
    v  = repmat(yy'-E, 1, M) - F*s_pred;
    Fe  = F*P_pred*F' + H;
    s_upd = s_pred + P_pred * F' * (Fe\v);
  
    P_upd = P_pred - P_pred * F' * (Fe\F) * P_pred;
    P_upd = nearestSPD(P_upd);
    
    s_fore       = s_upd + chol(P_upd)' * randn(ns,M);

    % Weights Calculation 
    perror  = repmat(yy'-E, 1, M) - F*s_fore;
    adjustment        = mvnpdf(s_fore',s_pred',P_pred) ./ mvnpdf(s_fore',s_upd',P_upd);
    density        = mvnpdf(perror',zeros(1,size(yy,2)), H).* adjustment; 
    
    % Sum of density (normalized weights)
    weights(tt,:) = density/mean(density); %IS THIS W_T OR W_T^{TILDE}???
    
    
    
   % Effective sample size
   Neff(tt,1) = (M^2)/sum(weights(tt,:).^2); 
   
   
   if Neff(tt,1)>=M/2
       s_up = s_fore; 
       
   else
  
       s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');  % Resampling if ESS falls below a threshold
   
   end    
%    end
    % Store results
    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;
    
    if tt>1
    for j=1:numVar
     particle_states(tt,j)=mean(all_s_up(tt,j,:));
     %[ alpha1(j),beta1(j,j)] = learning_update( particle_states(1:tt,j),tt); 
     beta1(j,j)=ywl_function( particle_states(1:tt,j));
     learning(j,tt,:)=[alpha1(j),beta1(j,j)];
      end
    end
    
end
%   for i=1:T for j=1:5
%           particle_states(i,j)=mean(all_s_up(i,j,:));
%           particle_var(i,j)=var(all_s_up(i,j,:));
%       end;end
% if ind==1
%     lik=-Inf;
lik=-sum(lik(5:T,1));
%     for i=1:T for j=1:5
%           particle_states(i,j)=mean(all_s_up(i,j,:));
%           particle_var(i,j)=var(all_s_up(i,j,:));
%       end;end
  
%   figure;
% subplot(5,1,1);
% plot(particle_states(1:end,1),'lineWidth',3);
% subplot(5,1,2);
% plot(particle_states(1:end,2),'lineWidth',3);
% subplot(5,1,3);
% plot(particle_states(1:end,3),'lineWidth',3);
% subplot(5,1,4);
% plot(particle_states(1:end,4),'lineWidth',3);
% subplot(5,1,5);
% plot(particle_states(1:end,5),'lineWidth',3);
% 
% figure;
% subplot(2,2,1);
% plot(reshape(learning(1,:,1),[T 1]),'lineWidth',3);
% title('alpha 1');
% hold on;
% subplot(2,2,2);
% plot(reshape(learning(2,:,1),[T 1]),'lineWidth',3);
% title('alpha 2');
% subplot(2,2,3);
% plot(reshape(learning(1,:,2),[T 1]),'lineWidth',3);
% title('beta 1');
% hold on;
% subplot(2,2,4);
% plot(reshape(learning(2,:,2),[T 1]),'lineWidth',3);
% title('beta 2');

% end
    
% if ind==0
% lik=-sum(lik(5:T,1));
% else lik=Inf;
% end
% index=150;
% for i=1:25;
%     
%     subplot(6,5,i)
%     plot(weights(index,:));
%     index=1+index;
%     
% end
% figure;
% plot(lik(2:end),'lineWidth',3,'LineStyle',':');hold on;
% plot(likl(2:end),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% [sum(lik(2:end)) sum(likl(2:end))]
% 
% 

% figure;
% subplot(5,1,1);
% plot(particle_states(:,1),'lineWidth',3,'LineStyle',':');hold on;plot(S(:,1),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% subplot(5,1,2);
% plot(particle_states(:,2),'lineWidth',3,'LineStyle',':');hold on;plot(S(:,2),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% subplot(5,1,3);
% plot(particle_states(:,3),'lineWidth',3,'LineStyle',':');hold on;plot(S(:,3),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% subplot(5,1,4);
% plot(particle_states(:,4),'lineWidth',3,'LineStyle',':');hold on;plot(S(:,4),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% subplot(5,1,5);
% plot(particle_states(:,5),'lineWidth',3,'LineStyle',':');hold on;plot(S(:,5),'lineWidth',1,'LineStyle','-');legend('COPF','KF');
% 
% 
%     

% end
  