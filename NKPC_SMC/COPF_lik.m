function [lik] = COPF_lik(param,dataset)
% clear;clc;close all;
%   param=[.5474    0.6441      0.8557      0.0264     3.0583       1.4879... 
%    0.4660     0.3993      0.3255      0.8718     0.3     0.1     0.1 ];
% param=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
 rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
 eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 
M=500; ind=0;
% load('full_dataset.mat');
%  first_obs=120;initLearn=1;last_obs=length(gap_hp);
% dataset=[gap_hp pinfobs robs];
% dataset=dataset(first_obs:last_obs,:);
l=3;N=length(dataset);numVar=5;
[A B C D E F G]=NKPC_sysmat(param);



Sigma=diag([eps_y^2;eps_pinf^2;eps_r^2]);

beta1=diag([0.8642;0.7897;0;0;0]);
% beta1=0.5*eye(numVar);
alpha1=0*ones(numVar,1);Ainv=A^(-1);

S0=zeros(numVar,1);P0=0.1*eye(numVar);
% H=diag([0.2943;0.1139;0.1805]);
H=0.05*diag(var(dataset));

gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;


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
weights=nan(T,M);
resample=0;

for tt=1:1:T
    
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
    weights(tt,:) = density/mean(density); 
    
       if sum(isnan(weights(tt,:)))>0
        
         warning('SOMETHING IS WRONG (POTENTIAL DEGENERACY)');
%                  weights(tt,:)=1/M*ones(M,1);
lik=-Inf;
break;
       end

    
    
   % Effective sample size
   Neff(tt,1) = (M^2)/sum(weights(tt,:).^2);
   
   
   if Neff(tt,1)>=M/2
       s_up = s_fore; 
       
   else
%         if norm(weights(tt,:))>0
        [id, m] = systematic_resampling(weights(tt,:)');
        s_up = squeeze(s_fore(:,id));
     weights(tt,:)=(1/M)*ones(M,1);
%        s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');
% weights(tt,:)=(1/M)*ones(M,1);
%        s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');
       resample=resample+1;% Resampling if ESS falls below a threshold
%         else ind=1;
%         end    
   end
    % Store results
    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;

    
end
%   for i=1:T for j=1:5
%           particle_states(i,j)=mean(all_s_up(i,j,:));
%           particle_var(i,j)=var(all_s_up(i,j,:));
%       end;end
% if ind==1
%     lik=-Inf;
  lik=-sum(lik(5:end));
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

end
  