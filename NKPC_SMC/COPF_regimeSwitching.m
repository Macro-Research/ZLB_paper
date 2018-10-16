clear;clc;%close all;
% param=[0 0 0 0.01 3 1.5 0.5 0.25 0.25 0.9 0.7 0.3 0.3 ];
% param1=[.18 0.91 0.25 0.038 2.2373 1.5189 0.67017 0.5 0.5 0.90779 0.49265  0.12 0.16 ];
% param2=[.18 0.91 0 0.038 2.2373 0 0 0.5 0.5 0 0.49265  0.12 0.02 ];
param1=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
param2=[0 0 0 0.01 3 0 0 0.5 0.5 0 0.7 0.3 0.03 ];
M=10000;
load('simulated_dataset.mat');first_obs=1;initLearn=1;last_obs=length(gap_hp);
% load('simulated_data.mat');first_obs=1;initLearn=1;last_obs=length(gap_hp);
dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);l=3;N=length(dataset);numVar=5;
p_LN=0.15*ones(M,N);p_NL=0.05*ones(M,N);
r_bar=param1(3);

[A1 B1 C1 D1 E1 F1 G1]= NKPC_sysmat_regime1(param1);
[A2 B2 C2 D2 E2 F2 G2]= NKPC_sysmat_regime2(param2);
A1_inv=A1^(-1);A2_inv=A2^(-1);

eps_y1=param1(end-2);eps_pinf1=param1(end-1);eps_r1=param1(end);
Sigma1=diag([eps_y1^2;eps_pinf1^2;eps_r1^2]);


eps_y2=param2(end-2);eps_pinf2=param2(end-1);eps_r2=param2(end);
Sigma2=diag([eps_y2^2;eps_pinf2^2;eps_r2^2]);


 beta1=diag([0.8642;0.7897;0;0;0]);
alpha1=0*ones(numVar,1);Ainv=A1^(-1);
gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A1_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A1_inv*D2;
S0=zeros(numVar,1);P0=1*eye(numVar);
% H=diag([0.2943;0.1139;0.1805]);
H=0.1*diag(var(dataset));


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
markov_state=nan(T,M);
PHI=2000;
for tt=1:1:T
     disp(tt)
 shadowMatrix=gamma1_1*s_up;
 shadowRate(tt,:)=shadowMatrix(3,:)';
    for mm=1:M
%   p_LN(mm,tt)=1/(1+exp(-PHI*(shadowRate(tt,mm)+r_bar))); 
%  p_NL(mm,tt)=1/(1+exp(+PHI*(shadowRate(tt,mm))+r_bar)); 
 diff=shadowRate(tt,mm)-(-r_bar);
 if diff<0
     p_LN(mm,tt)=0;
     p_NL(mm,tt)=1;
 else 
     p_LN(mm,tt)=1;
     p_NL(mm,tt)=0;
 end
%  
 markov_state(tt,mm)=p_NL(mm,tt);
%  
%         if tt>1
%     markov_state(tt,mm)=findRegime(markov_state(tt-1,mm),p_LN(mm,tt),p_NL(mm,tt));
% % %  if shadowRate(tt,mm)<-r_bar
% % %      markov_state(tt,mm)=1;
% %  else markov_state(tt,mm)=0;
% % %  end
% % % 
%         else markov_state(tt,mm)=0;
%         end
    end
gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A2_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A2_inv*D2;

    yy = dataset(tt,:);

    markov_matrix1=repmat(markov_state(tt,:),[numVar 1]);
    markov_matrix2=repmat(markov_state(tt,:),[3 1]);
    %S_pred here is conditional the markov state for each particle
    %but why does p_pred not depend on p_{t-1}|{t-1}
    s_pred           =markov_matrix1.*(gamma1_2*s_up)+(ones(numVar,M)-markov_matrix1).*(gamma1_1*s_up);
    
    P_pred1           =gamma3_1*Sigma1*gamma3_1';
    P_pred1           = nearestSPD(P_pred1);
    
    P_pred2           =gamma3_2*Sigma2*gamma3_2';
    P_pred2           = nearestSPD(P_pred2);
    
    v1  = repmat(yy'-E1, 1, M) - F1*s_pred;
    v2  = repmat(yy'-E2, 1, M) - F2*s_pred;
    v=markov_matrix2.*v2+(ones(3,M)-markov_matrix2).*v1;
    Fe1  = F1*P_pred1*F1' + H;
    Fe2  = F2*P_pred2*F2' + H;
    S_upd1 = s_pred + P_pred1 * F1' * (Fe1\v);
    S_upd2 = s_pred + P_pred2 * F2' * (Fe2\v);
    S_upd=(ones(5,M)-markov_matrix1).*S_upd1+(markov_matrix1).*S_upd2;
  
    P_upd1 = P_pred1 - P_pred1 * F1' * (Fe1\F1) * P_pred1;
    P_upd2 = P_pred2 - P_pred2 * F1' * (Fe2\F1) * P_pred2;
    P_upd1 = nearestSPD(P_upd1);
    P_upd2 = nearestSPD(P_upd2);
    
    s_fore       =(ones(5,M)-markov_matrix1).*(S_upd + chol(P_upd1)' * randn(ns,M))+...
       (markov_matrix1).*(S_upd + chol(P_upd2)' * randn(ns,M)) ; 

    % Weights Calculation 
    perror  = (ones(3,M)-markov_matrix2).*(repmat(yy'-E1, 1, M) - F1*s_fore)+...
       (markov_matrix2).*(repmat(yy'-E2, 1, M) - F2*s_fore) ;
    adjustment        =(ones(1,M)-markov_state(tt,:))'.* (mvnpdf(s_fore',s_pred',P_pred1) ./ mvnpdf(s_fore',S_upd',P_upd1))+...
        (markov_state(tt,:))'.* (mvnpdf(s_fore',s_pred',P_pred2) ./ mvnpdf(s_fore',S_upd',P_upd2));
    
    %density=omega_tilde_t*W_{t-1}
    density        = mvnpdf(perror',zeros(1,size(yy,2)), H).* adjustment; 
    
    % Sum of density (normalized weights)
    weights(tt,:) = density/mean(density); %normalized weights
    
    
    
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
    
%        if tt>5
    for j=1:numVar
%      particle_states(tt,j)=mean(all_s_up(tt,j,:));
% %      [ alpha1(j),beta1(j,j)] = learning_update( particle_states(1:tt,j),tt); 
%     beta1(j,j)=ywl_function( particle_states(1:tt,j));
     learning(j,tt,:)=[alpha1(j),beta1(j,j)];
%       end
    end

  predictiveDensity(tt,:,:)=((ones(3,M)-markov_matrix2).*(repmat(E1, 1, M)+F1*s_pred)+...
      markov_matrix2.*(repmat(E2, 1, M)+F2*s_pred));
end
   lik=-sum(lik(25:end))
  
    for i=1:T for j=1:5
          particle_states(i,j)=mean(all_s_up(i,j,:));
          particle_states2(i,j)=median(all_s_up(i,j,:));
          particle_var(i,j)=var(all_s_up(i,j,:));
      end;end
  
  figure;
  grid minor
  for jj=1:numVar
subplot(5,1,jj);
plot(particle_states(:,jj),'lineWidth',3);
hold on;
plot(particle_states2(:,jj),'lineWidth',3);
  end

figure;

markov_state1=mean(markov_state');
p_LN1=mean(p_LN);
p_NL1=mean(p_NL);
subplot(3,1,1);
grid minor
plot(markov_state1,'lineWidth',3);
title('markov state')
subplot(3,1,2);
grid minor
plot(p_LN1,'lineWidth',3);
title('prob. from binding to not binding')
subplot(3,1,3);
grid minor
plot(p_NL1,'lineWidth',3);
title('prob. from not binding to binding');
%[mean(p_LN1) mean(p_NL1)]
[mean(markov_state1) 1-mean(markov_state1)]

figure;
grid minor
subplot(2,2,1);
plot(reshape(learning(1,:,1),[T 1]),'lineWidth',3);
title('alpha 1');
hold on;
subplot(2,2,2);
plot(reshape(learning(2,:,1),[T 1]),'lineWidth',3);
title('alpha 2');
subplot(2,2,3);
plot(reshape(learning(1,:,2),[T 1]),'lineWidth',3);
title('beta 1');
hold on;
subplot(2,2,4);
plot(reshape(learning(2,:,2),[T 1]),'lineWidth',3);
title('beta 2');

for tt=1:T
    for i=1:3
forecastMean(tt,i)=mean(predictiveDensity(tt,i,:));
forecastMedian(tt,i)=median(predictiveDensity(tt,i,:));
predictiveDensity(tt,i,:)=sort(predictiveDensity(tt,i,:));
aux=reshape(predictiveDensity(tt,i,:),[M 1]);
forecastUpper(tt,i)=aux(round(0.95*M));
forecastLower(tt,i)=aux(round(0.05*M));
    end
end

figure;
for jj=1:3
subplot(3,1,jj);
plot(forecastMean(:,jj),'*')
hold on;
plot(dataset(:,jj),'lineWidth',2);
hold on
plot(forecastUpper(:,jj),'--');
hold on
plot(forecastLower(:,jj),'--');
hold on
plot(forecastMedian(:,jj),':','lineWidth',2);
legend('mean Forecast','actual','upper','lower','median');
grid minor
end

figure;
subplot(3,2,1)
hist(reshape(predictiveDensity(15,3,:),[M 1]));
title('Predictive Density for Observable 1');
subplot(3,2,2)
hist(reshape(predictiveDensity(25,3,:),[M 1]));
title('Predictive Density for Observable 2');
subplot(3,2,3)
hist(reshape(all_s_up(15,3,:),[M 1]));
title('Predictive Density for Latent Variable 1');
subplot(3,2,4)
hist(reshape(all_s_up(25,3,:),[M 1]));
title('Predictive Density for Latent Variable 2');

figure;
index=1;
for jj=5:5:200
    subplot(8,5,index)
   index=index+1;
  hist(weights(jj,:),50)
end
