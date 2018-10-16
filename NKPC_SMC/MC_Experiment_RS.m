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
N=200;
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
p_LN=nan*ones(M,N);p_NL=nan*ones(M,N);
p_LN_upd=nan*ones(M,N);p_NL_upd=nan*ones(M,N);
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
     disp(tt)
 shadowMatrix=gamma1_1*s_up;
 shadowRate1(tt,:)=shadowMatrix(3,:)';
    for mm=1:M

 diff=shadowRate1(tt,mm)-(-r_bar);
 if diff<0
     p_LN(mm,tt)=0;
     p_NL(mm,tt)=1;
 else 
     p_LN(mm,tt)=1;
     p_NL(mm,tt)=0;
 end
 markov_state(tt,mm)=p_NL(mm,tt);
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
        resample=resample+1;
%         [id(resample,:), m] = systematic_resampling(weights(tt,:)');

%         s_up = squeeze(s_fore(:,id(resample,:)));
     s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');
%               weights(tt,:)=(1/M)*ones(M,1);
% s_up=multivariate_smooth_resampling(s_fore',weights(tt,:));
%          s_aux = squeeze(s_fore(:,id(resample,:)));
%       s_up = s_aux(:,randsample(M,M,true,weights(tt,:))');
         
  end
  
  for mm=1:M
  diff=s_up(3,mm)-(-r_bar);
 if diff<0
     p_LN_upd(mm,tt)=0;
     p_NL_upd(mm,tt)=1;
 else 
     p_LN_upd(mm,tt)=1;
     p_NL_upd(mm,tt)=0;
 end
 

  end
   markovState_upd(tt)=mean(p_NL_upd(:,tt));
  %UPDATE MARKOV PROBABILITIES
% s_aux=(weights(tt,:).*s_up(3,:))/sum(weights(tt,:));
 
  
  
  

    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;
   
    
    
%

%   predictiveDensity(tt,:,:)=((ones(3,M)-markov_matrix2).*(repmat(E1, 1, M)+F1*s_pred)+...
%       markov_matrix2.*(repmat(E2, 1, M)+F2*s_pred));
end
markov_statePrediction=mean(markov_state');
% markov_stateFiltered=mean(markov_state_upd');
markov_stateFiltered=markovState_upd;
figure
plot(regime,'lineWidth',3);
hold on;
plot(markov_stateFiltered,'lineWidth',3);
hold on;
plot(markov_statePrediction,'lineWidth',3);
legend('realized','filtered','predicted');
  figure;
  plot(Neff/M,'*');
  hold on;
  plot(0.5*ones(N,1),'lineWidth',3);
  title('Neff (as percentage of total # of particles)');

filteredShocks1=squeeze(all_s_up(:,4,:));
filteredShocks1=mean(filteredShocks1')';

filteredShocks2=squeeze(all_s_up(:,5,:));
filteredShocks2=mean(filteredShocks2')';
  
  figure;
  title('realized and filtered shocks');
  subplot(2,1,1);
  plot(X(4,:),'lineWidth',3);hold on;plot(filteredShocks1,'*');
  subplot(2,1,2);
  plot(X(5,:),'lineWidth',3);hold on;plot(filteredShocks2,'*');
  
  
difference=regime-markov_stateFiltered';
difference=round(difference);
correctPredictions=sum(difference(:) == 0);
falseFlags=sum(difference(:) == -1);
notPredicted=sum(difference(:) == 1);
tempText1=['CORRECT PREDICTION: ', num2str(100*correctPredictions/N),' %'];
tempText2=['FALSE FLAGS : ', num2str(100*falseFlags/N),' %'];
tempText3=['NOT PREDICTED : ', num2str(100*notPredicted/N),' %'];
tempText4=['RESAMPLE PERCENTAGE :', num2str(100*resample/N),'%'];
disp(tempText1);disp(tempText2);disp(tempText3);disp(tempText4);