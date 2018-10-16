%----------------SIMULATION
clear;clc;close all;
tic
param=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.2 0.1 0.1 ];
sigma_y = param(end-2);
sigma_pinf=param(end-1);
sigma_r=param(end);
[A B C D E F G]=NKPC_sysmat(param);
N=1000;numVar=5;
X(:,1)=rand(5,1);
alpha1=zeros(5,1);
beta1=0.5*eye(numVar);
eps_y = normrnd(0,sigma_y,[N,1]);
eps_pinf = normrnd(0,sigma_pinf,[N,1]);    
eps_r = normrnd(0,sigma_r,[N,1]);
errors=[eps_y eps_pinf eps_r]' ; 
Ainv = A^(-1);
for t=2:N
X(:,t) =  Ainv * ( B*X(:,t-1)+...
    C*(alpha1+beta1^2*(X(:,t-1)-alpha1))+D*errors(:,t)) ;
end
%--------PARTICLE FILTER
M=400; ind=0;
dataset=X(1:3,:)';
l=3;N=length(dataset);
Sigma=diag([sigma_y^2;sigma_pinf^2;sigma_r^2]);
S0=zeros(numVar,1);P0=0.1*eye(numVar);
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
    
        weights(tt,:) = density/mean(density); 
    
    if sum(isnan(weights(tt,:)))>0
         warning('SOMETHING IS WRONG (POTENTIAL DEGENERACY)');
                 weights(tt,:)=1/M*ones(M,1);
       end
       
          % Effective sample size
   Neff(tt,1) = (M^2)/sum(weights(tt,:).^2);
   
   
   if Neff(tt,1)>=M/2
       s_up = s_fore; 
       
   else
        resample=resample+1;
%         s_up=residual_resampling(s_fore',weights(tt,:),rand(M,1)');
%         s_up=traditional_resampling(s_fore',weights(tt,:),rand(M,1)')';
%         
%         s_up=multivariate_smooth_resampling(s_fore',weights(tt,:));
        [id(resample,:), m] = systematic_resampling(weights(tt,:)');
%          weights(tt,:)=(1/M)*ones(M,1);
%         s_up = squeeze(s_fore(:,id(resample,:)));
%      s_up = s_fore(:,randsample(M,M,true,weights(tt,:))');
         s_aux = squeeze(s_fore(:,id(resample,:)));
      s_up = s_aux(:,randsample(M,M,true,weights(tt,:))');
   end
   
       % Store results
    lik(tt,1)        = log(mean(density));
    all_s_up(tt,:,:) = s_up;

    
end

particle_states=nan(T,numVar);
  lik=-sum(lik(5:end));
    for i=1:T for j=1:numVar
          particle_states(i,j)=mean(all_s_up(i,j,:));
      end;end

  for jj=4:5
      figure;
  plot(X(jj,:),'lineWidth',2);
  hold on;
  plot(particle_states(:,jj),'lineWidth',2,'lineStyle','--');
  legend('actual','filtered');
  end
   tempText=['RESAMPLE PERCENTAGE:', num2str(resample/N)];
  disp(tempText);
  figure;
  plot(Neff/M,'*');
  hold on;
  plot(0.5*ones(N,1),'lineWidth',3);
  title('Neff (as percentage of total # of particles)');
  
  [norm(abs(X(4,:)-particle_states(:,4)')) norm(abs(X(5,:)-particle_states(:,5)'))]
toc
  