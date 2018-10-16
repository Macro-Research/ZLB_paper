function [likl] =kalmanNKPC(param,dataset)
y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
 rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
 eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 
% param=[0.8 0.7 0.9 0.1 5 1.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5];   
% param=[0.5475;0.6441;0.8556;0.0264;3.0583;1.4879;0.4660;0.3993;0.3255;0.8718;0.7289;0.2993;0.2944];
% load('full_dataset.mat');
first_obs=1;initLearn=1;burnIn=5;
% dataset=[gap_hp pinfobs robs];
% dataset=dataset(first_obs:end,:);
l=3;N=length(dataset);numVar=5;
%FL-variable indices [1 2]
%observable indices: y1 pinf2 robs3
%variable order in matrices: y pinf r u_y u_pi
[A B C D E F G]=NKPC_sysmat(param);
% gain=0.05;
beta1=diag([0.8642;0.7897;0;0;0]);
alpha1=0*ones(numVar,1);Ainv=A^(-1);
r_auxiliary=0*ones(numVar,1);
sigma=diag([eps_y^2;eps_pinf^2;eps_r^2]);
% sigmaInv=inv(sigma);
% load('init_var.mat');S(1,:)=init_state';P(1,:,:)=init_var;
% S(1,:)=[0.7810;0.0318;0.2243;0.6094;-0.0163];
% load('beta_previous.mat');
% 
% beta_previous=beta1;save beta_previous.mat beta_previous;
% beta1=diag([0.8726;0.8050;0.9084;0.3993;0.3255]);

% S(1,:)=zeros(numVar,1);
% S(1,:)=[1.3284;0.6759;1.08;0.6094;-0.0163];
% P(1,:,:)=1*eye(numVar);

% learning(1,1:5,2)=0.8762;learning(2,1:5,2)=0.8050;
likl=zeros(N,1);

% learning=nan(numVar,N,2);
% learning(:,initLearn,2)=diag(beta1);

% auxiliary_function=@(beta) function_g(beta,A,B,C,D,sigma);
% options=optimoptions('lsqnonlin',...
%     'MaxFunEvals',99999,'maxIter',9999,...
% 'tolFun',10e-5);
% % beta1=fsolve(auxiliary_function,0.1*ones(24,1),options);
% beta1=lsqnonlin(auxiliary_function,0.95*ones(numVar,1),...
%     zeros(numVar,1),0.9999*ones(numVar,1),options);
% beta1=diag(beta1);
% beta1=diag(beta1);
%    beta1=0.5*eye(numVar);
    gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;
 %   S(1,:)=zeros(numVar,1);
 
%      [ S(1,:),P(1,:,:)]=kalman_init(gamma1,gamma2,gamma3,sigma);
     S(1,:)=zeros(numVar,1);P(1,:,:)=1*eye(numVar);
%  beta1=NKPC_fixedPoint(A,B,C,D,sigma);

% P(1,:,:)=99*eye(numVar);
for i=2:N
    gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;
    Sprime=gamma1*S(i-1,:)'+gamma2;
    Pprime=gamma1*reshape(P(i-1,:,:),[numVar,numVar])*gamma1'+gamma3*sigma*gamma3';
     kGain=Pprime*F'*(F*Pprime*F')^(-1);
       v=(F*gamma3)^(-1)*(dataset(i,:)'-E-F*Sprime);%
    S(i,:)=Sprime+kGain*(dataset(i,:)'-E-F*Sprime);
%       S(i,3)=max(-r_bar,S(i,3));%ZLB constraint
    P(i,:,:)=Pprime-kGain*(F*Pprime);
  
%     likl(i)=-0.5*l*log(2*pi)-0.5*log(det(sigma))-0.5*v'*(sigmaInv*v);
likl(i)=-0.5*l*log(2*pi)-0.5*log(det(sigma))-0.5*v'*((sigma)\v);

% beta1=beta_update(A,B,C,D,sigma,beta1);
% learning(i-1,:)=diag(beta1);

%       
%       if i>1
%       for j=1:numVar
%           [ alpha1(j),beta1(j,j),r_auxiliary(j)] = recursive_update( S(1:i,j),i,alpha1(j),beta1(j,j),r_auxiliary(j));
% %          [ alpha1(j),beta1(j,j)] = learning_update( S(1:i,j),i);
% %          alpha1(j)=0;
% % alpha1(j)=0;
% % alpha1(j)=ywl_function_alpha(S(1:i,j),alpha1(j));
% % beta1(j,j)=ywl_function(S(1:i,j));
% %       [ alpha1(j),beta1(j,j)] = cg_learning( S(1:i,j),i);
% %         beta1(j,j)=ywl_function( S(1:i,j));
% % [alpha1(j),beta1(j,j)]=cgl_alpha_beta(S(1:i,j),alpha1(j));% beta1(j,j)=0;
%          learning(j,i,:)=[alpha1(j),beta1(j,j)];
% % 
%       end
%       end
%     beta1=beta_update(A,B,C,D,sigma,beta1);
%           learning(i,:)=diag(beta1);
end
i=1;
   v=(F*Ainv*D)^(-1)*(dataset(i,:)'-E);
likl(1)=-0.5*l*log(2*pi)-0.5*log(det(sigma))-0.5*v'*((sigma)\v);
 likl=-sum(likl(burnIn:N));

% disp(diag(beta1));
 
 


% 
end