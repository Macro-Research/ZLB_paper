function [likl S] =kalmanNKPC_measError(param,dataset)
y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
 rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
 eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 
% clear;clc;%close all;
% param=[.03 0.82 1.19 0.0196 2.59 1.43 0.31 0.42 0.31 0.89 0.75 0.29 0.29 ];
load('full_dataset.mat');
first_obs=1;last_obs=length(gap_hp);
% dataset=[gap_hp pinfobs robs];
% dataset=dataset(first_obs:last_obs,:);
l=3;N=length(dataset);numVar=5;
%FL-variable indices [1 2]
%observable indices: y1 pinf2 robs3
%variable order in matrices: y pinf r u_y u_pi
[A B C D E F G]=NKPC_sysmat(param);


Sigma=diag([eps_y^2;eps_pinf^2;eps_r^2]);

beta1=diag([0.8642;0.7897;0;0;0]);
alpha1=0*ones(numVar,1);Ainv=A^(-1);

likl=zeros(N,1);

S(1,:)=zeros(numVar,1);P(1,:,:)=1*eye(numVar);
% H=diag([0.2943;0.1139;0.1805]);
H=0.25*diag(var(dataset));
gamma1=Ainv*(B+C*beta1^2);gamma2=Ainv*C*(eye(numVar)-beta1^2)*alpha1;gamma3=Ainv*D;
for i=2:N
    
    Sprime=gamma1*S(i-1,:)';
    
   Pprime=gamma1*reshape(P(i-1,:,:),[numVar,numVar])*gamma1'+gamma3*Sigma*gamma3';
   yprediction =E+F*Sprime;
  v=(dataset(i,:)'-yprediction);%
    kGain=Pprime*F'*(F*Pprime*F'+H)^(-1);
    S(i,:)=Sprime+kGain*(v);
    P(i,:,:)=Pprime-kGain*(F*Pprime);


likl(i)=-0.5*l*log(2*pi)-0.5*log(det(F*Pprime*F'+H))-0.5*v'*inv(F*Pprime*F'+H)*v;

%      if i>1
%       for j=1:numVar
%           [ alpha1(j),beta1(j,j),r_auxiliary(j)] = recursive_update( S(1:i,j),i,alpha1(j),beta1(j,j),r_auxiliary(j));
% %          [ alpha1(j),beta1(j,j)] = learning_update( S(1:i,j),i);
% %          alpha1(j)=0;
% % alpha1(j)=0;
% % alpha1(j)=ywl_function_alpha(S(1:i,j),alpha1(j));
% % beta1(j,j)=ywl_function(S(1:i,j));
% %       [ alpha1(j),beta1(j,j)] = cg_learning( S(1:i,j),i);
% %         beta1(j,j)=ywl_function( S(1:i,j));
% [alpha1(j),beta1(j,j)]=cgl_alpha_beta(S(1:i,j),alpha1(j));% beta1(j,j)=0;
% 
% % 
%       end
%       end



end
i=1;
   v=(F*Ainv*D)^(-1)*(dataset(1,:)'-E);
likl(1)=-0.5*l*log(2*pi)-0.5*log(det(F*Pprime*F'+H))-0.5*v'*inv(F*Pprime*F'+H)*v;
likl=-sum(likl(5:N));
end

