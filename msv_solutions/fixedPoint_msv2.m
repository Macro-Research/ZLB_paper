%clear;
%clc;

parameters=[0.58 0.66 0.80 0.0121 2.6414 1.5574 0.2889 0.3742 0.3256 0.9073 0.7396 0.2986 0.2991 ];


varCovar=[parameters(end-2)^2,0,0;0,parameters(end-1)^2,0;0,0,parameters(end)^2];
varCovar_vec=reshape(varCovar,[length(varCovar)^2,1]);
numVar=5;

[AA, BB, CC, DD EE FF GG]=NKPC_sysmat(parameters);
gamma1=AA^(-1)*BB;
gamma2=AA^(-1)*CC;
gamma3=AA^(-1)*DD;

func= @(beta)G_mapping(beta,gamma1,gamma2,gamma3,varCovar_vec,numVar);
func2=@(beta)G_mapping(beta,gamma1,gamma2,gamma3,varCovar_vec,numVar)-beta;
N=1000;
M=100;
gain=0.1;
beta11=zeros(numVar,numVar,N,M);
beta22=zeros(numVar,numVar,M);
options = optimset('Display','off');
% for jj=1:M
%     disp(jj)
% iteration=rand(numVar,numVar);
% iteration2=zeros(numVar,numVar);
% for ii=1:N
   
    % iteration=diag(iteration);
   % iteration=diag(diag(iteration));
%      iteration=iteration+gain*(func(iteration)-iteration);
%   iteration(iteration>1)=0;
%   iteration(iteration<-1)=0;
    
    
   
%       beta11(:,:,ii,jj)=iteration;
%       
% end
init=diag(0.85*ones(5,1));
 iteration2=lsqnonlin(func2,init,...
     -9999*ones(numVar,numVar),9999*ones(numVar,numVar),options)
abs( eigs(iteration2))
% beta22(:,:,jj)=iteration2;
% disp(iteration)
% display(iteration2)
% end
% disp(iteration)

% figure;

% for jj=1:2
%     for ii=1:5
% subplot(5,2,(jj-1)*5+ii);
% hist(squeeze(beta11(jj,ii,end,:)),50);
%     end
% end


% figure;
% 
% for jj=1:2
%     for ii=1:5
%         subplot(5,2,(jj-1)*5+ii);
%         hist(squeeze(beta22(ii,jj,:)),50);
%     end
% end
















%func=@(beta) gamma1+gamma2*diag(beta)^2-diag(beta);
 %options=optimoptions('lsqnonlin',...
  %  'MaxFunEvals',999999,'maxIter',999999,...
 %'tolFun',10e-6);


%fixedPoint1=lsqnonlin(func,rand(numVar,numVar),-9999*ones(numVar,numVar),9999*ones(numVar,numVar),options);
%fixedPoint2=fsolve(func,rand(numVar,1),options)
%disp(fixedPoint1');disp(fixedPoint2');




%fixedPoint=fsolve(func,(rand(5,1)))

 % beta1=fsolve(auxiliary_function,0.1*ones(numVar,1));
%beta1=lsqnonlin(auxiliary_function,0.95*ones(numVar,1),zeros(numVar,1),0.9999*ones(numVar,1),options);
% beta1=diag(beta1)

% MM= @(beta) gamma1+gamma2*beta^2;
% 
% vec0= @(beta) (eye(numVar^2)- kron(MM(beta),MM(beta)))^(-1)*...
%     (kron(gamma3,gamma3))*varCovar_vec;
% 
% vec1=@(beta) (kron(eye(numVar),gamma1)+kron(eye(numVar),gamma2*beta^2))*vec0(beta);
% 
% G = @(beta) reshape(vec1(beta)./vec0(beta),[numVar,numVar]);

%fixedPoint= @(beta) G(beta)-beta;

%fsolve(fixedPoint,diag(ones(5,1)));

 

