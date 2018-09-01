 clear;clc;close all;
% seed=round(1000*rand);
% rng(seed);
parameters=[0.58 0.66 0.80 0.0121 2.6414 1.5574 0.2889 0.3742 0.3256 0 0.7396 0.2986 0.2991 ];
N=2000;
transition=0;
numVar=5;

sigma_y = parameters(end-2);
sigma_pinf=parameters(end-1);
sigma_r=parameters(end);
[A B C D]= NKPC_matrixConverter(parameters);

X=zeros(numVar,N);

alphaTotal=zeros(numVar,N);
betaTotal=zeros(numVar,numVar,N);
rrTotal=zeros(numVar+1,numVar+1,N);
alpha_t=0.0*ones(numVar,1);
beta_t=0*ones(numVar,numVar);
r_t=0.5*eye(numVar+1);
       
eps_y = normrnd(0,sigma_y,[N,1]);
eps_pinf = normrnd(0,sigma_pinf,[N,1]);    
eps_r = normrnd(0,sigma_r,[N,1]);

errors=[eps_y eps_pinf eps_r]' ; 

A_inv = A^(-1);
largestEig=zeros(N,1);
warning('initialize')
 for tt=2:N

gain=0.01;
disp(tt);



      X(:,tt) =  A_inv * ( B*X(:,tt-1)+C*((eye(numVar)+beta_t)*alpha_t+beta_t^2*(X(:,tt-1)))+...
      D*errors(:,tt)     ) ;
 % largestEig(tt)=eigs(A^(-1)*(B+C*beta_t^2),1);

%       X(:,tt) =  A_inv * ( B*X(:,tt-1)+C*(alpha_t+beta_t*(X(:,tt-1)))+...
%        D*errors(:,tt)     ) ;


[alpha_t,beta_t,r_t,largestEig(tt),projectionFlag(tt)]=...
    l_LS(X(:,tt),[1;X(:,tt-1)],alpha_t,beta_t,r_t,gain);

beta_aux=beta_t;
beta_t=zeros(numVar,numVar);
beta_t(1:3,1:3)=beta_aux(1:3,1:3);


alphaTotal(:,tt)=alpha_t;
betaTotal(:,:,tt)=beta_t;
rrTotal(:,:,tt)=r_t;
    

[warnmsg,msgid]=lastwarn;
if strcmp(msgid,'MATLAB:singularMatrix')==1
    disp('BAD BAD MATRIX')
    break;
end
end
 
 figure;
 for jj=1:numVar
 subplot(5,1,jj);
 plot(alphaTotal(jj,transition+1:end),'lineWidth',3);
 hold on;
 plot(zeros(N-transition,1),'lineStyle','--','color','red');
 title('alphas');
 end
 
 %transition matrix taken from dynare:
 policyMatrix=[0,0,0,0.49,-0.97;0,0,0,0.009,0.97;0,0,0,0.26,0.97;0,0,0,0.5,0;0,0,0,0,0.5];
 figure;
 for jj=1:numVar
     for ii=1:numVar
 subplot(numVar,numVar,(jj-1)*numVar+ii);
 plot(squeeze(betaTotal(jj,ii,transition+1:end)),'lineWidth',3);
 hold on;
 plot(policyMatrix(jj,ii)*ones(N-transition,1),'lineStyle','--','color','red');
     end
     title('betas');
 end

%  figure;
%   for jj=1:numVar+1
%       for ii=1:numVar+1
%           subplot(numVar+1,numVar+1,(jj-1)*(numVar+1)+ii);
%           plot(squeeze(rrTotal(ii,jj,transition+1:end)),'lineWidth',3);
%           hold on;
%           plot(zeros(N-transition,1),'--','color','red');
%       end
%   end
 
 figure;
 subplot(2,1,1);
 plot(X(1,transition+1:end));
 subplot(2,1,2);
 plot(X(2,transition+1:end));
 
 
   
 disp(squeeze(betaTotal(:,:,end)))
figure;
plot(largestEig,'lineWidth',3);
hold on;
plot(projectionFlag,'lineWidth',3);
legend('largest eigenvalue','projection flag');