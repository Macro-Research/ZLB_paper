%Simulation of 3-equation NKPC with Markov-switching, MSV-learning with
%least squares.
clear;clc;close all;tic
addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
%------------------SIMULATION
seed=round(1000*rand);
param1=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0.9 0.7 0.3 0.3 ];
param2=[0 0 0 0.01 3 0 0 0.5 0.5 0 0.7 0.3 0.03 ];
numEndo=3;numExo=3;N=10000;numSimul=500;
numVar=5;

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);


[A1 B1 C1 D1]= NKPC_sysmat(param1);
[A2 B2 C2 D2]= NKPC_sysmat(param2);
RHO=zeros(numExo,numExo);RHO(1:2,1:2)=B1(numEndo+1:end,numEndo+1:end);
A1=A1(1:numEndo,1:numEndo);B1=B1(1:numEndo,1:numEndo);C1=C1(1:numEndo,1:numEndo);


A2=A2(1:numEndo,1:numEndo);B2=B2(1:numEndo,1:numEndo);C2=C2(1:numEndo,1:numEndo);



A1_inv=A1^(-1);A2_inv=A2^(-1);
p_11=0.99;p_22=0.9; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
d_REE= (eye(size(C1,1)^2)-kron(RHO',(ergodic_states(1)*...
A1_inv*C1+ergodic_states(2)*A2_inv*C2)))^(-1)*vec(A1_inv*ergodic_states(1)+A2_inv*ergodic_states(2));
d_REE=reshape(d_REE,[size(C1,1),size(C1,2)]);
d_regime1= (eye(size(C1,1)^2)-kron(RHO',(ergodic_states(1)*...
A1_inv*C1)))^(-1)*vec(A1_inv);
d_regime1=reshape(d_regime1,[size(C1,1),size(C1,2)]);
learning=zeros(numSimul,6,3);
largestEig=zeros(numSimul,N);
projectionFlag=zeros(numSimul,N);
for mm=1:numSimul
rng(mm);

X=zeros(numEndo,N);X(:,1)=zeros(numEndo,1);
EPS1=zeros(numExo,N);EPS1(:,1)=zeros(numExo,1);
EPS2=zeros(numExo,N);EPS2(:,1)=zeros(numExo,1);
EPS=zeros(numExo,N);EPS(:,1)=zeros(numExo,1);
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


aa_tt=zeros(numEndo,1);%constant. coef for learning
cc_tt=zeros(numEndo,numEndo);%coef. on lagged endo variables
dd_tt=rand(numEndo,numExo);%coef on shocks
rr_tt=10*eye(6);%auxiliary learning matrix
expectations=zeros(numEndo,N);

disp(mm);
toc

% try
for tt=2:N
    gain=0.01;


EPS1(:,tt)=RHO*EPS1(:,tt-1)+errors1(:,tt);
EPS2(:,tt)=RHO*EPS2(:,tt-1)+errors2(:,tt);
regime(tt)=findRegime(regime(tt-1),p_11,p_22);
EPS(:,tt)=EPS1(:,tt)*regime(tt)+EPS2(:,tt)*(1-regime(tt));


expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*EPS(:,tt);

    
X(:,tt) =  regime(tt)*(A1_inv*(B1*X(:,tt-1)+C1*expectations(:,tt)+EPS1(:,tt)))+...;%state realized
       (1-regime(tt))*(A2_inv*(B2*X(:,tt-1)+C2*expectations(:,tt)+EPS2(:,tt)));

%thetaOld=[aa_tt cc_tt dd_tt];
thetaOld=[aa_tt cc_tt dd_tt(1:3,1:2)];
[theta rr_tt] =l_LS_version2(X(:,tt),[1;X(:,tt-1);EPS(1:2,tt)],thetaOld,rr_tt,gain,[1 3 3]);
aa_tt=theta(1,:)';cc_tt=theta(2:4,:)';dd_tt(1:3,1:2)=theta(5:6,:)';

% aa(:,tt)=aa_tt;
% cc(tt,:,:)=cc_tt;
% dd(tt,:,:)=dd_tt;
% learningCovariance(tt,:,:)=rr_tt;


%aa(:,tt)=aa_tt;
%cc(tt,:,:)=cc_tt;
%dd(tt,:,:)=dd_tt;
%learningCovariance(tt,:,:)=rr_tt;

end
    
    learning(mm,:,:)=theta;
      % largestEig(mm,:)=nan;
   % projectionFlag(mm,:)=ones(N,1);
    
% catch
%     learning(mm,:,:)=nan(3,4);
%    % largestEig(mm,:)=nan;
%    % projectionFlag(mm,:)=ones(N,1);
% end

end

save MC_MS_MSV_results.mat learning largestEig projectionFlag d_REE d_regime1;





