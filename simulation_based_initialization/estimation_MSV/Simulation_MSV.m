 function[aa_init cc_init dd_init rr_init]=simulation_MSV(parameters)
%clear;clc;%close all;
rng(1);
%parameters=[-0.28 0.65 0.86 0.0052 2.9 1.56 0.43 0.8964 0.8686 0.8923 0.1474 0.0386 .3174 0.0333 0.0102 0.0169 0.1291 0.0173];    
gain=parameters(18);
 
param1=parameters(1:13);
param2=parameters(1:13);
param2(3)=parameters(14);
param2(13)=parameters(15);
param2(6)=0;param2(7)=0;param2(10)=0;
q_11=1-parameters(16);q_22=1-parameters(17);
numEndo=3;numExo=3;N=10000;
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

Q=[q_11,1-q_11;1-q_22,q_22];
ergodic_states=[(1-q_22)/(2-q_11-q_22);(1-q_11)/(2-q_11-q_22)];

aa_tt=zeros(3,1);%constant. coef for learning
cc_tt=zeros(3,3);%coef. on lagged endo variables
dd_tt=rand(3,3);%coef on shocks
rr_tt=10*eye(4);%auxiliary learning matrix



for tt=2:N
    
% gain=1/tt;

EPS1(:,tt)=RHO*EPS1(:,tt-1)+errors1(:,tt);
EPS2(:,tt)=RHO*EPS2(:,tt-1)+errors2(:,tt);
regime(tt)=findRegime(regime(tt-1),q_11,q_22);
EPS(:,tt)=EPS1(:,tt)*regime(tt)+EPS2(:,tt)*(1-regime(tt));


expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*EPS(:,tt);

X(:,tt) =  regime(tt)*(A1_inv*(B1*X(:,tt-1)+C1*expectations(:,tt)+EPS1(:,tt)))+...;%state realized
       (1-regime(tt))*(A2_inv*(B2*X(:,tt-1)+C2*expectations(:,tt)+EPS2(:,tt)));

thetaOld=[aa_tt cc_tt(:,3) dd_tt(:,1:2)];
[theta rr_tt] =msv_learning2(X(:,tt),[1;X(3,tt-1);EPS(1:2,tt)],thetaOld,rr_tt,gain);
aa_tt=theta(1,:)';cc_tt(:,3)=theta(2,:)';dd_tt(:,1:2)=theta(3:4,:)';

aa(:,tt)=aa_tt;
cc(tt,:)=cc_tt(:,3);
dd(tt,:,:)=dd_tt(:,1:2);
learningCovariance(tt,:,:)=rr_tt;

end

ll=round(N/10);
aa_init=zeros(3,1);
cc_init=mean(cc(end-ll:end,:))';
dd_init=squeeze(mean(dd(end-ll:end,:,:)));
rr_init=squeeze(mean(learningCovariance(end-ll:end,:,:)));

 end

% end
% figure;
% subplot(3,1,1);
% plot(aa(1,:));
% subplot(3,1,2);
% plot(aa(2,:));
% subplot(3,1,3);
% plot(aa(3,:));
% 
% figure;
% subplot(3,1,1);
% plot(cc(:,1));
% subplot(3,1,2);
% plot(cc(:,2));
% subplot(3,1,3);
% plot(cc(:,3));
% 
% figure;
% subplot(3,2,1);
% plot(dd(:,1,1));
% subplot(3,2,2);
% plot(dd(:,1,2));
% subplot(3,2,3);
% plot(dd(:,2,1));
% subplot(3,2,4);
% plot(dd(:,2,2));
% subplot(3,2,5);
% plot(dd(:,3,1));
% subplot(3,2,6);
% plot(dd(:,3,2));
