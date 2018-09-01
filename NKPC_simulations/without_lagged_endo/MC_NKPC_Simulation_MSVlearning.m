clear;clc;close all;
%benchmark 3-equation model with MSV-least squares learning.
dynare NKPC_REE_Simulation;policyMatrix=oo_.dr.ghu(4:5,1:2)';
delete NKPC_REE_Simulation.log NKPC_REE_Simulation_results.mat;

clearvars -except policyMatrix;
tic
phi_y=0.5;phi_pinf=1.5;tau=3;gamma=0.01;lambda=0.99;
sigma_y=0.7;sigma_pinf=0.3;rho_pinf=0.5;rho_y=0.5;
A=[1+phi_y/tau,phi_pinf/tau;-gamma,1];
B=[1,1/tau;0,lambda];
C=[1,0;0,1];

N=10000;
numSimul=500;
transition=0;
numVar=4;numEndo=2;numExo=2;

RHO=diag([rho_y;rho_pinf]);

X=zeros(numEndo,N);
eps=zeros(2,N);

theta=zeros(3,2);
learning=zeros(numSimul,3,2);
largestEig=zeros(numSimul,N);
projectionFlag=zeros(numSimul,N);

 for mm=1:numSimul
    rng(mm); 
    
     eta_y = normrnd(0,sigma_y,[N,1]);
    eta_pinf = normrnd(0,sigma_pinf,[N,1]);    


errors=[eta_y eta_pinf]' ; 

A_inv = A^(-1);
aa_tt=zeros(numEndo,1);%constant. coef for learning
cc_tt=zeros(numEndo,numEndo);%coef. on lagged endo variables
dd_tt=rand(numEndo,numExo);%coef on shocks
rr_tt=10*eye(3);%auxiliary learning matrix
theta=zeros(numVar+1,numEndo);thetaOld=zeros(numEndo,numVar+1);

expectations=zeros(numEndo,N);
        
gain=0.01;
disp(mm);
toc

try
     for tt=2:N
eps(:,tt)=RHO*eps(:,tt-1)+errors(:,tt);%shocks realized

expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*eps(:,tt);


X(:,tt) =  A_inv*(B*expectations(:,tt)+C*eps(:,tt));%state realized

%expectations updated
thetaOld=[aa_tt dd_tt(1:2,1:2)];
[theta rr_tt largestEig(mm,tt) projectionFlag(mm,tt)] =l_LS_version2(X(:,tt),[1;eps(:,tt)],thetaOld,rr_tt,gain,[1 2 2]);
aa_tt=theta(1,:)';
%cc_tt=theta(2:3,:)';
dd_tt(1:2,1:2)=theta(2:3,:)';



     end
    
     learning(mm,:,:)=theta;
     
catch
     learning(mm,:,:)=nan(size(theta,1),size(theta,2));
     largestEig(mm,:)=nan;
    projectionFlag(mm,:)=ones(N,1);
end

 end
 
 save MC_MSV_results.mat learning learningEig projectionFlag policyMatrix;

 
