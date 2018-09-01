clear;clc;close all;
%benchmark 3-equation model with MSV-least squares learning.
dynare NKPC_REE_Simulation;policyMatrix=oo_.dr.ghu(4:5,1:2)';
delete NKPC_REE_Simulation.log NKPC_REE_Simulation_results.mat;

clearvars -except policyMatrix;
phi_y=0.5;phi_pinf=1.5;tau=3;gamma=0.01;lambda=0.99;
sigma_y=0.7;sigma_pinf=0.3;rho_pinf=0.9;rho_y=0.9;
A=[1+phi_y/tau,phi_pinf/tau;-gamma,1];
B=[1,1/tau;0,lambda];
C=[1,0;0,1];

N=2000;
transition=0;
numVar=4;numEndo=2;numExo=2;

RHO=diag([rho_y;rho_pinf]);

X=zeros(numEndo,N);
eps=zeros(2,N);



eta_y = normrnd(0,sigma_y,[N,1]);
eta_pinf = normrnd(0,sigma_pinf,[N,1]);    


errors=[eta_y eta_pinf]' ; 

A_inv = A^(-1);
largestEig=zeros(N,1);
aa_tt=zeros(numEndo,1);%constant. coef for learning
cc_tt=zeros(numEndo,numEndo);%coef. on lagged endo variables
dd_tt=zeros(numEndo,numExo);%coef on shocks
rr_tt=1*eye(numVar+1);%auxiliary learning matrix
theta=zeros(numVar+1,numEndo);thetaOld=zeros(numEndo,numVar+1);
learning=zeros(N,size(theta,1),size(theta,2));
expectations=zeros(numEndo,N);
    
for tt=2:N
gain=0.01;
%disp(tt);
     
eps(:,tt)=RHO*eps(:,tt-1)+errors(:,tt);%shocks realized

expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*eps(:,tt);


X(:,tt) =  A_inv*(B*expectations(:,tt)+C*eps(:,tt));%state realized

%expectations updated
thetaOld=[aa_tt cc_tt dd_tt];
[theta rr_tt] =l_LS_version2(X(:,tt),[1;X(:,tt-1);eps(:,tt)],thetaOld,rr_tt,gain,[1 2 2]);
aa_tt=theta(1,:)';cc_tt=theta(2:3,:)';dd_tt=theta(4:5,:)';

learning(tt,:,:)=theta;

 end

 
figure;
subplot(2,1,1);
plot(X(1,:),'lineWidth',3);
subplot(2,1,2);
plot(X(2,:),'lineWidth',3);

figure;
subplot(2,1,1);
plot(learning(:,1,1));
subplot(2,1,2);
plot(learning(:,1,2));
title('learning-intercepts');

% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'MS_simulation_intercept','-dpdf');

figure;
subplot(2,2,1);
plot(learning(:,2,1));
subplot(2,2,2);
plot(learning(:,2,2));
subplot(2,2,3);
plot(learning(:,3,1));
subplot(2,2,4);
plot(learning(:,3,2));
title('learning-lagged variables:these should be zeros');

% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'MS_simulation_betas','-dpdf');

figure;
subplot(2,2,1);
plot(learning(:,4,1));
hold on;plot(ones(N,1)*policyMatrix(1,1));
subplot(2,2,2);
plot(learning(:,4,2));
hold on;plot(ones(N,1)*policyMatrix(1,2));
subplot(2,2,3);
plot(learning(:,5,1));
hold on;plot(ones(N,1)*policyMatrix(2,1));
subplot(2,2,4);
plot(learning(:,5,2));
hold on;plot(ones(N,1)*policyMatrix(2,2));
title('learning-shocks');

% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'MS_simulation_shockCoef','-dpdf');