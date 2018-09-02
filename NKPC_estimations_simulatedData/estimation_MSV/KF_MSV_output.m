% clear;clc;close all;
% parameters=[0 0 0 0.01 2 1.5 0.5 0.5 0.5 0.9 0.3 0.3 0.3 0 0.01 0.01 0.1 0.035];      
load('estimation_results.mat');
parameters=x;
gain=parameters(18);
burn_in=1;

load('simulated_dataset.mat');
dataset=simulated_dataset;


param(:,1)=parameters(1:13);
param(:,2)=parameters(1:13);
param(3,2)=parameters(14);

param(13,2)=parameters(15);

param(6,2)=0;param(7,2)=0;param(10,2)=0;
q_11=1-parameters(16);q_22=1-parameters(17);

Q=[q_11,1-q_11;1-q_22,q_22];
% q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
numVar=5;numExo=2;numObs=3;numRegimes=2;numEndo=3;
T=size(dataset,1);l=3;
alpha1=0*ones(3,1);beta1=0*eye(3);cc1=0*ones(3,2);rr=eye(4);

load('msv_initialBeliefs.mat');

alpha1=alpha_init;
beta1(:,3)=beta_init;
cc1=cc_init;
rr=rr_init;


[AA1 BB1 CC1 DD1 EE1 FF1 rho1 E1 F1 G1]=NKPC_sysmat_MSV(param(:,1));
[AA2 BB2 CC2 DD2 EE2 FF2 rho2 E2 F2 G2]=NKPC_sysmat_MSV(param(:,2));
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);


Sigma1=diag([param(end-2,1)^2;param(end-1,1)^2;param(end,1)^2]);
Sigma2=diag([param(end-2,2)^2;param(end-1,2)^2;param(end,2)^2]);
gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;


H1=zeros(3,3);H2=zeros(3,3);



S_fore11=zeros(numVar,1);S_fore12=zeros(numVar,1);S_fore21=zeros(numVar,1);S_fore22=zeros(numVar,1);
P_fore11=zeros(numVar,numVar);P_fore12=zeros(numVar,numVar);P_fore21=zeros(numVar,numVar);P_fore22=zeros(numVar,numVar);
v11=zeros(numObs,1);v12=zeros(numObs,1);v21=zeros(numObs,1);v22=zeros(numObs,1);
Fe11=zeros(numObs,numObs);Fe12=zeros(numObs,numObs);Fe21=zeros(numObs,numObs);Fe22=zeros(numObs,numObs);
S_upd11=zeros(numVar,1);S_upd12=zeros(numVar,1);S_upd21=zeros(numVar,1);S_upd22=zeros(numVar,1);
P_upd11=zeros(numVar,numVar);P_upd12=zeros(numVar,numVar);P_upd21=zeros(numVar,numVar);P_upd22=zeros(numVar,numVar);
ml11=0;ml12=0;ml21=0;ml22=0;
likl=zeros(T,1);

pp_fore11=1;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_upd11=1;pp_upd12=0;pp_upd21=0;pp_upd22=0;

pp_collapse1=0.5;pp_collapse2=0.5;
S_collapse1=zeros(numVar,1);S_collapse2=zeros(numVar,1);
P_collapse1=eye(numVar);P_collapse2=eye(numVar);
q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
S_filtered=zeros(T,numVar);

pp_filtered=ones(T,1);


%----------------------------------------------------------------------------

for tt=2:T
gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(3)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(3)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;

gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho1];
gamma2_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_1_tilde;zeros(numExo,1)];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[EE1;FF1];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho2];
gamma2_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_2_tilde;zeros(numExo,1)];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[EE2;FF2];


    x_tt=dataset(tt,:)';

%kalman block
S_fore11=gamma1_1*S_collapse1+gamma2_1;%
S_fore12=gamma1_2*S_collapse1+gamma2_2;%
S_fore21=gamma1_1*S_collapse2+gamma2_1;%
S_fore22=gamma1_2*S_collapse2+gamma2_2;%
%
P_fore11=gamma1_1*P_collapse1*gamma1_1'+gamma3_1*Sigma1*gamma3_1';
P_fore12=gamma1_2*P_collapse1*gamma1_2'+gamma3_2*Sigma2*gamma3_2';  
P_fore21=gamma1_1*P_collapse2*gamma1_1'+gamma3_1*Sigma1*gamma3_1';
P_fore22=gamma1_2*P_collapse2*gamma1_2'+gamma3_2*Sigma2*gamma3_2';
%    
v11=x_tt-E1-F1*S_fore11;
v12=x_tt-E2-F2*S_fore12;
v21=x_tt-E1-F1*S_fore21;
v22=x_tt-E2-F2*S_fore22;
%
Fe11=F1*P_fore11*F1'+H1;
Fe12=F2*P_fore12*F2'+H2;
Fe21=F1*P_fore21*F1'+H1;
Fe22=F2*P_fore22*F2'+H2;
%
S_upd11=S_fore11+P_fore11*F1'*(Fe11^(-1))*v11;
S_upd12=S_fore12+P_fore12*F2'*(Fe12^(-1))*v12;
S_upd21=S_fore21+P_fore21*F1'*(Fe21^(-1))*v21;
S_upd22=S_fore22+P_fore22*F2'*(Fe22^(-1))*v22;
%
P_upd11=(eye(numVar)-P_fore11*F1'*Fe11^(-1)*F1)*P_fore11;    
P_upd12=(eye(numVar)-P_fore12*F2'*Fe12^(-1)*F2)*P_fore12;  
P_upd21=(eye(numVar)-P_fore21*F1'*Fe21^(-1)*F1)*P_fore21;  
P_upd22=(eye(numVar)-P_fore22*F2'*Fe22^(-1)*F2)*P_fore22;  
%
pp_fore11=q_11*pp_collapse1;
pp_fore12=q_12*pp_collapse1;
pp_fore21=q_21*pp_collapse2;
pp_fore22=q_22*pp_collapse2;
%
ml11=(2*pi)^(-l/2)*det(Fe11)^(-0.5)*exp(-0.5*v11'*Fe11^(-1)*v11);
ml12=(2*pi)^(-l/2)*det(Fe12)^(-0.5)*exp(-0.5*v12'*Fe12^(-1)*v12);
ml21=(2*pi)^(-l/2)*det(Fe21)^(-0.5)*exp(-0.5*v21'*Fe21^(-1)*v21);
ml22=(2*pi)^(-l/2)*det(Fe22)^(-0.5)*exp(-0.5*v22'*Fe22^(-1)*v22);
% 
likl(tt)=ml11*pp_fore11+ml12*pp_fore12+ml21*pp_fore21+ml22*pp_fore22;
%
pp_upd11=(ml11*pp_fore11)/likl(tt);
pp_upd12=(ml12*pp_fore12)/likl(tt);
pp_upd21=(ml21*pp_fore21)/likl(tt);
pp_upd22=(ml22*pp_fore22)/likl(tt);
%
pp_collapse1=pp_upd11+pp_upd21;
pp_collapse2=pp_upd12+pp_upd22;
%
if pp_collapse1>10e-5;
    S_collapse1=(pp_upd11*S_upd11+pp_upd21*S_upd21)/pp_collapse1;
P_collapse1=(pp_upd11*(P_upd11+(S_collapse1-S_upd11)*(S_collapse1-S_upd11)')+...
    pp_upd21*(P_upd21+(S_collapse1-S_upd21)*(S_collapse1-S_upd21)'))/pp_collapse1;
else S_collapse1=zeros(numVar,1);
    P_collapse1=eye(numVar);
end

if pp_collapse2>10e-5
P_collapse2=(pp_upd12*(P_upd12+(S_collapse2-S_upd12)*(S_collapse2-S_upd12)')+...
    pp_upd22*(P_upd22+(S_collapse2-S_upd22)*(S_collapse2-S_upd22)'))/pp_collapse2;
S_collapse2=(pp_upd12*S_upd12+pp_upd22*S_upd22)/pp_collapse2;
else S_collapse2=zeros(numVar,1);
    P_collapse2=eye(numVar);
end

%
S_filtered(tt,:)=pp_collapse1*S_collapse1+pp_collapse2*S_collapse2;
% if pp_collapse1>pp_collapse2
%     S_filtered(tt,:)=S_collapse1;
% else S_filtered(tt,:)=S_collapse2;
% end

pp_filtered(tt)=pp_collapse1;




thetaOld=[alpha1 beta1(:,3) cc1];
[theta rr largestEig(tt) pr_flag(tt)] =msv_learning2(S_filtered(tt,1:numEndo)',[1,S_filtered(tt-1,3),S_filtered(tt,4:5)]',thetaOld,rr,gain);
alpha1=theta(1,:)';beta1(:,3)=theta(2,:)';cc1=theta(3:4,:)';

alpha_tt(:,tt)=alpha1(1:3);
beta_tt(:,tt)=beta1(1:3,3);
cc_tt(:,:,tt)=cc1;


error11(:,tt)=v11;
error12(:,tt)=v12;
error21(:,tt)=v21;
error22(:,tt)=v22;
errors_tot(:,tt)=error11(:,tt)*pp_upd11+error22(:,tt)*pp_upd22+error12(:,tt)*pp_upd12+error21(:,tt)*pp_upd21;

 end

likl=-sum(log(likl(burn_in+1:end)));



figure('Name','Interest Rate and Regime Probability');;
 plot(dataset(2:end,3),'--');
plot(ones(T-1,1)-pp_filtered(2:end),'lineWidth',3);
legend('true','filtered');


figure('Name','Regime Probability');
area(1-pp_filtered(2:end));
hold on;
plot(1-simulated_regime,'lineWidth',3);
legend('filtered','true');


figure('Name','Learning: Intercepts');
subplot(3,1,1);
plot(alpha_tt(1,2:end),'lineWidth',3);
hold on;
plot(alpha_simulated(1,2:end));
legend('filtered','true');
subplot(3,1,2);
plot(alpha_tt(2,2:end),'lineWidth',3);
hold on;
plot(alpha_simulated(2,2:end));
legend('filtered','true');
subplot(3,1,3);
plot(alpha_tt(3,2:end),'lineWidth',3);
hold on;
plot(alpha_simulated(3,2:end));
legend('filtered','true');
%------------------------------------

figure('Name','Learning: First-order autocorrelation');
subplot(3,1,1);
plot(squeeze(beta_tt(1,2:end)),'lineWidth',3);
hold on;
plot(squeeze(beta_simulated(1,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,1,2);
plot(squeeze(beta_tt(2,2:end)),'lineWidth',3);
hold on;
plot(squeeze(beta_simulated(2,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,1,3);
plot(squeeze(beta_tt(3,2:end)),'lineWidth',3);
hold on;
plot(squeeze(beta_simulated(3,2:end)),'lineWidth',3);
legend('filtered','true');
%-----------------------------------------

figure('Name','Learning: Shock coefficients');
subplot(3,2,1);
plot(squeeze(cc_tt(1,1,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(1,1,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,2,2);
plot(squeeze(cc_tt(1,2,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(1,2,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,2,3);
plot(squeeze(cc_tt(2,1,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(2,1,2:end)),'lineWidth',3);
legend('filtered','true');



subplot(3,2,4);
plot(squeeze(cc_tt(2,2,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(2,2,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,2,5);
plot(squeeze(cc_tt(3,1,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(3,1,2:end)),'lineWidth',3);
legend('filtered','true');

subplot(3,2,6);
plot(squeeze(cc_tt(3,2,2:end)),'lineWidth',3);
hold on;
plot(squeeze(cc_simulated(3,2,2:end)),'lineWidth',3);
legend('filtered','true');


figure('Name','model errors');
subplot(4,3,1);
plot(error11(1,:));
subplot(4,3,2);
plot(error11(2,:));
subplot(4,3,3);
plot(error11(3,:));

subplot(4,3,4);
plot(error12(1,:));
subplot(4,3,5);
plot(error12(2,:));
subplot(4,3,6);
plot(error12(3,:));

subplot(4,3,7);
plot(error21(1,:));
subplot(4,3,8);
plot(error21(2,:));
subplot(4,3,9);
plot(error21(3,:));

subplot(4,3,10);
plot(error22(1,:));
subplot(4,3,11);
plot(error22(2,:));
subplot(4,3,12);
plot(error22(3,:));


figure('Name','model errors');
subplot(3,1,1);
plot(errors_tot(1,:));
subplot(3,1,2);
plot(errors_tot(2,:));
subplot(3,1,3);
plot(errors_tot(3,:));
