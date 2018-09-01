 clear;clc;close all;
 tic
load('estimation_results.mat');
%parameters=x;
%clear;clc;close all;
parameters=[-0.3 0.64 0.86 0.004 2.6 1.55 0.43 0.9 0.87 0.89 0.15 0.04 0.32 0.03 0.01 0.02 0.13 0.0173];    
   gain=parameters(18);
%  q_11=0.99;q_22=0.85; 

load('us_dataset.mat');
first_obs=24;burn_in=20;
dataset=[gap_cbo,pinfobs,robs];
dataset=dataset(first_obs:end,:);
startDate=datenum('01-01-1966');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,length(dataset(2:end,3)));

param(:,1)=parameters(1:13);
param(:,2)=parameters(1:13);
param(3,2)=parameters(14);
% param(12,1)=0.3;param(13,1)=0.3;
param(13,2)=parameters(15);
% param(13,2)=0.03;
param(6,2)=0;param(7,2)=0;param(10,2)=0;
q_11=1-parameters(16);q_22=1-parameters(17);

Q=[q_11,1-q_11;1-q_22,q_22];
% q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
numVar=5;numExo=2;numObs=3;numRegimes=2;numEndo=3;
T=size(dataset,1);l=3;
alpha1=0*ones(3,1);beta1=0*eye(3);cc1=0*ones(3,2);rr=eye(4);
 [alpha_init,beta_init,cc_init,rr_init]=Simulation_MSV2(parameters);
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


H1=0*diag(var(dataset(1:150,:)));H2=0*diag(var(dataset(150:end,:)));
%H2(1,1)=0;H2(2,2)=0;



S_fore11=zeros(numVar,1);S_fore12=zeros(numVar,1);S_fore21=zeros(numVar,1);S_fore22=zeros(numVar,1);
P_fore11=zeros(numVar,numVar);P_fore12=zeros(numVar,numVar);P_fore21=zeros(numVar,numVar);P_fore22=zeros(numVar,numVar);
v11=zeros(numObs,1);v12=zeros(numObs,1);v21=zeros(numObs,1);v22=zeros(numObs,1);
Fe11=zeros(numObs,numObs);Fe12=zeros(numObs,numObs);Fe21=zeros(numObs,numObs);Fe22=zeros(numObs,numObs);
S_upd11=zeros(numVar,1);S_upd12=zeros(numVar,1);S_upd21=zeros(numVar,1);S_upd22=zeros(numVar,1);
P_upd11=zeros(numVar,numVar);P_upd12=zeros(numVar,numVar);P_upd21=zeros(numVar,numVar);P_upd22=zeros(numVar,numVar);
ml11=0;ml12=0;ml21=0;ml22=0;
likl=zeros(T,1);
%pp_fore11=0;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_fore11=1;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_upd11=1;pp_upd12=0;pp_upd21=0;pp_upd22=0;
%pp_collapse1=regime(1);pp_collapse2=1-regime(1);
pp_collapse1=0.5;pp_collapse2=0.5;
S_collapse1=zeros(numVar,1);S_collapse2=zeros(numVar,1);
P_collapse1=eye(numVar);P_collapse2=eye(numVar);
q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
S_filtered=zeros(T,numVar);
%pp_filtered=zeros(T,1);
pp_filtered=ones(T,1);


%----------------------------------------------------------------------------

for tt=2:T
%    gain=1/tt;
gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;

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
pp_filtered(tt)=pp_collapse1;

thetaOld=[alpha1 beta1(:,3) cc1];
[theta rr largestEig(tt) pr_flag(tt)] =msv_learning2(S_filtered(tt,1:numEndo)',[1,S_filtered(tt-1,3),S_filtered(tt,4:5)]',thetaOld,rr,gain);
alpha1=theta(1,:)';beta1(:,3)=theta(2,:)';cc1=theta(3:4,:)';

alpha_tt(:,tt)=alpha1(1:3);
beta_tt(:,tt)=beta1(:,3);
cc_tt(:,:,tt)=cc1;


 end






likl=-sum(log(likl(burn_in+1:end)));

toc

figure('Name','Interest Rate and Regime Probability');;
 plot(Date,dataset(2:end,3),'--');
hold on;
plot(Date,ones(T-1,1)-pp_filtered(2:end),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_filter_init_MSV_regime','-dpdf');

figure('Name','Regime Probability');
area(Date,1-pp_filtered(2:end));
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_filter_init_MSV_regimeProb','-dpdf');

figure('Name','Learning: Intercepts');
subplot(3,1,1);
plot(Date,alpha_tt(1,2:end),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
subplot(3,1,2);
plot(Date,alpha_tt(2,2:end),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
subplot(3,1,3);
plot(Date,alpha_tt(2,2:end),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_filter_init_MSV_alphas','-dpdf');

figure('Name','Learning: First-order autocorrelation');
subplot(3,1,1);
plot(Date,squeeze(beta_tt(1,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');

subplot(3,1,2);
plot(Date,squeeze(beta_tt(2,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');

subplot(3,1,3);
plot(Date,squeeze(beta_tt(3,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');



fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_filter_init_MSV_betas','-dpdf');



figure('Name','Learning: Shock coefficients');
subplot(2,2,1);
plot(Date,squeeze(cc_tt(1,1,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');

subplot(2,2,2);
plot(Date,squeeze(cc_tt(1,2,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');

subplot(2,2,3);
plot(Date,squeeze(cc_tt(2,1,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');



subplot(2,2,4);
plot(Date,squeeze(cc_tt(2,2,2:end)),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_filter_init_MSV_shockCoef','-dpdf');


figure;
plot(largestEig);
hold on;
plot(pr_flag);