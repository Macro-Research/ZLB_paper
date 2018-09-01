function[likl]=KF_MS(parameters)
   gain=parameters(18);
%  q_11=0.99;q_22=0.85; 

load('us_dataset.mat');
first_obs=40;burn_in=6;
dataset=[gap_hp,pinfobs,robs];
dataset=dataset(first_obs:end,:);
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
alpha1=0*ones(numVar-numExo,1);beta1=0.75*eye(numVar-numExo);cc1=zeros(numEndo,numExo);
beta1=[0,0,-0.636;0,0,-0.016;0,0,0.727];
cc1=[4.112,0.167,0.442;-1.399,6.261,1.667]';
% cc1(1:3,1:2)=ones(3,2);
% acf1=autocorr(dataset(:,1),1);acf2=autocorr(dataset(:,2),1);
% beta1(1,1)=acf1(2);
% beta1(2,2)=acf2(2);
% dataset_lag=lagmatrix(dataset,1);
% dataset_lag=dataset_lag(2:end,:);
%beta1(1:3,1:3)=0.5*corr(dataset(2:end,:),dataset_lag);
% beta1(1:3,1:3)=(dataset_lag'*dataset_lag)^(-1)*dataset_lag'*dataset(2:end,:);
 rr=eye(6);
% beta1(1,1)=0.5;beta1(2,2)=0.5;
%rr=zeros(2,2,2);rr(:,:,1)=0.2*eye(2);rr(:,:,2)=0.2*eye(2);%when all variables are included in plm
%rr=eye(numVar-numExo);%when exogenous variables are unobserved

%rr(:,:,1)=eye(2);rr(:,:,2)=eye(2);
% [A1 B1 C1 D1, E1 F1 G1]=NKPC_sysmat_regime1(param(:,1));
% [A2 B2 C2 D2 E2 F2 G2]=NKPC_sysmat_regime1(param(:,2));
[AA1 BB1 CC1 DD1 EE1 FF1 rho1 E1 F1 G1]=NKPC_sysmat_MSV(param(:,1));
[AA2 BB2 CC2 DD2 EE2 FF2 rho2 E2 F2 G2]=NKPC_sysmat_MSV(param(:,2));
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);


Sigma1=diag([param(end-2,1)^2;param(end-1,1)^2;param(end,1)^2]);
Sigma2=diag([param(end-2,2)^2;param(end-1,2)^2;param(end,2)^2]);
gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;


H1=0*diag(var(dataset(1:150,:)));H2=0.1*diag(var(dataset(150:end,:)));



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
S_collapse1=(pp_upd11*S_upd11+pp_upd21*S_upd21)/pp_collapse1;
S_collapse2=(pp_upd12*S_upd12+pp_upd22*S_upd22)/pp_collapse2;
%
P_collapse1=(pp_upd11*(P_upd11+(S_collapse1-S_upd11)*(S_collapse1-S_upd11)')+...
    pp_upd21*(P_upd21+(S_collapse1-S_upd21)*(S_collapse1-S_upd21)'))/pp_collapse1;

P_collapse2=(pp_upd12*(P_upd12+(S_collapse2-S_upd12)*(S_collapse2-S_upd12)')+...
    pp_upd22*(P_upd22+(S_collapse2-S_upd22)*(S_collapse2-S_upd22)'))/pp_collapse2;

%
S_filtered(tt,:)=pp_collapse1*S_collapse1+pp_collapse2*S_collapse2;
pp_filtered(tt)=pp_collapse1;

thetaOld=[alpha1 beta1 cc1];
[theta rr] =msv_learning2(S_filtered(tt,1:numEndo)',[1,S_filtered(tt-1,1:numEndo),S_filtered(tt,numEndo+1:end)]',thetaOld,rr,gain);
alpha1=theta(1,:)';beta1=theta(2:4,:)';cc1=theta(5:6,:)';

% [alpha1(1:3) beta1(1:3,1:3) rr(1:4,1:4)] =msv_learning(S_filtered(tt,1:3)',[1,S_filtered(tt-1,1:3)]'...
% ,alpha1(1:3),beta1(1:3,1:3),rr(1:4,1:4),gain);

% 
%  for jj=1:2
% %     
%   [alpha1(jj) beta1(jj,jj) rr(:,:,jj)] =msv_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj)]'...
%   ,alpha1(jj),beta1(jj,jj),rr(:,:,jj),gain);
%          learning_filtered(jj,tt-1,:)=[alpha1(jj),beta1(jj,jj)];
%  end




     for jj=1:2
    % [alpha1(jj) beta1(jj,jj) rr(:,:,jj)] =...
     %     msv_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj)]',...
      %alpha1(jj),beta1(jj,jj),rr(:,:,jj),gain);
% % % %        [alpha1(jj),beta1(jj,jj)]=l_SAC_CGL_nonRecursive(S_filtered(1:tt,jj),alpha1(jj),gain);
   %[alpha1(jj) beta1(jj,jj), rr(jj,jj)]=l_SAC_CGL...
   %(S_filtered(tt,jj),S_filtered(tt-1,jj),alpha1(jj),beta1(jj,jj),rr(jj,jj),gain);
% %  [alpha1(jj) beta1(jj,jj)]=l_SAC_DGL(S_filtered(1:tt,jj),tt);
      learning_filtered_alpha1(tt-1,:)=alpha1;
      learning_filtered_beta1(tt-1,:,:)=beta1;
      learning_filtered_cc1(tt-1,:,:)=cc1;
      end

%     for jj=1:2
%        [alpha1(jj) beta1(jj,jj) rr(jj)] =...
%            msv_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj)]',...
%        alpha1(jj),beta1(jj,jj),rr(jj),gain);
% % % %        [alpha1(jj),beta1(jj,jj)]=l_SAC_CGL_nonRecursive(S_filtered(1:tt,jj),alpha1(jj),gain);
% %   [alpha1(jj) beta1(jj,jj), rr(jj,jj)]=l_SAC_CGL...
% %       (S_filtered(tt,jj),S_filtered(tt-1,jj),alpha1(jj),beta1(jj,jj),rr(jj,jj),gain);
% %  [alpha1(jj) beta1(jj,jj)]=l_SAC_DGL(S_filtered(1:tt,jj),tt);
% % %     learning_filtered(jj,tt-1,:)=[alpha1(jj),beta1(jj,jj)];
%      end
% % if tt>10
%  [alpha1,beta1]=l_YW_CGL(S_filtered(1:tt,:)',alpha1,gain,numVar);
% end
%   [alpha1 beta1 rr] =msv_learning(S_filtered(tt,:)',[1,S_filtered(tt-1,:)]'...
%       ,alpha1,beta1,rr,gain);
%   [alpha1(1:3) beta1(1:3,1:3) rr] =...
%       msv_learning(S_filtered(tt,1:3)',[1,S_filtered(tt-1,1:3)]'...
%      ,alpha1(1:3),beta1(1:3,1:3),rr,gain);
% beta1=diag(diag(beta1));

 end

likl=-sum(log(likl(burn_in+1:end)));

 end

% figure;
% subplot(3,1,1);
% plot(learning_filtered_alpha1(:,1));
% subplot(3,1,2);
% plot(learning_filtered_alpha1(:,2));
% subplot(3,1,3);
% plot(learning_filtered_alpha1(:,3));
% 
% figure;
% subplot(2,2,1);
% plot(learning_filtered_beta1(:,1,1));
% subplot(2,2,2);
% plot(learning_filtered_beta1(:,1,2));
% subplot(2,2,3);
% plot(learning_filtered_beta1(:,2,1));
% subplot(2,2,4);
% plot(learning_filtered_beta1(:,2,2));
% 
% figure;
% subplot(2,2,1);
% plot(learning_filtered_cc1(:,1,1));
% subplot(2,2,2);
% plot(learning_filtered_cc1(:,1,2));
% subplot(2,2,3);
% plot(learning_filtered_cc1(:,2,1));
% subplot(2,2,4);
% plot(learning_filtered_cc1(:,2,2));
% 
% figure;
% plot(1-pp_filtered,'lineWidth',3);
% hold on;
% plot(dataset(:,3),'--');
%  grid on;
%  grid minor;