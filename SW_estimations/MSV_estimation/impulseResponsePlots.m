clear;clc;%close all;
load('kf_output_workspace.mat');
periods=100;
numPeriods=50;
numVar=24;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%mc zcap rk k1   q c inve y lab pinf w r kp 
   
beta_tt=0*eye(numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=0*ones(numEndo,numShocks);
rr_tt=5*eye(numBackward+numExo+1);
load('initial_beliefs_msv.mat');
beta_tt(forward_indices,backward_indices)=beta_init;
cc_tt(forward_indices,:)=cc_init;
%impulse order: period-dataset, period-impulse,numVars,numShocks);

for jj=2:T

gamma1_1_tilde=AA1_inv*(BB1+CC1*beta_tt^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_1_tilde=(AA1_inv*CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta_tt^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_2_tilde=(AA2_inv*CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2_inv*DD2;


gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1];
gamma2_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_1_tilde;zeros(numExo,1)];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA1_inv*EE1;FF1];
gamma4_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(numEndo,numExo);GG1];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2];
gamma2_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_2_tilde;zeros(numExo,1)];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA2_inv*EE2;FF2];
gamma4_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(numEndo,numExo);GG2];
   
imp11(jj,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(jj,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);

%    alpha_old=alpha_tt;beta_old=beta_tt;cc_old=cc_tt;    
%    thetaOld=[alpha_tt(forward_indices) beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];
%    [theta,rr_tt ,largestEig(tt),pr_flag(tt)] =msv_learning2(S_filtered(tt,forward_indices)',[1;S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,numBackward,backward_indices);
% % % 
% alpha_tt(forward_indices)=theta(1,:)';
%  beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
%  cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
% 
%             try
% largest_eig2(tt)=max(abs(eigs(gamma1_1,1)),abs(eigs(gamma1_2,1)));
%             catch
%                 largest_eig2(tt)=1.01;
%             end
% %             
%     if largest_eig2(tt)>0.99
%        alpha_tt=alpha_old;
%     beta_tt=beta_old;
%     cc_tt=cc_old;
%     pr_flag(tt)=1;
%     end
%  
%  learning_filtered(tt,:,:)=theta;

end

%%%%%%%%%%
load('impREE.mat');

figure('Name','impulse responses-consumption');
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp11(56,2:end,6,jj),'--');
hold on;
plot(imp22(127,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE(:,1,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');

figure('Name','impulse responses-investment');
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp11(56,2:end,7,jj),'--');
hold on;
plot(imp22(127,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE(:,2,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');


figure('Name','impulse responses-output');
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp11(56,2:end,8,jj),'--');
hold on;
plot(imp22(127,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE(:,3,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');

figure('Name','impulse responses-inflation');
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp11(56,2:end,10,jj),'--');
hold on;
plot(imp22(115,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE(:,4,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
