clear;clc;%close all;
load('kf_output_workspace.mat');
periods=200;
numPeriods=50;
numVar=24;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%mc zcap rk k1   q c inve y lab pinf w r kp 
   
alpha1=0*ones(numVar,1);
beta1=0*eye(numVar);
rr=repmat(eye(2),[1 1 numVar]);
load('AR1_initial_beliefs.mat');
%impulse order: period-dataset, period-impulse,numVars,numShocks);
for jj=1:length(beta_init)
    beta1(forward_indices(jj),forward_indices(jj))=beta_init(jj);
    rr(:,:,forward_indices(jj))=rr_init(:,:,jj);
end
for jj=2:T

gamma1_1=AA1_inv*(BB1+CC1*beta1^2);gamma2_1=AA1_inv*CC1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=AA1_inv*DD1;
gamma1_2=AA2_inv*(BB2+CC2*beta1^2);gamma2_2=AA2_inv*CC2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=AA1_inv*DD2;


imp11(jj,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(jj,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);

  for kk=[forward_indices]

beta1(kk,kk)=learning_filtered(kk,jj-1,2);
  end

end

figure;
for jj=1:7
subplot(4,2,jj);
plot(imp11(56,2:end,15,jj),'--');
hold on;
plot(imp22(127,2:end,15,jj),'lineWidth',3);
end
legend('normal','zlb');



