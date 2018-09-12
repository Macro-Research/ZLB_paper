clear;clc;%close all;
load('kf_output_workspace.mat');
periods=40;
numPeriods=50;
numVar=17;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%vars={'y','pi','r','u_y','u_{\pi}'};
alpha1=0*ones(numVar,1);
beta1=0*eye(numVar);
rr=repmat(eye(2),[1 1 numVar]);
% beta1(3,3)=0.9;beta1(5,5)=0.9 ;beta1(6,6)=0.9;beta1(7,7)=0.9;beta1(9,9)=0.9;
% beta1(10,10)=0.9;beta1(11,11)=0.9;
load('AR1_initial_beliefs.mat');
for jj=1:length(beta_init)
    beta1(forward_indices(jj),forward_indices(jj))=beta_init(jj);
    rr(:,:,forward_indices(jj))=rr_init(:,:,jj);
end
for jj=T-numPeriods:T-1

gamma1_1_tilde=AA1_inv*(BB1+CC1*beta1^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta1)*alpha1;gamma3_1_tilde=(AA1_inv*CC1)*(beta1*cc1+cc1*rho1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta1^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta1)*alpha1;gamma3_2_tilde=(AA2_inv*CC2)*(beta1*cc1+cc1*rho2)+AA2_inv*DD2;

gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho1];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA1_inv*EE1;FF1];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),rho2];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA2_inv*EE2;FF2];



imp11(index,:,:,:)=impulse_response(parameters,gamma1_1,gamma3_1,periods);
imp22(index,:,:,:)=impulse_response(parameters,gamma1_2,gamma3_2,periods);
end

refLine1=linspace(1,periods,periods);refLine1=repmat(refLine1,[numPeriods 1]);
refLine2=linspace(1,numPeriods,numPeriods);refLine2=repmat(refLine2,[periods 1])';




