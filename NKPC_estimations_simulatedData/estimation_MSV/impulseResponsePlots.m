clear;clc;%close all;
load('kf_output_workspace.mat');
periods=40;
numPeriods=50;
numVar=5;numShocks=3;
shocks={'\eta_y','\eta_{pi}','\eta_r'};
vars={'y','pi','r','u_y','u_{\pi}'};
beta1=zeros(3,3);
index=0;
for jj=T-numPeriods:T-1
    index=index+1;
    beta1(:,3)=beta_tt(:,jj);
    cc1=cc_tt(:,:,jj);
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




index=0;
 figure('Name','Impulse Responses-MSV Model');
       legend('normal regime','zlb regime'); 
for ii=1:2
    for jj=1:2

        index=index+1;
         subplot(2,2,index);
%         subplot(numVar-2,numShocks,index);
%  plot3(refLine1',refLine2',imp11(:,:,ii,jj)','color','red','lineWidth',0.5);
%  hold on;
%  plot3(refLine1',refLine2',imp22(:,:,ii,jj)','color','black','lineWidth',0.5);
plot(imp11(5,:,ii,jj),'lineWidth',3,'color','blue');
hold on;
plot(imp22(5,:,ii,jj),'lineWidth',3,'color','green');
       title([vars(ii) shocks(jj)]);
    end
end


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_ree_init_MSV_IR','-dpdf');

% index=0;
% for ii=1:numVar-2
%     for jj=1:numShocks
%         figure;
%         index=index+1;
%         subplot(numVar-2,numShocks,index);
% plot3(refLine1',refLine2',imp33(:,:,ii,jj)','k','lineWidth',1);
% legend('weighted average');
%        title([vars(ii) shocks(jj)]);
%     end
% end

