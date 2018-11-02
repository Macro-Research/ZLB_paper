clear;clc;%close all;
load('kf_output_workspace_october22.mat');
periods=50;
numPeriods=T;
numVar=24;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%mc zcap rk k1   q c inve y lab pinf w r kp 
   
alpha1=0*ones(numVar,1);
beta1=0*eye(numVar);
rr=repmat(eye(2),[1 1 numVar]);
load('AR1_initial_beliefs.mat');
%impulse order: (period-dataset, period-impulse,numVars,numShocks);
%  beta1(3,3)=0.98;beta1(5,5)=0.2 ;beta1(6,6)=0.98;beta1(7,7)=0.98;beta1(9,9)=0.98;
%  beta1(10,10)=0.6;beta1(11,11)=0.98;
for jj=1:length(beta_init)
    beta1(forward_indices(jj),forward_indices(jj))=beta_init(jj);
    rr(:,:,forward_indices(jj))=rr_init(:,:,jj);
end

gamma1_1=AA1_inv*(BB1+CC1*beta1^2);
gamma2_1=AA1_inv*CC1*(eye(numVar)+beta1)*alpha1;
gamma3_1=AA1_inv*DD1;%gamma4_1=AA1_inv*EE1;

gamma1_2=AA2_inv*(BB2+CC2*beta1^2);
gamma2_2=AA2_inv*CC2*(eye(numVar)+beta1)*alpha1;
gamma3_2=AA2_inv*DD2;%gamma4_2=AA2_inv*EE2;

for tt=2:T



imp11(tt,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(tt,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);
imp_averaged(tt,:,:,:)=imp11(tt,:,:,:)*pp_filtered(tt)+imp22(tt,:,:,:)*(1-pp_filtered(tt));

%      pr_flag(tt)=0;
        beta_old=beta1;
        alpha_old=alpha1;
        rr_old=rr;
     for jj=[forward_indices]
 
 [alpha1(jj) beta1(jj,jj) rr(:,:,jj)] =...
           msv_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj)]',...
        alpha1(jj),beta1(jj,jj),rr(:,:,jj),gain);
     end
   
gamma1_1=AA1_inv*(BB1+CC1*beta1^2);
gamma1_2=AA2_inv*(BB2+CC2*beta1^2);


largest_eig1(tt)=abs(eigs(gamma1_1,1));
largest_eig2(tt)=abs(eigs(gamma1_2,1));
average_eig(tt)=ergodic_states(1)*largest_eig1(tt)+ergodic_states(2)*largest_eig2(tt);

if largest_eig1(tt)>1
    pr_flag(tt)=1;
    alpha1=alpha_old;
    beta1=beta_old;
    rr=rr_old;
% elseif largest_eig2(tt)>1
%         pr_flag(tt)=1;
%     alpha1=alpha_old;
%     beta1=beta_old;
    
elseif average_eig(tt)>1
    pr_flag(tt)=1;
    alpha1=alpha_old;
    beta1=beta_old;
    rr=rr_old;
end

gamma1_1=AA1_inv*(BB1+CC1*beta1^2);
gamma2_1=AA1_inv*CC1*(eye(numVar)+beta1)*alpha1;
gamma3_1=AA1_inv*DD1;%gamma4_1=AA1_inv*EE1;

gamma1_2=AA2_inv*(BB2+CC2*beta1^2);
gamma2_2=AA2_inv*CC2*(eye(numVar)+beta1)*alpha1;
gamma3_2=AA2_inv*DD2;%gamma4_2=AA2_inv*EE2;

end


%%%%%%%%%%
%compare impulse responses under REE with learning, zlb & normal regimes
% load('impREE.mat');
% impREE_stdev=impREE_stdev(1:periods,:,:);
% impREE_unitShock=impREE_unitShock(1:periods,:,:);
% 
% figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,1:periods,6,jj),'--');
% hold on;
% plot(imp_averaged(115,1:periods,6,jj),'lineWidth',3);
% hold on;
% plot(impREE_unitShock(1:periods,1,index),'*');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_cons_reeComp','-dpdf');
% 
% figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,1:periods,7,jj),'--');
% hold on;
% plot(imp_averaged(115,1:periods,7,jj),'lineWidth',3);
% hold on;
% plot(impREE_unitShock(1:periods,2,index),'*');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_inve_reeComp','-dpdf');
% 
% 
% figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,1:periods,8,jj),'--');
% hold on;
% plot(imp_averaged(115,1:periods,8,jj),'lineWidth',3);
% hold on;
% plot(impREE_unitShock(1:periods,3,index),'*');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_output_reeComp','-dpdf');
% 
% 
% figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,1:periods,10,jj),'--');
% hold on;
% plot(imp_averaged(115,1:periods,10,jj),'lineWidth',3);
% hold on;
% plot(impREE_unitShock(1:periods,4,index),'*');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_pinf_reeComp','-dpdf');
% 
% 
% %-----------------time-varying plots-----------------
% 
% refLine1=linspace(1,periods,periods);refLine1=repmat(refLine1,[numPeriods 1]);
% %refLine2=linspace(1,numPeriods,numPeriods);refLine2=repmat(refLine2,[periods 1])';
% refLine2=linspace(startDate,endDate,numPeriods);refLine2=repmat(refLine2,[periods 1])';
% 
% figure('Name','time varying impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,6,jj));
%   ylim([startDate endDate])
%    datetick('y','yyyy','keeplimits');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_cons_3d','-dpdf');
% 
% figure('Name','time varying impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,7,jj));
%   ylim([startDate endDate])
%    datetick('y','yyyy','keeplimits');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_inv_3d','-dpdf');
% 
% 
% figure('Name','time varying impulse responses-output','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,8,jj));
%   ylim([startDate endDate])
%    datetick('y','yyyy','keeplimits');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_output_3d','-dpdf');
% 
% 
% figure('Name','time varying impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,10,jj));
%   ylim([startDate endDate])
%    datetick('y','yyyy','keeplimits');
% title(shocks(jj));
% end
% legend('normal','zlb','ree');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_pinf_3d','-dpdf');


%------------------------------compare with rise output
%compare impulse responses under REE-MS with learning, zlb & normal regimes
load('impREE_subSample_rise.mat');
impREE_normal=impREE_normal(2:end,:,:);
impREE_zlb=impREE_zlb(2:end,:,:);

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,6,jj),'--');
hold on;
plot(imp_averaged(115,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE_normal(1:periods,1,index),'*');
hold on;
plot(impREE_zlb(1:periods,1,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_cons_riseComp','-dpdf');


figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,7,jj),'--');
hold on;
plot(imp_averaged(115,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE_normal(1:periods,2,index),'*');
hold on;
plot(impREE_zlb(1:periods,2,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_inv_riseComp','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,8,jj),'--');
hold on;
plot(imp_averaged(115,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE_normal(1:periods,3,index),'*');
hold on;
plot(impREE_zlb(1:periods,3,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_output_riseComp','-dpdf');


figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,10,jj),'--');
hold on;
plot(imp_averaged(115,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE_normal(1:periods,4,index),'*');
hold on;
plot(impREE_zlb(1:periods,4,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_pinf_riseComp','-dpdf');


%---------------------------------------------------------
%compare with RISE---cumulative responses


% figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% cumul=cumsum(imp_averaged(90,2:end,6,jj));
% plot(cumul,'--');
% hold on;
% cumul=cumsum(imp_averaged(115,2:end,6,jj));
% plot(cumul,'lineWidth',3);
% hold on;
% cumul=cumsum(impREE_normal(2:end,1,index));
% plot(cumul,'*');
% hold on;
% cumul=cumsum(impREE_zlb(2:end,1,index));
% plot(cumul,'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_cons_riseComp_cumul','-dpdf');
% 
% 
% figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% cumul=cumsum(imp_averaged(90,2:end,7,jj));
% plot(cumul,'--');
% hold on;
% cumul=cumsum(imp_averaged(115,2:end,7,jj));
% plot(cumul,'lineWidth',3);
% hold on;
% cumul=cumsum(impREE_normal(2:end,2,index));
% plot(cumul,'*');
% hold on;
% cumul=cumsum(impREE_zlb(2:end,2,index));
% plot(cumul,'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_inv_riseComp_cumul','-dpdf');
% 
% 
% figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% cumul=cumsum(imp_averaged(90,2:end,8,jj));
% plot(cumul,'--');
% hold on;
% cumul=cumsum(imp_averaged(115,2:end,8,jj));
% plot(cumul,'lineWidth',3);
% hold on;
% cumul=cumsum(impREE_normal(2:end,3,index));
% plot(cumul,'*');
% hold on;
% cumul=cumsum(impREE_zlb(2:end,3,index));
% plot(cumul,'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_output_riseComp_cumul','-dpdf');
% 
% 
% figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% cumul=cumsum(imp_averaged(90,2:end,10,jj));
% plot(cumul,'--');
% hold on;
% cumul=cumsum(imp_averaged(115,2:end,10,jj));
% plot(cumul,'lineWidth',3);
% hold on;
% cumul=cumsum(impREE_normal(2:end,4,index));
% plot(cumul,'*');
% hold on;
% cumul=cumsum(impREE_zlb(2:end,4,index));
% plot(cumul,'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_pinf_riseComp_cumul','-dpdf');
% 
% 
% 
% 
