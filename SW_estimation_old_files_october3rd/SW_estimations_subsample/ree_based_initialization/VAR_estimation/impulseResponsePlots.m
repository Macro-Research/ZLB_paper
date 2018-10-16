clear;clc;%close all;
load('kf_output_workspace_30sept.mat');
periods=100;
numPeriods=50;
numVar=24;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%mc zcap rk k1   q c inve y lab pinf w r kp 
   
alpha1=0*ones(numVar,1);
beta1=0*eye(numVar);
rr=repmat(eye(2),[1 1 numVar]);
load('AR1_initial_beliefs.mat');
%impulse order: (period-dataset, period-impulse,numVars,numShocks);
 beta1(3,3)=0.98;beta1(5,5)=0.2 ;beta1(6,6)=0.98;beta1(7,7)=0.98;beta1(9,9)=0.98;
 beta1(10,10)=0.6;beta1(11,11)=0.98;
for jj=1:length(beta_init)
    %beta1(forward_indices(jj),forward_indices(jj))=beta_init(jj);
    rr(:,:,forward_indices(jj))=rr_init(:,:,jj);
end
for jj=2:T

gamma1_1=AA1_inv*(BB1+CC1*beta1^2);gamma2_1=AA1_inv*CC1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=AA1_inv*DD1;
gamma1_2=AA2_inv*(BB2+CC2*beta1^2);gamma2_2=AA2_inv*CC2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=AA1_inv*DD2;


imp11(jj,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(jj,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);
imp_averaged(jj,:,:,:)=imp11(jj,:,:,:)*pp_filtered(jj)+imp22(jj,:,:,:)*(1-pp_filtered(jj));

     pr_flag(tt)=0;
        beta_old=beta1;
        alpha_old=alpha1;
     for jj=[forward_indices]
 
 [alpha1(jj) beta1(jj,jj) rr(:,:,jj)] =...
           msv_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj)]',...
        alpha1(jj),beta1(jj,jj),rr(:,:,jj),gain);
   
% 
        if abs(beta1(jj,jj))>.99999
         pr_flag(tt)=1;
         beta1(jj,jj)=beta_old(jj,jj);
       
% 
        end
        
    end

end


%%%%%%%%%%
%compare impulse responses under REE with learning, zlb & normal regimes
load('impREE.mat');

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(160,2:end,6,jj),'--');
hold on;
plot(imp_averaged(190,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,1,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_cons_reeComp','-dpdf');

figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(160,2:end,7,jj),'--');
hold on;
plot(imp_averaged(190,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,2,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_inve_reeComp','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(160,2:end,8,jj),'--');
hold on;
plot(imp_averaged(190,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,3,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_output_reeComp','-dpdf');


figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(160,2:end,10,jj),'--');
hold on;
plot(imp_averaged(190,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,4,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'AR1_impresp_pinf_reeComp','-dpdf');


%%%%%%%%%%
%compare impulse responses under REE-MS with learning, zlb & normal regimes
% load('impREE_rise.mat');
% 
% figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,2:end,6,jj),'--');
% hold on;
% plot(imp_averaged(130,2:end,6,jj),'lineWidth',3);
% hold on;
% plot(impREE_normal(2:end,1,index),'*');
% hold on;
% plot(impREE_zlb(2:end,1,index),'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_cons','-dpdf');
% 
% figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,2:end,7,jj),'--');
% hold on;
% plot(imp_averaged(130,2:end,7,jj),'lineWidth',3);
% hold on;
% plot(impREE_normal(2:end,2,index),'*');
% hold on;
% plot(impREE_zlb(2:end,2,index),'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_inve','-dpdf');
% 
% 
% figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,2:end,8,jj),'--');
% hold on;
% plot(imp_averaged(130,2:end,8,jj),'lineWidth',3);
% hold on;
% plot(impREE_normal(2:end,3,index),'*');
% hold on;
% plot(impREE_zlb(2:end,3,index),'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_output','-dpdf');
% 
% 
% figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(90,2:end,10,jj),'--');
% hold on;
% plot(imp_averaged(130,2:end,10,jj),'lineWidth',3);
% hold on;
% plot(impREE_normal(2:end,4,index),'*');
% hold on;
% plot(impREE_zlb(2:end,4,index),'.-');
% title(shocks(jj));
% end
% legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'AR1_impresp_pinf','-dpdf');
% 
