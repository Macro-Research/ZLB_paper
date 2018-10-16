clear;clc;%close all;
load('kf_output_workspace_30sept.mat');
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
rr_tt=rr_init;
%impulse order: period-dataset, period-impulse,numVars,numShocks);

update_matrices;

for jj=2:T

    
update_beliefs;
update_matrices;
check_eigenvalue;
revert_if_explosive;
    
update_matrices;
check_eigenvalue; 

imp11(jj,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(jj,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);
imp_averaged(jj,:,:,:)=pp_filtered(jj)*imp11(jj,:,:,:)+(1-pp_filtered(jj))*imp22(jj,:,:,:);

end

%%%%%%%%%%
load('impREE.mat');

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(80,2:end,6,jj),'--');
hold on;
plot(imp_averaged(190,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE(:,1,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_cons_reeComp','-dpdf');



figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(80,2:end,7,jj),'--');
hold on;
plot(imp_averaged(190,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE(:,2,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_inve_reeComp','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(80,2:end,8,jj),'--');
hold on;
plot(imp_averaged(190,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE(:,3,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_output_reeComp','-dpdf');




figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(80,2:end,10,jj),'--');
hold on;
plot(imp_averaged(190,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE(:,4,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_pinf_reeComp','-dpdf');









%% %%%%%%%%%%
% %compare impulse responses under REE-MS with learning, zlb & normal regimes
% load('impREE_rise.mat');
% 
% figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
% 
% index=0;
% for jj=[1 2 3 4 6 7]
%     index=index+1;
% subplot(3,2,index);
% plot(imp_averaged(99,2:end,6,jj),'--');
% hold on;
% plot(imp_averaged(102,2:end,6,jj),'lineWidth',3);
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
% plot(imp_averaged(99,2:end,7,jj),'--');
% hold on;
% plot(imp_averaged(102,2:end,7,jj),'lineWidth',3);
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
% plot(imp_averaged(99,2:end,8,jj),'--');
% hold on;
% plot(imp_averaged(102,2:end,8,jj),'lineWidth',3);
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
% plot(imp_averaged(99,2:end,10,jj),'--');
% hold on;
% plot(imp_averaged(102,2:end,10,jj),'lineWidth',3);
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


