clear;clc;%close all;
load('kf_output_workspace_october15.mat');
periods=100;
numPeriods=T;
numVar=24;numShocks=7;
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
%mc zcap rk k1   q c inve y lab pinf w r kp 
   
beta_tt=0*eye(numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=0*ones(numEndo,numShocks);
  %rr_tt=5*ones(numBackward+numExo+1);
 load('initial_beliefs_msv.mat');
  beta_tt(forward_indices,backward_indices)=beta_init;
   cc_tt(forward_indices,:)=cc_init;
   rr_tt=rr_init;
  rr_tt=triu(rr_tt);
  rr_tt=(rr_tt+rr_tt')/2;
%impulse order: period-dataset, period-impulse,numVars,numShocks);

update_matrices;

for tt=2:T

%     
 update_beliefs;
 update_matrices;
% check_eigenvalue;
% revert_if_explosive;
%     
% update_matrices;
% check_eigenvalue; 

imp11(tt,:,:,:)=impulse_response(parameters(:,1),gamma1_1,gamma3_1,periods);
imp22(tt,:,:,:)=impulse_response(parameters(:,2),gamma1_2,gamma3_2,periods);
imp_averaged(tt,:,:,:)=pp_filtered(tt)*imp11(tt,:,:,:)+(1-pp_filtered(tt))*imp22(tt,:,:,:);

end

%%%%%%%%%%
%compare impulse responses under REE with learning, zlb & normal regimes
load('impREE.mat');
impREE_stdev=impREE_stdev(1:periods,:,:);
impREE_unitShock=impREE_unitShock(1:periods,:,:);

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,6,jj),'--');
hold on;
plot(imp_averaged(115,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,1,index),'*');
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
plot(imp_averaged(90,2:end,7,jj),'--');
hold on;
plot(imp_averaged(115,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,2,index),'*');
title(shocks(jj));
end
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
plot(imp_averaged(90,2:end,8,jj),'--');
hold on;
plot(imp_averaged(115,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,3,index),'*');
title(shocks(jj));
end
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
plot(imp_averaged(90,2:end,10,jj),'--');
hold on;
plot(imp_averaged(115,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE_unitShock(:,4,index),'*');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_pinf_reeComp','-dpdf');


%-----------------time-varying plots-----------------

refLine1=linspace(1,periods,periods);refLine1=repmat(refLine1,[numPeriods 1]);
%refLine2=linspace(1,numPeriods,numPeriods);refLine2=repmat(refLine2,[periods 1])';
refLine2=linspace(startDate,endDate,numPeriods);refLine2=repmat(refLine2,[periods 1])';

figure('Name','time varying impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,6,jj));
  ylim([startDate endDate])
   datetick('y','yyyy','keeplimits');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_cons_3d','-dpdf');

figure('Name','time varying impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,7,jj));
  ylim([startDate endDate])
   datetick('y','yyyy','keeplimits');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_inv_3d','-dpdf');


figure('Name','time varying impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,8,jj));
  ylim([startDate endDate])
   datetick('y','yyyy','keeplimits');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_output_3d','-dpdf');


figure('Name','time varying impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot3(refLine1(2:end,2:end)',refLine2(2:end,2:end)',imp_averaged(2:end,2:end,10,jj));
  ylim([startDate endDate])
   datetick('y','yyyy','keeplimits');
title(shocks(jj));
end
legend('normal','zlb','ree');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_pinf_3d','-dpdf');


%------------------------------compare with rise output
%compare impulse responses under REE-MS with learning, zlb & normal regimes
load('impREE_subSample_rise.mat');

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,6,jj),'--');
hold on;
plot(imp_averaged(115,2:end,6,jj),'lineWidth',3);
hold on;
plot(impREE_normal(2:end,1,index),'*');
hold on;
plot(impREE_zlb(2:end,1,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_cons_riseComp','-dpdf');


figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,7,jj),'--');
hold on;
plot(imp_averaged(115,2:end,7,jj),'lineWidth',3);
hold on;
plot(impREE_normal(2:end,2,index),'*');
hold on;
plot(impREE_zlb(2:end,2,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_inv_riseComp','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,8,jj),'--');
hold on;
plot(imp_averaged(115,2:end,8,jj),'lineWidth',3);
hold on;
plot(impREE_normal(2:end,3,index),'*');
hold on;
plot(impREE_zlb(2:end,3,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_output_riseComp','-dpdf');


figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(imp_averaged(90,2:end,10,jj),'--');
hold on;
plot(imp_averaged(115,2:end,10,jj),'lineWidth',3);
hold on;
plot(impREE_normal(2:end,4,index),'*');
hold on;
plot(impREE_zlb(2:end,4,index),'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_pinf_riseComp','-dpdf');

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
cumul=cumsum(imp_averaged(90,2:end,6,jj));
plot(cumul,'--');
hold on;
cumul=cumsum(imp_averaged(115,2:end,6,jj));
plot(cumul,'lineWidth',3);
hold on;
cumul=cumsum(impREE_normal(2:end,1,index));
plot(cumul,'*');
hold on;
cumul=cumsum(impREE_zlb(2:end,1,index));
plot(cumul,'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_cons_riseComp_cumul','-dpdf');


figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
cumul=cumsum(imp_averaged(90,2:end,7,jj));
plot(cumul,'--');
hold on;
cumul=cumsum(imp_averaged(115,2:end,7,jj));
plot(cumul,'lineWidth',3);
hold on;
cumul=cumsum(impREE_normal(2:end,2,index));
plot(cumul,'*');
hold on;
cumul=cumsum(impREE_zlb(2:end,2,index));
plot(cumul,'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_inv_riseComp_cumul','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);

index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
cumul=cumsum(imp_averaged(90,2:end,8,jj));
plot(cumul,'--');
hold on;
cumul=cumsum(imp_averaged(115,2:end,8,jj));
plot(cumul,'lineWidth',3);
hold on;
cumul=cumsum(impREE_normal(2:end,3,index));
plot(cumul,'*');
hold on;
cumul=cumsum(impREE_zlb(2:end,3,index));
plot(cumul,'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_output_riseComp_cumul','-dpdf');


figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
cumul=cumsum(imp_averaged(90,2:end,10,jj));
plot(cumul,'--');
hold on;
cumul=cumsum(imp_averaged(115,2:end,10,jj));
plot(cumul,'lineWidth',3);
hold on;
cumul=cumsum(impREE_normal(2:end,4,index));
plot(cumul,'*');
hold on;
cumul=cumsum(impREE_zlb(2:end,4,index));
plot(cumul,'.-');
title(shocks(jj));
end
legend('learning-normal','learning-zlb','ree-normal','ree-zlb');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MSV_impresp_pinf_riseComp_cumul','-dpdf');





