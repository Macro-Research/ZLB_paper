clear;clc;close all;
load('impREE_rise.mat');

figure('Name','impulse responses-consumption','units','normalized','outerposition',[0 0 1 1]);
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(impREE_normal(2:end,1,index),'*');
hold on;
plot(impREE_zlb(2:end,1,index),'lineWidth',3);
title(shocks(jj));
end
legend('normal','zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'rise_impresp_cons','-dpdf');

figure('Name','impulse responses-investment','units','normalized','outerposition',[0 0 1 1]);
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(impREE_normal(2:end,2,index),'*');
hold on;
plot(impREE_zlb(2:end,2,index),'lineWidth',3);
title(shocks(jj));
end
legend('normal','zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'rise_impresp_inve','-dpdf');


figure('Name','impulse responses-output','units','normalized','outerposition',[0 0 1 1]);
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(impREE_normal(2:end,3,index),'*');
hold on;
plot(impREE_zlb(2:end,3,index),'lineWidth',3);
title(shocks(jj));
end
legend('normal','zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'rise_impresp_output','-dpdf');

figure('Name','impulse responses-inflation','units','normalized','outerposition',[0 0 1 1]);
shocks={'\eta_a','\eta_b','\eta_g','\eta_i','\eta_r','\eta_p','\eta_w'};
index=0;
for jj=[1 2 3 4 6 7]
    index=index+1;
subplot(3,2,index);
plot(impREE_normal(2:end,4,index),'*');
hold on;
plot(impREE_zlb(2:end,4,index),'lineWidth',3);
title(shocks(jj));
end
legend('normal','zlb');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'rise_impresp_pinf','-dpdf');
