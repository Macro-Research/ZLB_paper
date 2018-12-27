clear;clc;close all;
load newest_dataset.mat;
figure('Name','Point of Departure','units','normalized','outerposition',[0 0 1 1]);
T=108;
before_GFC=mean(dy(140:214));
after_GFC=mean(dy(218:end));
startDate=datenum('01-01-1990');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,T);
tickNum=10;
tickSize=(Date(end)-Date(1))/tickNum;

subplot(2,1,1);
plot(Date,dy(end-T+1:end),'color','black','lineWidth',5);
hold on;
% hold on;
% plot(Date,before_GFC*ones(T,1),'--','lineWidth',5,'color','red');
% hold on;
% plot(Date,after_GFC*ones(T,1),'--','lineWidth',5,'color','red');
xlim([startDate endDate]);
datetick('x','yy','keeplimits');
line([Date(1),Date(75)],[before_GFC,before_GFC],'lineStyle','--','lineWidth',5','Color','blue');
hold on;
line([Date(79),Date(end)],[after_GFC,after_GFC],'lineStyle','--','lineWidth',5','Color','red');
title('Output Growth','FontSize',35);
set(gca,'FontSize',20);

% subplot(3,1,2);
% plot(Date,pinfobs(end-T+1:end),'color','black','lineWidth',5);
% xlim([startDate endDate]);
% datetick('x','yyyy','keeplimits');
% title('Inflation','FontSize',35);
% set(gca,'FontSize',35);

subplot(2,1,2);
plot(Date,robs(end-T+1:end),'color','black','lineWidth',5);
xlim([startDate endDate]);
datetick('x','yy','keeplimits');
title('Interest Rates','FontSize',35);
set(gca,'FontSize',20);
hold on;
yl=ylim;
line([Date(end-31),Date(end-31)],[yl(1),yl(2)],'lineStyle','--','lineWidth',5','Color','red');
hold on;
line([Date(end-2),Date(end-2)],[yl(1),yl(2)],'lineStyle','--','lineWidth',5','Color','red');


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'point_of_departure','-dpdf'); 