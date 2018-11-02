clear;clc;close all;
load REE_filtered_vars.mat;
T=size(dataset,1);numObs=7;
startDate=datenum('01-01-1966');
endDate = datenum('01-12-2016');
obs_indices=[14 15 16 17 10 12 9];
Date=linspace(startDate,endDate,T);
yylim=[min(dataset)',max(dataset)'];
obs_names=[{'Output Growth','Consumption Growth','Real Investment Growth','Real Wage Growth',...
    'Inflation','Fed Funds Rate','Hours Worked'}];

figure('Name','Forecast Errors','units','normalized','outerposition',[0 0 1 1]);
for jj=1:length(obs_indices)

subplot(7,1,jj);
plot(Date,REE_filtered_vars(:,jj),'-','color','black','lineWidth',1);
hold on;
plot(Date,dataset(:,jj),'-','color','red','lineWidth',1);
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
title(obs_names(jj));
grid on;
grid minor;
ylim([yylim(jj,:)]);
if jj==length(obs_indices)
    legend('forecast','observable');
end


end
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'REE_forecast_errors','-dpdf'); 