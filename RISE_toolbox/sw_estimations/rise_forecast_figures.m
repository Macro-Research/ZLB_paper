clear;clc;close all;
load estimation_results_rs_fullsample.mat;
load yylim.mat


rise_forecasts=double([sw.filtering.Expected_filtered_variables.dy...
    sw.filtering.Expected_filtered_variables.dc...
    sw.filtering.Expected_filtered_variables.dinve...
    sw.filtering.Expected_filtered_variables.dw...
    sw.filtering.Expected_filtered_variables.pinfobs...
    sw.filtering.Expected_filtered_variables.robs...
    sw.filtering.Expected_filtered_variables.labobs]);
T=length(rise_forecasts);
startDate=datenum('01-01-1966');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,T);

obs=dataset(end-length(rise_forecasts)+1:end,:);
%inflation & hours worked mixed, switch them 
obs_aux=obs;
obs_aux(:,5)=obs(:,7);
obs_aux(:,7)=obs(:,5);
obs=obs_aux;

%clearvars -except rise_forecasts obs;
obs_names=[{'Output Growth','Consumption Growth','Real Investment Growth','Real Wage Growth',...
    'Inflation','Fed Funds Rate','Hours Worked'}];
yylim=[min(obs)',max(obs)'];
figure('Name','Rise Forecast Errors','units','normalized','outerposition',[0 0 1 1]);
% for jj=1:length(obs(1,:));
for jj=1:4;%only growth rates
subplot(4,1,jj);
plot(Date(2:end),rise_forecasts(1:end-1,jj),'-','color','black','lineWidth',1);
hold on;
plot(Date(2:end),obs(2:end,jj),'-','color','red','lineWidth',2);
title(obs_names(jj),'FontSize',35);
ylim([yylim(jj,:)]);
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
end
leg=legend('forecast','observable');
leg.FontSize=35;
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'rise_forecast_errors','-dpdf');

