clear;clc;close all;
load('sw_estimation_results.mat');
startDate=datenum('01-01-1964');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,208);

regime_prob=double(sw.filtering.smoothed_regime_probabilities.regime_2);

figure;
subplot(2,1,1);
area(Date,regime_prob(1:end));
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
subplot(2,1,2);
plot(Date,robs(end-207:end));
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig,'SW_sigmaPoint_regimeProb','-dpdf');