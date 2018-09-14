clear;clc;close all;

load('raf_dataset.mat');
vnames=fieldnames(dataset);
dataset=[dy dc dinve dw labobs robs pinfobs];

data_start='1947q1';
data=ts(data_start,dataset,...
    {'dy','dc','dinve','dw','labobs','robs','pinfobs'});

% sw=rise('usmodel_tolga_switching','data',data,...
%     'estim_start_date',obs2date(data_start,147),...
%     'kf_presample',4);

sw=rise('usmodel_tolga_switching','data',data,...
    'estim_start_date',obs2date(data_start,151));

sw=estimate(sw);

load('sw_estimation_results.mat');