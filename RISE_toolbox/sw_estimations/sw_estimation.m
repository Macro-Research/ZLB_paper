clear;clc;close all;

load('raf_dataset.mat');
vnames=fieldnames(dataset);
dataset=[dy dc dinve dw labobs robs pinfobs];

data_start='1966q1';
data=ts(data_start,dataset,...
    {'dy','dc','dinve','dw','labobs','robs','pinfobs'});

sw=rise('usmodel_tolga_switching','data',data,...
    'estim_start_date',obs2date(data_start,71),...
    'kf_presample',4);

sw=estimate(sw);