clear;clc;close all;
load('raf_dataset.mat');
first_obs=180;


recursive_mean=zeros(7,length(dy)-first_obs);
for jj=1:length(dy)-first_obs;
    recursive_mean(:,jj)=mean([dy(first_obs:first_obs+jj)...
    dc(first_obs:first_obs+jj)...
    dinve(first_obs:first_obs+jj)...
    dw(first_obs:first_obs+jj)...
    labobs(first_obs:first_obs+jj)...
    pinfobs(first_obs:first_obs+jj)...
    robs(first_obs:first_obs+jj)]);
    
end


for jj=1:7
    figure;
    %subplot(7,1,jj);
    plot(recursive_mean(jj,:));
end

figure('Name','Inflation');
plot(pinfobs(first_obs:end));