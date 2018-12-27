clear;clc;close all;
load us_dataset.mat;
lagged_r=robs(2:247);
pinf=pinfobs(1:246);

window=20;
coef=zeros(length(lagged_r)-window);
for i=1:length(pinf)-window;
    
   coef(i)=corr(lagged_r(1+i:window+i),pinf(1+i:window+i));
    
end

figure;
plot(coef,'lineWidth',5);