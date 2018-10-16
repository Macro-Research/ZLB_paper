clear;clc;close all;
tic
N=1e6;
draws=nan(N,1);
draws=gpuArray(draws);
parfor i=1:N
    disp(i)
    draws(i)=rand;
end
toc