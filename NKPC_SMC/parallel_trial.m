clear;clc;close all;
tic
N=10e7;
draw=nan(N,1);

for i=1:N
    draw(i)=rand;
end
toc