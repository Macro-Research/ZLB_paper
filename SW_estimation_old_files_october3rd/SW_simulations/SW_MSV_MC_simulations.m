clear;clc;close all;

numSimul=100;
pr_flag_mean=nan(numSimul,1);
theta_size=[7 15];
theta_collection=nan(numSimul,theta_size(1),theta_size(2));

for nn=1:numSimul
    disp(nn);
     if nn>1
    disp(['previous pr_flag:',num2str(pr_flag_mean(nn-1))]);
     end
   [pr_flag_mean(nn),theta_collection(nn,:,:)]=SW_MSV_simulation_function(); 
    
end