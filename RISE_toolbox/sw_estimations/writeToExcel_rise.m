clear;clc;close all;
load('estimation_results_rs.mat');
xx=sw.estimation.posterior_maximization.mode;

output_file='Estimation_Results.csv';
 output_sheet='Estimation_Results';
 output=round(xx',2);
 laplace1_=oo_.MarginalDensity.LaplaceApproximation;
 
%  xlswrite(output_file,priorMean,output_sheet,'c15');
%  xlswrite(output_file,priorStd,output_sheet,'d15');
 xlswrite(output_file,output,output_sheet,'f3');
 xlswrite(output_file,laplace1_,output_sheet,'f45');