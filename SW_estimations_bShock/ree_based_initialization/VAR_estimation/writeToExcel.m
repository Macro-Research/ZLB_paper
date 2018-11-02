clear;clc;close all;
load('estimation_results_30sept.mat');

  output_file='Estimation_Results.csv';
 output_sheet='Estimation_Results';
 output=round(x,2);
%  xlswrite(output_file,priorMean,output_sheet,'c15');
%  xlswrite(output_file,priorStd,output_sheet,'d15');
 xlswrite(output_file,output,output_sheet,'e3');
 xlswrite(output_file,-laplace1_,output_sheet,'e45');