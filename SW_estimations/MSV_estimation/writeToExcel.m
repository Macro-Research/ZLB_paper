

  output_file='Estimation_Results.csv';
 output_sheet='Estimation_Results';
%  xlswrite(output_file,priorMean,output_sheet,'c15');
%  xlswrite(output_file,priorStd,output_sheet,'d15');
 xlswrite(output_file,output,output_sheet,'e3');
 xlswrite(output_file,laplace_,output_sheet,'e22');