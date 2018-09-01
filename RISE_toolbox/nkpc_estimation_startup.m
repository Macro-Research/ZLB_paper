clear classes
clc
close all

mm=rise('nkpc_estimation_ms.rs');

%mm=solve(mm);

%mm.print_solution

load('full_dataset.mat'); 
dataset=[gap_cbo,pinfobs,robs];
dataset=dataset(45:end,:);
startdate='1966q2';
dataset=ts(startdate,dataset,{'gap_cbo','pinfobs','robs'});
enddate=obs2date(startdate,size(dataset,1));
vnames=dataset.varnames;
mm=set(mm,'data',dataset,'estim_end_date',enddate);


mm=estimate(mm);


plot_priors(mm)
% plot_posteriors(m,Results.pop)
% plot_priors_and_posteriors(m,Results.pop)



save('nkpc_estimation_results.mat');