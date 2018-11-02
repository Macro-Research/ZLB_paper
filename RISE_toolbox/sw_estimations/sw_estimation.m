clear;clc;%close all;
% rng(11);
load('raf_dataset.mat');
vnames=fieldnames(dataset);
dataset=[dy dc dinve dw labobs robs pinfobs];

data_start='1947q1';
% data=ts(data_start,dataset,...
%     {'dy','dc','dinve','dw','labobs','robs','pinfobs', 'GS2'});

data=ts(data_start,dataset,...
    {'dy','dc','dinve','dw','labobs','robs','pinfobs'});

% sw=rise('usmodel_tolga_switching','data',data,...
%     'estim_start_date',obs2date(data_start,147),...
%     'kf_presample',4);

% sw=rise('usmodel_tolga_switching','data',data,...
%     'estim_start_date',obs2date(data_start,147),...
%     'kf_presample',4,'kf_init_variance',10,...
%     'steady_state_unique',false,'steady_state_imposed',true);

sw=rise('usmodel_GS2_f20_RS_ZLB_Rcte_1stst_nobreak_eme','data',data,...
    'estim_start_date',obs2date(data_start,147),...
    'kf_presample',4,'kf_init_variance',10,...
    'steady_state_imposed',true,'steady_state_unique',false);

% sw=rise('usmodel_switching_bshock','data',data,...
%     'estim_start_date',obs2date(data_start,147),...
%     'kf_presample',4,'kf_init_variance',10,...
%     'steady_state_unique',false,'steady_state_imposed',false);


% sw=rise('usmodel_switching_bshock','data',data,...
%     'estim_start_date',obs2date(data_start,147),...
%     'kf_presample',4,'kf_init_variance',10);


% sw=rise('usmodel_tolga_switching','data',data,...
%     'estim_start_date',obs2date(data_start,147),'kf_presample',4,'kf_init_variance',10,'steady_state_imposed',true);


%% solve the model
[sw,retcode]=solve(sw,'fix_point_maxiter',10000,'fix_point_TolFun',1e-10);
[sw,loglik,likt_rs_orig]=filter(sw);

%% print results
%sw.print_solution()
%figure;
%subplot(2,1,1)
%plot(sw.filtering.filtered_regime_probabilities.regime_2)
%subplot(2,1,2)
%plot(sw.filtering.smoothed_regime_probabilities.regime_2)

%% sw=estimate(sw,'optimizer',@fmincon);
%sw=estimate(sw,'optimizer',@fminunc);
  %sw=estimate(sw,'kf_filtering_level',0,'solve_check_stability',false);

figure('Name','filtered and updated regime');
subplot(2,1,1)
plot(sw.filtering.filtered_regime_probabilities.regime_2)
subplot(2,1,2)
plot(sw.filtering.smoothed_regime_probabilities.regime_2);

figure('Name','filtered interest rate and b-shock');
subplot(2,1,1);
plot(sw.filtering.filtered_variables.r(:,1));
hold on;
plot(sw.filtering.filtered_variables.r(:,2));
legend('normal regime','zlb regime');
subplot(2,1,2);
plot(sw.filtering.filtered_variables.b(:,1));
hold on;
plot(sw.filtering.filtered_variables.b(:,2));
legend('normal regime','zlb regime');


rise_forecasts=double([sw.filtering.Expected_filtered_variables.dy...
    sw.filtering.Expected_filtered_variables.dc...
    sw.filtering.Expected_filtered_variables.dinve...
    sw.filtering.Expected_filtered_variables.dw...
    sw.filtering.Expected_filtered_variables.pinfobs...
    sw.filtering.Expected_filtered_variables.robs...
    sw.filtering.Expected_filtered_variables.labobs]);


T=length(rise_forecasts);
startDate=datenum('01-01-1985');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,T);
obs=dataset(end-length(rise_forecasts)+1:end,:);
%inflation & hours worked mixed, switch them 
obs_aux=obs;
obs_aux(:,5)=obs(:,7);
obs_aux(:,7)=obs(:,5);
obs=obs_aux;

obs_names=[{'Output Growth','Consumption Growth','Real Investment Growth','Real Wage Growth',...
    'Inflation','Fed Funds Rate','Hours Worked'}];
yylim=[min(obs)',max(obs)'];
figure('Name','Rise Forecast Errors');
%
for jj=1:4
%for jj=1:length(obs(1,:));
subplot(4,1,jj);
plot(Date(2:end),rise_forecasts(1:end-1,jj),'-','color','black','lineWidth',1);
hold on;
plot(Date(2:end),obs(2:end,jj),'-','color','red','lineWidth',1);
title(obs_names(jj));
ylim([yylim(jj,:)]);
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
end
legend('forecast','observable');

steady_state=table(sw.endogenous.name',full(sw.solution.ss{1,1}),full(sw.solution.ss{1,2}));
disp(steady_state);