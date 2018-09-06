% PROGRAM FOR GENERATING CUMULATIVE INFLATION FORECAST ERRORS PLOT.

% WITH LEVEL VARIABLES, STANDARD DYNARE RE CODE DOESN'T PREDICT OBSERVABLE 
% VARIABLES SUCH AS pinfobs BUT ITS MODEL COUNTERPARTS SUCH AS pinf, WHICH 
% MEANS THAT I NEED TO ADD THE ESTIMATED CONSTANT TO THE VARIABLE.

% FILES T_AR2_no4_oldfile_newdata_176_cont_MCMC_results AND
% T_RE_oldfile_newdata_176obs_MCMC_results ARE RESULT FILES PRODUCED BY OUR
% ADAPTIVE LEARNING PACKAGE. IN ORDER TO MINIMIZE THE SIZE OF THE CODE 
% RELEASE, ALL RELEVANT VARIABLES ARE STORED INSTEAD IN A SINGLE DATA FILE.
% EXCEL FILES WITH SPF AND REAL TIME DATA ARE AVAILABLE FROM THE
% PHILADELPHIA FED SITE:
% http://www.philadelphiafed.org/research-and-data/real-time-center/real-ti
% me-data/data-files/files/PQvQd.xls 
%                              AND
% http://www.philadelphiafed.org/research-and-data/real-time-center/survey-
% of-professional-forecasters/historical-data/medianLevel.xls

t_in = 1947.5 + 70 / 4;
t_b = (t_in:0.25:2008.75)';

%%%%% START OF THE COMMENTED BLOCK

% % Preparing the variables: ADAPTIVE LEARNING
% 
% load T_AR2_no4_oldfile_newdata_176_cont_MCMC_results
% 
% % Constructing variables: INFLATION
% pinfobs_AL = oo_.FilteredVariables.pinfobs;
% pinfobs_AL_fcast = oo_.ForecastVariables.pinfobs;
% 
% % Preparing the variables: AVERAGE MODEL PREDICTION
% % Given that only one model is present in the Baseline, average model
% % prediction is just AR(2) prediction
% pinf_AVM_fcast = oo_.exps(4,:)' + oo_.posterior_mode.parameters.constepinf;
% 
% % Preparing the variables: RATIONAL EXPECTATIONS
% % DiffuseKalmanSmoother1 has been rewritten to resemble the code for
% % adaptive learning; however, standard DYNARE manner still remains whereby
% % all observables are saved without corresponding constants, which explains
% % the difference in treatment
% 
% load T_RE_oldfile_newdata_176obs_MCMC_results
% 
% % Constructing variables: INFLATION
% pinfobs_RE = oo_.FilteredVariables.pinfobs + oo_.posterior_mode.parameters.constepinf;
% pinfobs_RE_fcast = oo_.ForecastVariables.pinfobs;
% [m,n] = size(pinfobs_RE_fcast);
% pinfobs_RE_fcast = pinfobs_RE_fcast + oo_.posterior_mode.parameters.constepinf * ones(m,n);
% 
% % Preparing the variables: SPF
% % Variables are taken from the Median SPF file
% % Fifth column contains VAR3 which is 1Q ahead forecast. Fourth column is
% % VAR2 which is current quarter's nowcast.
% pgdp_SPF = xlsread('medianLevel','PGDP');
% pinf_SPF_fcast = 100 * log(pgdp_SPF(:,5) ./ pgdp_SPF(:,4));
% pinf_SPF_ncast = 100 * log(pgdp_SPF(:,4) ./ pgdp_SPF(:,3));
% % The first time real-time data for 1969:Q1 to 1968:Q4 inflation is available is
% % 1969:Q2. Construction: take two latest points from every vintage
% % % The construction below results in the pattern that's very close to the
% % one using final data; apparently, first revision contains all the
% % information already
% pgdp_RT = xlsread('PQvQd','P','p90:ft250');
% pgdp_RTl = xlsread('PQvQd','P','p89:ft249');
% pgdp_RT = diag(pgdp_RT);
% pgdp_RTl = diag(pgdp_RTl);
% % In 1996:Q1, one diagonal point is not available
% pinf_RT = 100 * log(pgdp_RT(1:107) ./ pgdp_RTl(1:107));
% pinf_RT = [pinf_RT; pinf_RT(end)];
% pinf_RT = [pinf_RT; 100 * log(pgdp_RT(109:end) ./ pgdp_RTl(109:end))];

%%%%% END OF THE COMMENTED BLOCK

% Loading the saved variables instead of constructing them
load Fig_3_data

% Every variable has its starting date, I need to put all of them into the
% same format. Pad by observables when SPF forecast is not available
% INFLATION: starts in 1968:Q4 (first forecast for 1969:Q1) vs 1965:Q1 for my sample, need to add 16 Q
% For nowcast, I need to start from 1969:Q2 which is nowcast for 1969:Q1
pinf_SPF_fcast = [pinfobs_AL(1:16); pinf_SPF_fcast];
pinf_SPF_fcast(length(pinfobs_AL)+1:end) = [];
pinf_SPF_ncast = [pinfobs_AL(1:16); pinf_SPF_ncast(2:end)];
pinf_SPF_ncast(length(pinfobs_AL)+1:end) = [];
pinf_RT = [pinfobs_AL(1:16); pinf_RT];
pinf_RT(length(pinfobs_AL)+1:end) = [];

tmp_AVM = pinf_AVM_fcast(4+1:end,1) - pinfobs_AL(4+1:end);
tmp_RE  = pinfobs_RE_fcast(4+1:end-4,1) - pinfobs_AL(4+1:end);    
tmp_RTf   = pinf_SPF_fcast(4+1:end) - pinf_RT(4+1:end);
plots = [tmp_RE tmp_AVM tmp_RTf];
leg = {'RE';'AR(2)';'SPF\_fcast'};

figure(33)
plot(datenum(t_b(17:end),1,1),cumsum(plots(13:end,:)));
datetick('x',17,'keeplimits');
legend(leg)
title('Cumulative forecast errors, relative to 1^{st} release, inflation')
