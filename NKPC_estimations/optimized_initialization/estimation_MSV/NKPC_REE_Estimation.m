%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'NKPC_REE_Estimation';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('NKPC_REE_Estimation.log');
M_.exo_names = 'eps_y';
M_.exo_names_tex = 'eps\_y';
M_.exo_names_long = 'eps_y';
M_.exo_names = char(M_.exo_names, 'eps_pi');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_pi');
M_.exo_names_long = char(M_.exo_names_long, 'eps_pi');
M_.exo_names = char(M_.exo_names, 'eps_r');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_r');
M_.exo_names_long = char(M_.exo_names_long, 'eps_r');
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'u_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'u\_y');
M_.endo_names_long = char(M_.endo_names_long, 'u_y');
M_.endo_names = char(M_.endo_names, 'u_pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'u\_pi');
M_.endo_names_long = char(M_.endo_names_long, 'u_pi');
M_.endo_names = char(M_.endo_names, 'gap_cbo');
M_.endo_names_tex = char(M_.endo_names_tex, 'gap\_cbo');
M_.endo_names_long = char(M_.endo_names_long, 'gap_cbo');
M_.endo_names = char(M_.endo_names, 'pinfobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'pinfobs');
M_.endo_names_long = char(M_.endo_names_long, 'pinfobs');
M_.endo_names = char(M_.endo_names, 'robs');
M_.endo_names_tex = char(M_.endo_names_tex, 'robs');
M_.endo_names_long = char(M_.endo_names_long, 'robs');
M_.endo_partitions = struct();
M_.param_names = 'y_bar';
M_.param_names_tex = 'y\_bar';
M_.param_names_long = 'y_bar';
M_.param_names = char(M_.param_names, 'pi_bar');
M_.param_names_tex = char(M_.param_names_tex, 'pi\_bar');
M_.param_names_long = char(M_.param_names_long, 'pi_bar');
M_.param_names = char(M_.param_names, 'r_bar');
M_.param_names_tex = char(M_.param_names_tex, 'r\_bar');
M_.param_names_long = char(M_.param_names_long, 'r_bar');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'tau');
M_.param_names_tex = char(M_.param_names_tex, 'tau');
M_.param_names_long = char(M_.param_names_long, 'tau');
M_.param_names = char(M_.param_names, 'phi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pi');
M_.param_names_long = char(M_.param_names_long, 'phi_pi');
M_.param_names = char(M_.param_names, 'phi_y');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_y');
M_.param_names_long = char(M_.param_names_long, 'phi_y');
M_.param_names = char(M_.param_names, 'rho_y');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_y');
M_.param_names_long = char(M_.param_names_long, 'rho_y');
M_.param_names = char(M_.param_names, 'rho_pi');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_pi');
M_.param_names_long = char(M_.param_names_long, 'rho_pi');
M_.param_names = char(M_.param_names, 'rho_r');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_r');
M_.param_names_long = char(M_.param_names_long, 'rho_r');
M_.param_names = char(M_.param_names, 'a_1');
M_.param_names_tex = char(M_.param_names_tex, 'a\_1');
M_.param_names_long = char(M_.param_names_long, 'a_1');
M_.param_names = char(M_.param_names, 'a_2');
M_.param_names_tex = char(M_.param_names_tex, 'a\_2');
M_.param_names_long = char(M_.param_names_long, 'a_2');
M_.param_names = char(M_.param_names, 'a_3');
M_.param_names_tex = char(M_.param_names_tex, 'a\_3');
M_.param_names_long = char(M_.param_names_long, 'a_3');
M_.param_names = char(M_.param_names, 'b_11');
M_.param_names_tex = char(M_.param_names_tex, 'b\_11');
M_.param_names_long = char(M_.param_names_long, 'b_11');
M_.param_names = char(M_.param_names, 'b_12');
M_.param_names_tex = char(M_.param_names_tex, 'b\_12');
M_.param_names_long = char(M_.param_names_long, 'b_12');
M_.param_names = char(M_.param_names, 'b_13');
M_.param_names_tex = char(M_.param_names_tex, 'b\_13');
M_.param_names_long = char(M_.param_names_long, 'b_13');
M_.param_names = char(M_.param_names, 'b_21');
M_.param_names_tex = char(M_.param_names_tex, 'b\_21');
M_.param_names_long = char(M_.param_names_long, 'b_21');
M_.param_names = char(M_.param_names, 'b_22');
M_.param_names_tex = char(M_.param_names_tex, 'b\_22');
M_.param_names_long = char(M_.param_names_long, 'b_22');
M_.param_names = char(M_.param_names, 'b_23');
M_.param_names_tex = char(M_.param_names_tex, 'b\_23');
M_.param_names_long = char(M_.param_names_long, 'b_23');
M_.param_names = char(M_.param_names, 'b_31');
M_.param_names_tex = char(M_.param_names_tex, 'b\_31');
M_.param_names_long = char(M_.param_names_long, 'b_31');
M_.param_names = char(M_.param_names, 'b_32');
M_.param_names_tex = char(M_.param_names_tex, 'b\_32');
M_.param_names_long = char(M_.param_names_long, 'b_32');
M_.param_names = char(M_.param_names, 'b_33');
M_.param_names_tex = char(M_.param_names_tex, 'b\_33');
M_.param_names_long = char(M_.param_names_long, 'b_33');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 8;
M_.param_nbr = 22;
M_.orig_endo_nbr = 8;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'gap_cbo'};
options_.varobs(2)  = {'pinfobs'};
options_.varobs(3)  = {'robs'};
options_.varobs_id = [ 6 7 8  ];
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('NKPC_REE_Estimation_static');
erase_compiled_function('NKPC_REE_Estimation_dynamic');
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 4;
 0 5;
 1 6;
 2 7;
 3 8;
 0 9;
 0 10;
 0 11;]';
M_.nstatic = 5;
M_.nfwrd   = 0;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 3;
M_.ndynamic   = 3;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(22, 1);
M_.NNZDerivatives = [27; 0; -1];
M_.params( 1 ) = 0;
y_bar = M_.params( 1 );
M_.params( 2 ) = 0;
pi_bar = M_.params( 2 );
M_.params( 3 ) = 0;
r_bar = M_.params( 3 );
M_.params( 4 ) = 0.01;
kappa = M_.params( 4 );
M_.params( 5 ) = 3;
tau = M_.params( 5 );
M_.params( 6 ) = 1.5;
phi_pi = M_.params( 6 );
M_.params( 7 ) = 0.5;
phi_y = M_.params( 7 );
M_.params( 10 ) = 0;
rho_r = M_.params( 10 );
M_.params( 8 ) = 0.9;
rho_y = M_.params( 8 );
M_.params( 9 ) = 0.9;
rho_pi = M_.params( 9 );
M_.params( 11 ) = 0;
a_1 = M_.params( 11 );
M_.params( 12 ) = 0;
a_2 = M_.params( 12 );
M_.params( 13 ) = 0;
a_3 = M_.params( 13 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.7)^2;
M_.Sigma_e(2, 2) = (0.3)^2;
M_.Sigma_e(3, 3) = (0)^2;
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.param_vals = [estim_params_.param_vals; 1, 0.8, 0, 2, 3, 0.4, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, 0.7, (-2), 2, 2, 0.62, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 3, 0.9, (-2), 2, 2, 0.5, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, 0.1, 0, 1, 1, 0.3, 0.15, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, 5, 0, 10, 2, 2, 0.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, 1.5, 0, 10, 2, 1.5, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, 0.5, 0, 10, 2, 0.5, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, 0.5, 0, 1, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, 0.5, 0, 1, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 10, 0.5, 0, 1, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, 0.07, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, 1, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 3, 0.4, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 11, 0, (-5), 5, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 12, 0, (-5), 5, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 13, 0, (-5), 5, 3, 0, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 14, (-0.64), (-0.999), 0.999, 3, (-0.7), 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 17, (-0.02), (-0.999), 0.999, 3, (-0.02), 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 20, 0.73, (-0.999), 0.999, 3, 0.76, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 15, 5, (-20), 20, 3, 8.4, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 16, (-2), (-20), 20, 3, (-3.5), 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 18, 0.2, (-20), 20, 3, 0.48, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 19, 7, (-20), 20, 3, 7.8, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 21, 0.5, (-20), 20, 3, 0.56, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 22, 2, (-20), 20, 3, 2, 1, NaN, NaN, NaN ];
options_.mh_drop = 0.2;
options_.mh_jscale = 0.51;
options_.mh_nblck = 1;
options_.mh_replic = 0;
options_.mode_compute = 4;
options_.nodiagnostic = 1;
options_.nograph = 1;
options_.datafile = 'us_dataset';
options_.optim_opt = '''Algorithm'',''active-set''';
options_.first_obs = 44;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
options_.ar = 10;
options_.nograph = 1;
options_.periods = 10000;
var_list_ = char();
info = stoch_simul(var_list_);
save('NKPC_REE_Estimation_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('NKPC_REE_Estimation_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
