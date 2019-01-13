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
M_.fname = 'SW_BLE_Estimation_Mode';
M_.dynare_version = '4.5.4';
oo_.dynare_version = '4.5.4';
options_.dynare_version = '4.5.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('SW_BLE_Estimation_Mode.log');
M_.exo_names = 'eta_a';
M_.exo_names_tex = '\eta_a';
M_.exo_names_long = 'eta_a';
M_.exo_names = char(M_.exo_names, 'eta_b');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_b');
M_.exo_names_long = char(M_.exo_names_long, 'eta_b');
M_.exo_names = char(M_.exo_names, 'eta_g');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_g');
M_.exo_names_long = char(M_.exo_names_long, 'eta_g');
M_.exo_names = char(M_.exo_names, 'eta_i');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_i');
M_.exo_names_long = char(M_.exo_names_long, 'eta_i');
M_.exo_names = char(M_.exo_names, 'eta_r');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_r');
M_.exo_names_long = char(M_.exo_names_long, 'eta_r');
M_.exo_names = char(M_.exo_names, 'eta_p');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_p');
M_.exo_names_long = char(M_.exo_names_long, 'eta_p');
M_.exo_names = char(M_.exo_names, 'eta_w');
M_.exo_names_tex = char(M_.exo_names_tex, '\eta_w');
M_.exo_names_long = char(M_.exo_names_long, 'eta_w');
M_.endo_names = 'labobs';
M_.endo_names_tex = 'log(l^{obs})';
M_.endo_names_long = 'labobs';
M_.endo_names = char(M_.endo_names, 'robs');
M_.endo_names_tex = char(M_.endo_names_tex, 'd(log(r^{obs})) ');
M_.endo_names_long = char(M_.endo_names_long, 'robs');
M_.endo_names = char(M_.endo_names, 'pinfobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'd(log(\pi^{obs}))');
M_.endo_names_long = char(M_.endo_names_long, 'pinfobs');
M_.endo_names = char(M_.endo_names, 'dy');
M_.endo_names_tex = char(M_.endo_names_tex, ' d( log(y^{obs})) ');
M_.endo_names_long = char(M_.endo_names_long, 'dy');
M_.endo_names = char(M_.endo_names, 'dc');
M_.endo_names_tex = char(M_.endo_names_tex, ' d( log(c^{obs})) ');
M_.endo_names_long = char(M_.endo_names_long, 'dc');
M_.endo_names = char(M_.endo_names, 'dinve');
M_.endo_names_tex = char(M_.endo_names_tex, ' d( log(i^{obs})) ');
M_.endo_names_long = char(M_.endo_names_long, 'dinve');
M_.endo_names = char(M_.endo_names, 'dw');
M_.endo_names_tex = char(M_.endo_names_tex, ' d( log(w^{obs})) ');
M_.endo_names_long = char(M_.endo_names_long, 'dw');
M_.endo_names = char(M_.endo_names, 'mc');
M_.endo_names_tex = char(M_.endo_names_tex, 'mc');
M_.endo_names_long = char(M_.endo_names_long, 'mc');
M_.endo_names = char(M_.endo_names, 'zcap');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'zcap');
M_.endo_names = char(M_.endo_names, 'rk');
M_.endo_names_tex = char(M_.endo_names_tex, 'rk');
M_.endo_names_long = char(M_.endo_names_long, 'rk');
M_.endo_names = char(M_.endo_names, 'k_s');
M_.endo_names_tex = char(M_.endo_names_tex, 'ks');
M_.endo_names_long = char(M_.endo_names_long, 'k_s');
M_.endo_names = char(M_.endo_names, 'q');
M_.endo_names_tex = char(M_.endo_names_tex, 'q');
M_.endo_names_long = char(M_.endo_names_long, 'q');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'l');
M_.endo_names_tex = char(M_.endo_names_tex, 'l');
M_.endo_names_long = char(M_.endo_names_long, 'l');
M_.endo_names = char(M_.endo_names, 'pinf');
M_.endo_names_tex = char(M_.endo_names_tex, '\pi');
M_.endo_names_long = char(M_.endo_names_long, 'pinf');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'outputGap');
M_.endo_names_tex = char(M_.endo_names_tex, 'Gap');
M_.endo_names_long = char(M_.endo_names_long, 'outputGap');
M_.endo_names = char(M_.endo_names, 'eps_a');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_a');
M_.endo_names_long = char(M_.endo_names_long, 'eps_a');
M_.endo_names = char(M_.endo_names, 'eps_b');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_b');
M_.endo_names_long = char(M_.endo_names_long, 'eps_b');
M_.endo_names = char(M_.endo_names, 'eps_g');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_g');
M_.endo_names_long = char(M_.endo_names_long, 'eps_g');
M_.endo_names = char(M_.endo_names, 'eps_i');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_i');
M_.endo_names_long = char(M_.endo_names_long, 'eps_i');
M_.endo_names = char(M_.endo_names, 'eps_r');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_r');
M_.endo_names_long = char(M_.endo_names_long, 'eps_r');
M_.endo_names = char(M_.endo_names, 'eps_p');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_p');
M_.endo_names_long = char(M_.endo_names_long, 'eps_p');
M_.endo_names = char(M_.endo_names, 'eps_w');
M_.endo_names_tex = char(M_.endo_names_tex, '\eps_w');
M_.endo_names_long = char(M_.endo_names_long, 'eps_w');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'dx');
M_.endo_names_tex = char(M_.endo_names_tex, 'dx');
M_.endo_names_long = char(M_.endo_names_long, 'dx');
M_.endo_names = char(M_.endo_names, 'AUX_EXO_LAG_124_0');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_EXO\_LAG\_124\_0');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_EXO_LAG_124_0');
M_.endo_names = char(M_.endo_names, 'AUX_EXO_LAG_125_0');
M_.endo_names_tex = char(M_.endo_names_tex, 'AUX\_EXO\_LAG\_125\_0');
M_.endo_names_long = char(M_.endo_names_long, 'AUX_EXO_LAG_125_0');
M_.endo_partitions = struct();
M_.param_names = 'curv_w';
M_.param_names_tex = '\eps_w';
M_.param_names_long = 'curv_w';
M_.param_names = char(M_.param_names, 'curv_p');
M_.param_names_tex = char(M_.param_names_tex, '\eps_p');
M_.param_names_long = char(M_.param_names_long, 'curv_p');
M_.param_names = char(M_.param_names, 'G');
M_.param_names_tex = char(M_.param_names_tex, 'G');
M_.param_names_long = char(M_.param_names_long, 'G');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, '\delta');
M_.param_names_long = char(M_.param_names_long, 'delta');
M_.param_names = char(M_.param_names, 'phi_w');
M_.param_names_tex = char(M_.param_names_tex, '\phi_w');
M_.param_names_long = char(M_.param_names_long, 'phi_w');
M_.param_names = char(M_.param_names, 'l_bar');
M_.param_names_tex = char(M_.param_names_tex, '\bar{l}');
M_.param_names_long = char(M_.param_names_long, 'l_bar');
M_.param_names = char(M_.param_names, 'pi_bar');
M_.param_names_tex = char(M_.param_names_tex, '\bar{\pi}');
M_.param_names_long = char(M_.param_names_long, 'pi_bar');
M_.param_names = char(M_.param_names, 'beta_const');
M_.param_names_tex = char(M_.param_names_tex, '\bar{\beta}');
M_.param_names_long = char(M_.param_names_long, 'beta_const');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, '\alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'gamma_bar');
M_.param_names_tex = char(M_.param_names_tex, '\bar{\gamma}');
M_.param_names_long = char(M_.param_names_long, 'gamma_bar');
M_.param_names = char(M_.param_names, 'psi');
M_.param_names_tex = char(M_.param_names_tex, '\psi');
M_.param_names_long = char(M_.param_names_long, 'psi');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, '\phi');
M_.param_names_long = char(M_.param_names_long, 'phi');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, '\lambda');
M_.param_names_long = char(M_.param_names_long, 'lambda');
M_.param_names = char(M_.param_names, 'phi_p');
M_.param_names_tex = char(M_.param_names_tex, '\phi_p');
M_.param_names_long = char(M_.param_names_long, 'phi_p');
M_.param_names = char(M_.param_names, 'iota_w');
M_.param_names_tex = char(M_.param_names_tex, '\iota_w');
M_.param_names_long = char(M_.param_names_long, 'iota_w');
M_.param_names = char(M_.param_names, 'xi_w');
M_.param_names_tex = char(M_.param_names_tex, '\xi_w');
M_.param_names_long = char(M_.param_names_long, 'xi_w');
M_.param_names = char(M_.param_names, 'iota_p');
M_.param_names_tex = char(M_.param_names_tex, '\iota_p');
M_.param_names_long = char(M_.param_names_long, 'iota_p');
M_.param_names = char(M_.param_names, 'xi_p');
M_.param_names_tex = char(M_.param_names_tex, '\xi_p');
M_.param_names_long = char(M_.param_names_long, 'xi_p');
M_.param_names = char(M_.param_names, 'sigma_c');
M_.param_names_tex = char(M_.param_names_tex, '\sigma_c');
M_.param_names_long = char(M_.param_names_long, 'sigma_c');
M_.param_names = char(M_.param_names, 'sigma_l');
M_.param_names_tex = char(M_.param_names_tex, '\sigma_l');
M_.param_names_long = char(M_.param_names_long, 'sigma_l');
M_.param_names = char(M_.param_names, 'r_pi');
M_.param_names_tex = char(M_.param_names_tex, '\r_{\pi}');
M_.param_names_long = char(M_.param_names_long, 'r_pi');
M_.param_names = char(M_.param_names, 'r_dy');
M_.param_names_tex = char(M_.param_names_tex, 'r_{\delta y}');
M_.param_names_long = char(M_.param_names_long, 'r_dy');
M_.param_names = char(M_.param_names, 'r_y');
M_.param_names_tex = char(M_.param_names_tex, 'r_y');
M_.param_names_long = char(M_.param_names_long, 'r_y');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, '\rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, '\rho_a');
M_.param_names_long = char(M_.param_names_long, 'rho_a');
M_.param_names = char(M_.param_names, 'rho_b');
M_.param_names_tex = char(M_.param_names_tex, '\rho_b');
M_.param_names_long = char(M_.param_names_long, 'rho_b');
M_.param_names = char(M_.param_names, 'rho_p');
M_.param_names_tex = char(M_.param_names_tex, '\rho_p');
M_.param_names_long = char(M_.param_names_long, 'rho_p');
M_.param_names = char(M_.param_names, 'rho_w');
M_.param_names_tex = char(M_.param_names_tex, '\rho_w');
M_.param_names_long = char(M_.param_names_long, 'rho_w');
M_.param_names = char(M_.param_names, 'rho_i');
M_.param_names_tex = char(M_.param_names_tex, '\rho_i');
M_.param_names_long = char(M_.param_names_long, 'rho_i');
M_.param_names = char(M_.param_names, 'rho_r');
M_.param_names_tex = char(M_.param_names_tex, '\rho_r');
M_.param_names_long = char(M_.param_names_long, 'rho_r');
M_.param_names = char(M_.param_names, 'rho_ga');
M_.param_names_tex = char(M_.param_names_tex, '\rho_{\ga}');
M_.param_names_long = char(M_.param_names_long, 'rho_ga');
M_.param_names = char(M_.param_names, 'rho_g');
M_.param_names_tex = char(M_.param_names_tex, '\rho_g');
M_.param_names_long = char(M_.param_names_long, 'rho_g');
M_.param_names = char(M_.param_names, 'Mu_w');
M_.param_names_tex = char(M_.param_names_tex, '\mu_w');
M_.param_names_long = char(M_.param_names_long, 'Mu_w');
M_.param_names = char(M_.param_names, 'Mu_p');
M_.param_names_tex = char(M_.param_names_tex, '\mu_p');
M_.param_names_long = char(M_.param_names_long, 'Mu_p');
M_.param_names = char(M_.param_names, 'beta1_i');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_i');
M_.param_names_long = char(M_.param_names_long, 'beta1_i');
M_.param_names = char(M_.param_names, 'beta2_i');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_i');
M_.param_names_long = char(M_.param_names_long, 'beta2_i');
M_.param_names = char(M_.param_names, 'betar_i');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_i');
M_.param_names_long = char(M_.param_names_long, 'betar_i');
M_.param_names = char(M_.param_names, 'betac_i');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_i');
M_.param_names_long = char(M_.param_names_long, 'betac_i');
M_.param_names = char(M_.param_names, 'betainv_i');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_i');
M_.param_names_long = char(M_.param_names_long, 'betainv_i');
M_.param_names = char(M_.param_names, 'betay_i');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_i');
M_.param_names_long = char(M_.param_names_long, 'betay_i');
M_.param_names = char(M_.param_names, 'betapinf_i');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_i');
M_.param_names_long = char(M_.param_names_long, 'betapinf_i');
M_.param_names = char(M_.param_names, 'betaw_i');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_i');
M_.param_names_long = char(M_.param_names_long, 'betaw_i');
M_.param_names = char(M_.param_names, 'beta1_q');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_q');
M_.param_names_long = char(M_.param_names_long, 'beta1_q');
M_.param_names = char(M_.param_names, 'beta2_q');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_q');
M_.param_names_long = char(M_.param_names_long, 'beta2_q');
M_.param_names = char(M_.param_names, 'betar_q');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_q');
M_.param_names_long = char(M_.param_names_long, 'betar_q');
M_.param_names = char(M_.param_names, 'betac_q');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_q');
M_.param_names_long = char(M_.param_names_long, 'betac_q');
M_.param_names = char(M_.param_names, 'betainv_q');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_q');
M_.param_names_long = char(M_.param_names_long, 'betainv_q');
M_.param_names = char(M_.param_names, 'betay_q');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_q');
M_.param_names_long = char(M_.param_names_long, 'betay_q');
M_.param_names = char(M_.param_names, 'betapinf_q');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_q');
M_.param_names_long = char(M_.param_names_long, 'betapinf_q');
M_.param_names = char(M_.param_names, 'betaw_q');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_q');
M_.param_names_long = char(M_.param_names_long, 'betaw_q');
M_.param_names = char(M_.param_names, 'beta1_rk');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_rk');
M_.param_names_long = char(M_.param_names_long, 'beta1_rk');
M_.param_names = char(M_.param_names, 'beta2_rk');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_rk');
M_.param_names_long = char(M_.param_names_long, 'beta2_rk');
M_.param_names = char(M_.param_names, 'betar_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_rk');
M_.param_names_long = char(M_.param_names_long, 'betar_rk');
M_.param_names = char(M_.param_names, 'betac_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_rk');
M_.param_names_long = char(M_.param_names_long, 'betac_rk');
M_.param_names = char(M_.param_names, 'betainv_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_rk');
M_.param_names_long = char(M_.param_names_long, 'betainv_rk');
M_.param_names = char(M_.param_names, 'betay_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_rk');
M_.param_names_long = char(M_.param_names_long, 'betay_rk');
M_.param_names = char(M_.param_names, 'betapinf_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_rk');
M_.param_names_long = char(M_.param_names_long, 'betapinf_rk');
M_.param_names = char(M_.param_names, 'betaw_rk');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_rk');
M_.param_names_long = char(M_.param_names_long, 'betaw_rk');
M_.param_names = char(M_.param_names, 'beta1_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_pinf');
M_.param_names_long = char(M_.param_names_long, 'beta1_pinf');
M_.param_names = char(M_.param_names, 'beta2_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_pinf');
M_.param_names_long = char(M_.param_names_long, 'beta2_pinf');
M_.param_names = char(M_.param_names, 'betar_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betar_pinf');
M_.param_names = char(M_.param_names, 'betac_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betac_pinf');
M_.param_names = char(M_.param_names, 'betainv_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betainv_pinf');
M_.param_names = char(M_.param_names, 'betay_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betay_pinf');
M_.param_names = char(M_.param_names, 'betapinf_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betapinf_pinf');
M_.param_names = char(M_.param_names, 'betaw_pinf');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_pinf');
M_.param_names_long = char(M_.param_names_long, 'betaw_pinf');
M_.param_names = char(M_.param_names, 'beta1_c');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_c');
M_.param_names_long = char(M_.param_names_long, 'beta1_c');
M_.param_names = char(M_.param_names, 'beta2_c');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_c');
M_.param_names_long = char(M_.param_names_long, 'beta2_c');
M_.param_names = char(M_.param_names, 'betar_c');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_c');
M_.param_names_long = char(M_.param_names_long, 'betar_c');
M_.param_names = char(M_.param_names, 'betac_c');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_c');
M_.param_names_long = char(M_.param_names_long, 'betac_c');
M_.param_names = char(M_.param_names, 'betainv_c');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_c');
M_.param_names_long = char(M_.param_names_long, 'betainv_c');
M_.param_names = char(M_.param_names, 'betay_c');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_c');
M_.param_names_long = char(M_.param_names_long, 'betay_c');
M_.param_names = char(M_.param_names, 'betapinf_c');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_c');
M_.param_names_long = char(M_.param_names_long, 'betapinf_c');
M_.param_names = char(M_.param_names, 'betaw_c');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_c');
M_.param_names_long = char(M_.param_names_long, 'betaw_c');
M_.param_names = char(M_.param_names, 'beta1_l');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_l');
M_.param_names_long = char(M_.param_names_long, 'beta1_l');
M_.param_names = char(M_.param_names, 'beta2_l');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_l');
M_.param_names_long = char(M_.param_names_long, 'beta2_l');
M_.param_names = char(M_.param_names, 'betar_l');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_l');
M_.param_names_long = char(M_.param_names_long, 'betar_l');
M_.param_names = char(M_.param_names, 'betac_l');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_l');
M_.param_names_long = char(M_.param_names_long, 'betac_l');
M_.param_names = char(M_.param_names, 'betainv_l');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_l');
M_.param_names_long = char(M_.param_names_long, 'betainv_l');
M_.param_names = char(M_.param_names, 'betay_l');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_l');
M_.param_names_long = char(M_.param_names_long, 'betay_l');
M_.param_names = char(M_.param_names, 'betapinf_l');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_l');
M_.param_names_long = char(M_.param_names_long, 'betapinf_l');
M_.param_names = char(M_.param_names, 'betaw_l');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_l');
M_.param_names_long = char(M_.param_names_long, 'betaw_l');
M_.param_names = char(M_.param_names, 'beta1_w');
M_.param_names_tex = char(M_.param_names_tex, 'beta1\_w');
M_.param_names_long = char(M_.param_names_long, 'beta1_w');
M_.param_names = char(M_.param_names, 'beta2_w');
M_.param_names_tex = char(M_.param_names_tex, 'beta2\_w');
M_.param_names_long = char(M_.param_names_long, 'beta2_w');
M_.param_names = char(M_.param_names, 'betar_w');
M_.param_names_tex = char(M_.param_names_tex, 'betar\_w');
M_.param_names_long = char(M_.param_names_long, 'betar_w');
M_.param_names = char(M_.param_names, 'betac_w');
M_.param_names_tex = char(M_.param_names_tex, 'betac\_w');
M_.param_names_long = char(M_.param_names_long, 'betac_w');
M_.param_names = char(M_.param_names, 'betainv_w');
M_.param_names_tex = char(M_.param_names_tex, 'betainv\_w');
M_.param_names_long = char(M_.param_names_long, 'betainv_w');
M_.param_names = char(M_.param_names, 'betay_w');
M_.param_names_tex = char(M_.param_names_tex, 'betay\_w');
M_.param_names_long = char(M_.param_names_long, 'betay_w');
M_.param_names = char(M_.param_names, 'betapinf_w');
M_.param_names_tex = char(M_.param_names_tex, 'betapinf\_w');
M_.param_names_long = char(M_.param_names_long, 'betapinf_w');
M_.param_names = char(M_.param_names, 'betaw_w');
M_.param_names_tex = char(M_.param_names_tex, 'betaw\_w');
M_.param_names_long = char(M_.param_names_long, 'betaw_w');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 7;
M_.endo_nbr = 31;
M_.param_nbr = 90;
M_.orig_endo_nbr = 29;
M_.aux_vars(1).endo_index = 30;
M_.aux_vars(1).type = 3;
M_.aux_vars(1).orig_index = 6;
M_.aux_vars(1).orig_lead_lag = 0;
M_.aux_vars(2).endo_index = 31;
M_.aux_vars(2).type = 3;
M_.aux_vars(2).orig_index = 7;
M_.aux_vars(2).orig_lead_lag = 0;
options_.varobs = cell(1);
options_.varobs(1)  = {'dy'};
options_.varobs(2)  = {'dc'};
options_.varobs(3)  = {'dinve'};
options_.varobs(4)  = {'labobs'};
options_.varobs(5)  = {'pinfobs'};
options_.varobs(6)  = {'dw'};
options_.varobs(7)  = {'robs'};
options_.varobs_id = [ 4 5 6 1 3 7 2  ];
M_.Sigma_e = zeros(7, 7);
M_.Correlation_matrix = eye(7, 7);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('SW_BLE_Estimation_Mode_static');
erase_compiled_function('SW_BLE_Estimation_Mode_dynamic');
M_.orig_eq_nbr = 29;
M_.eq_nbr = 31;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 21;
 0 22;
 0 23;
 0 24;
 0 25;
 0 26;
 0 27;
 0 28;
 0 29;
 1 30;
 0 31;
 2 32;
 3 33;
 4 34;
 5 35;
 6 36;
 7 37;
 8 38;
 9 39;
 10 40;
 11 41;
 12 42;
 13 43;
 14 44;
 15 45;
 16 46;
 17 47;
 18 48;
 0 49;
 19 50;
 20 51;]';
M_.nstatic = 11;
M_.nfwrd   = 0;
M_.npred   = 20;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 20;
M_.ndynamic   = 20;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:7];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(31, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(7, 1);
M_.params = NaN(90, 1);
M_.NNZDerivatives = [137; 0; -1];
M_.params( 51 ) = 0.9;
beta1_rk = M_.params( 51 );
M_.params( 43 ) = 0.9;
beta1_q = M_.params( 43 );
M_.params( 67 ) = 0.9;
beta1_c = M_.params( 67 );
M_.params( 35 ) = 0.9;
beta1_i = M_.params( 35 );
M_.params( 75 ) = 0;
beta1_l = M_.params( 75 );
M_.params( 59 ) = 0.9;
beta1_pinf = M_.params( 59 );
M_.params( 83 ) = 0.9;
beta1_w = M_.params( 83 );
M_.params( 52 ) = 0;
beta2_rk = M_.params( 52 );
M_.params( 44 ) = 0;
beta2_q = M_.params( 44 );
M_.params( 68 ) = 0;
beta2_c = M_.params( 68 );
M_.params( 36 ) = 0;
beta2_i = M_.params( 36 );
M_.params( 76 ) = 0;
beta2_l = M_.params( 76 );
M_.params( 60 ) = 0;
beta2_pinf = M_.params( 60 );
M_.params( 84 ) = 0;
beta2_w = M_.params( 84 );
M_.params( 53 ) = 0;
betar_rk = M_.params( 53 );
M_.params( 45 ) = 0;
betar_q = M_.params( 45 );
M_.params( 69 ) = 0;
betar_c = M_.params( 69 );
M_.params( 37 ) = 0;
betar_i = M_.params( 37 );
M_.params( 77 ) = 0;
betar_l = M_.params( 77 );
M_.params( 61 ) = 0;
betar_pinf = M_.params( 61 );
M_.params( 85 ) = 0;
betar_w = M_.params( 85 );
M_.params( 54 ) = 0;
betac_rk = M_.params( 54 );
M_.params( 46 ) = 0;
betac_q = M_.params( 46 );
M_.params( 70 ) = 0;
betac_c = M_.params( 70 );
M_.params( 38 ) = 0;
betac_i = M_.params( 38 );
M_.params( 78 ) = 0;
betac_l = M_.params( 78 );
M_.params( 62 ) = 0;
betac_pinf = M_.params( 62 );
M_.params( 86 ) = 0;
betac_w = M_.params( 86 );
M_.params( 55 ) = 0;
betainv_rk = M_.params( 55 );
M_.params( 47 ) = 0;
betainv_q = M_.params( 47 );
M_.params( 71 ) = 0;
betainv_c = M_.params( 71 );
M_.params( 39 ) = 0;
betainv_i = M_.params( 39 );
M_.params( 79 ) = 0;
betainv_l = M_.params( 79 );
M_.params( 63 ) = 0;
betainv_pinf = M_.params( 63 );
M_.params( 87 ) = 0;
betainv_w = M_.params( 87 );
M_.params( 56 ) = 0;
betay_rk = M_.params( 56 );
M_.params( 48 ) = 0;
betay_q = M_.params( 48 );
M_.params( 72 ) = 0;
betay_c = M_.params( 72 );
M_.params( 40 ) = 0;
betay_i = M_.params( 40 );
M_.params( 80 ) = 0;
betay_l = M_.params( 80 );
M_.params( 64 ) = 0;
betay_pinf = M_.params( 64 );
M_.params( 88 ) = 0;
betay_w = M_.params( 88 );
M_.params( 57 ) = 0;
betapinf_rk = M_.params( 57 );
M_.params( 49 ) = 0;
betapinf_q = M_.params( 49 );
M_.params( 73 ) = 0;
betapinf_c = M_.params( 73 );
M_.params( 41 ) = 0;
betapinf_i = M_.params( 41 );
M_.params( 81 ) = 0;
betapinf_l = M_.params( 81 );
M_.params( 65 ) = 0;
betapinf_pinf = M_.params( 65 );
M_.params( 89 ) = 0;
betapinf_w = M_.params( 89 );
M_.params( 58 ) = 0;
betaw_rk = M_.params( 58 );
M_.params( 50 ) = 0;
betaw_q = M_.params( 50 );
M_.params( 74 ) = 0;
betaw_c = M_.params( 74 );
M_.params( 42 ) = 0;
betaw_i = M_.params( 42 );
M_.params( 82 ) = 0;
betaw_l = M_.params( 82 );
M_.params( 66 ) = 0;
betaw_pinf = M_.params( 66 );
M_.params( 90 ) = 0;
betaw_w = M_.params( 90 );
M_.params( 4 ) = 0.025;
delta = M_.params( 4 );
M_.params( 5 ) = 1.5;
phi_w = M_.params( 5 );
M_.params( 3 ) = 0.18;
G = M_.params( 3 );
M_.params( 2 ) = 10;
curv_p = M_.params( 2 );
M_.params( 1 ) = 10;
curv_w = M_.params( 1 );
M_.params( 12 ) = 4.3869;
phi = M_.params( 12 );
M_.params( 19 ) = 1.0249;
sigma_c = M_.params( 19 );
M_.params( 13 ) = 0.7922;
lambda = M_.params( 13 );
M_.params( 16 ) = 0.6673;
xi_w = M_.params( 16 );
M_.params( 20 ) = 2.3583;
sigma_l = M_.params( 20 );
M_.params( 18 ) = 0.7389;
xi_p = M_.params( 18 );
M_.params( 15 ) = 0;
iota_w = M_.params( 15 );
M_.params( 17 ) = 0;
iota_p = M_.params( 17 );
M_.params( 11 ) = 0.4552;
psi = M_.params( 11 );
M_.params( 14 ) = 1.5003;
phi_p = M_.params( 14 );
M_.params( 22 ) = 0.1;
r_dy = M_.params( 22 );
M_.params( 24 ) = 0.5;
rho = M_.params( 24 );
M_.params( 7 ) = 0.5745;
pi_bar = M_.params( 7 );
M_.params( 8 ) = 0.1756;
beta_const = M_.params( 8 );
M_.params( 6 ) = 1.8140;
l_bar = M_.params( 6 );
M_.params( 10 ) = 0.4168;
gamma_bar = M_.params( 10 );
M_.params( 9 ) = 0.1325;
alpha = M_.params( 9 );
M_.params( 25 ) = 0.9743;
rho_a = M_.params( 25 );
M_.params( 26 ) = 0.3339;
rho_b = M_.params( 26 );
M_.params( 32 ) = 0.9860;
rho_g = M_.params( 32 );
M_.params( 29 ) = 0.4576;
rho_i = M_.params( 29 );
M_.params( 30 ) = 0.3301;
rho_r = M_.params( 30 );
M_.params( 27 ) = 0.0981;
rho_p = M_.params( 27 );
M_.params( 28 ) = 0.0438;
rho_w = M_.params( 28 );
M_.params( 34 ) = 0;
Mu_p = M_.params( 34 );
M_.params( 33 ) = 0;
Mu_w = M_.params( 33 );
M_.params( 31 ) = 0;
rho_ga = M_.params( 31 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.4084)^2;
M_.Sigma_e(2, 2) = (0.5541)^2;
M_.Sigma_e(3, 3) = (0.3977)^2;
M_.Sigma_e(4, 4) = (1.2428)^2;
M_.Sigma_e(5, 5) = (0.1075)^2;
M_.Sigma_e(6, 6) = (0.1871)^2;
M_.Sigma_e(7, 7) = (0.8479)^2;
M_.params( 15 ) = 0;
iota_w = M_.params( 15 );
M_.params( 17 ) = 0;
iota_p = M_.params( 17 );
M_.params( 13 ) = 0;
lambda = M_.params( 13 );
M_.params( 34 ) = 0;
Mu_p = M_.params( 34 );
M_.params( 33 ) = 0;
Mu_w = M_.params( 33 );
M_.params( 31 ) = 0;
rho_ga = M_.params( 31 );
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.param_vals = [estim_params_.param_vals; 12, 2.5, 1.1, 15, 3, 4, 1.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 19, 1.5, 0.25, 3, 3, 1.50, 0.375, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 13, 0.5, 0.001, 0.99, 1, 0.7, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 16, 0.57, 0.3, 0.99, 1, 0.75, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 20, 0.96, 0.25, 10, 3, 2, 0.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 18, 0.9, 0.01, 0.99, 1, 0.75, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 11, 0.3, 0.01, 1, 1, 0.5, 0.15, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 14, 1.3, 1.0, 3, 3, 1.25, 0.125, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 21, 1.2, 1.0, 3, 3, 1.5, 0.25, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 24, 0.8258, 0.5, 0.975, 1, 0.75, 0.10, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 23, 0.13, 0.001, 0.5, 3, 0.125, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 22, 0.12, 0.001, 0.5, 3, 0.125, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, 0.7, 0.1, 2.0, 2, 0.625, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, 0.19, 0.01, 2.0, 2, 0.25, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, (-2.25), (-10.0), 10.0, 3, 0.0, 2.0, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 10, 0.4, 0.1, 0.8, 3, 0.4, 0.10, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, 0.16, 0.01, 1.0, 3, 0.3, 0.05, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 25, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 26, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 32, 0.5, 0, .9999, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 29, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 30, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 27, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 28, 0.5, 0, .99, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, 0.55, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, 0.77, 0.025, 5, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 3, 0.55, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 4, 1.3, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 5, 0.2, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 6, 0.4, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 7, 0.3, 0.01, 3, 4, 0.1, 2, NaN, NaN, NaN ];
options_.kalman_algo = 1;
options_.lik_init = 1;
options_.mh_drop = 0.2;
options_.mh_jscale = 0.35;
options_.mh_nblck = 1;
options_.mh_replic = 00000;
options_.mode_compute = 1;
options_.nodiagnostic = 1;
options_.nograph = 1;
options_.prefilter = 0;
options_.presample = 4;
options_.datafile = 'raf_dataset';
options_.optim_opt = '''MaxIter'',200,''Algorithm'',''active-set''';
options_.first_obs = 147;
options_.nobs = 60;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
options_.ar = 10;
options_.irf = 100;
var_list_ = char();
info = stoch_simul(var_list_);
save('SW_BLE_Estimation_Mode_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('SW_BLE_Estimation_Mode_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
