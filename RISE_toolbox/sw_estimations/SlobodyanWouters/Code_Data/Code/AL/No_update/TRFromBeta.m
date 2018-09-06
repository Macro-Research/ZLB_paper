function [T,R] = TRFromBeta(beta,T)

%  This program was written for Sergey Slobodyan and Raf Wouters DSGE
%  estimation under adaptive learning project.
%  It takes agents' beliefs and uses them to generate a new transmission 
%  mechanism of the model. 
%  Both constant gain and Kalman filter learning are supported.
%  Sergey.Slobodyan@cerge-ei.cz, August 18, 2011

global dr_ estim_params_ options_

%if options_.kalman_algo > 200
%    betamat(1,:) = [];
%end

ys_list = options_.ys_list;
shock_list = options_.shock_list;
reorder = options_.reorder;

% Kalman filter learning piece. This subroutine is used only BEFORE the
% estimation begins, when all weights are initialized at 1/NUM_MODELS. 
% Therefore, aggregate beliefs are a simple average of individual models' beliefs.  

if options_.kalman_algo > 500
    betama = zeros(size(options_.m(1).y_st_full,1),size(options_.m(1).y_st_full,2),options_.num_mod);
    for i = 1:options_.num_mod
        betam = zeros(size(options_.m(i).y_st_full));
        betam = options_.m(i).y_st_full;
        betam(betam ~= 0) = options_.m(i).beta;
        betama(:,:,i) = betam;
    end
    betamat = sum(betama,3) / options_.num_mod;
else
    betamat = beta;
end

n_y = length(ys_list);
vars = dr_.order_var(sort([ys_list; shock_list]));
fw = size(dr_.jacobian,2) - estim_params_.nvx;
t_vars = dr_.npred+1:fw-dr_.nsfwrd;

% Re-constructing the matrix which gives us E_t[y_{t+1}]. Columns are
% still in ENDOGENOUS, SHOCKS order
RO = T(shock_list,shock_list);
betama_ = [betamat(1:n_y,:)' betamat(n_y+1:end,:)' * RO];
% Transforming columns back into order_var order
betama_ = betama_(:,reorder);
% Generating predictions as a function of current vars, and padding it
% to the square matrix.
EXP = zeros(size(dr_.jacobian,1));
EXP(:,vars) = dr_.jacobian(:,fw-dr_.nsfwrd+1:fw) * betama_;
% Taking the lagged part of jacobia_ into full matrix. Order of
% variables is alpha throughout. In order_var, they are alpha in two
% groups separetely, that's why need sorting here.
LG = zeros(size(dr_.jacobian,1));
LG(:,sort(dr_.order_var(dr_.nstatic+1:dr_.nstatic+dr_.npred))) = ...
  dr_.jacobian(:,1:dr_.npred);
% Transition equations for purely backward-looking model with learning
T = - (dr_.jacobian(:,t_vars) + EXP) \ LG;
R = - (dr_.jacobian(:,t_vars) + EXP) \ dr_.jacobian(:,fw+1:end);
% Transformation of T: both rows and columns should be put in
% order_var. R: rows should be put in order_var.
T = T(dr_.order_var,dr_.order_var);
R = R(dr_.order_var,:);
