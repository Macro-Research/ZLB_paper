function dynare_estimation(varlist)
global options_ bayestopt_ estim_params_ lgy_ lgx_ Sigma_e_ exo_nbr endo_nbr
global lgy_ lgx_ fname_ oo_  trend_coeff_  iy_ ykmin_ lgx_TeX_ lgy_TeX_

%  The program was modified by Sergey Slobodyan and Raf Wouters to allow for
%  adaptive learning, both constant gain and Kalman filter
%  Sergey.Slobodyan@cerge-ei.cz, August 18, 2011

%disp('no hessian');
%disp('no prior graphs') % line 75
%disp('filtered shocks')

% Posterior maximization could be improved if every Projection Facility hit
% is penalized. Try to use 10 if maximization is repeatedly stuck in corner
% If maximization has found a good starting point, typically you won't need
% this parameter (set it equal to 0)
LIK_penalty = 0;

% Updating of expectations could start later than the Kalman filter (UPD_start > 0).
% Typically, we use simultaneous start as below
UPD_start = 0;

% This allows the Projection Facility to be used. Corr = 0 keeps
% track of PF hits but does not implement the PF. 
Corr = 1;

% Ridge regression is useful when the second moments matrix R has
% eigenvalues close to zero. It might happen when initial beliefs,
% including R, are derived from some REE. When information set less than
% MSV variables is used, typically this will never be used. 
% Ridge regression is mostly relevant for MSV constant gain learning. 
Ridge = 1.e-05;

% If this equal 1, ridge adjustment is used, no used if set to 0. For
% Kalman filter learning, Ridge and Ridge_adj are ignored
Ridge_adj = 1;

% Smoothing of the projection facility, somethimes useful
PJ_threshold = 0.995; % the point after which the projection facility starts

%%%%% CONSTANT GAIN LEARNING SETTINGS
%%%%%

% This variable includes all endogenous variables that are present on the
% RHS of the agents' forecasting equations.
% As given, 'states' and 'shocks' correspond to the MSV learning for the
% model.
states = {'c','inve','kp','pinf','r','w','a','al','b','g','ms','qs','spinf','sw'};

% It is assumed that the agents observe exogenous processes at time t, and 
% also know the transition matrix for the shocks. If the exogenous
% processes' names are included into the 'states' variable, the correlation
% coefficients will be, effectively, learned.
% If you don't want to use shocks in your learning scheme, use shocks = {};
% shocks = {'a','al','b','g','ms','qs','spinf','sw'};
shocks = {};

% Very important - this variable should include ALL forward-looking
% variables of the model. In other words, ALL expectations MUST be
% adaptive, it is impossible to have some of them remaining rational.
forwards = {'c','inve','lab','pk','rk','w','pinf'};

% THESE VARIABLES CONTAIN ALL THE MSV STATES AND SHOCKS; THEY ARE USED IN 
% SOME SPECIFICATIONS OF THE LEARNING ALGORITHM 
states_MSV = {'c','inve','kp','pinf','r','w','y'};
shocks_MSV = {'a','b','epinfma','ewma','g','ms','qs','spinf','sw'};

% EXAMPLES OF OTHER SETS OF VARIABLES FOR THE FORECASTING FUNCTIONS
% % This is the VAR LEARNING set
% states = {'c','inve','lab','pinf','y','r','w'};
% shocks = {};

%%%%% KALMAN FILTER LEARNING SETTINGS
%%%%%

% Number of different forecasting models the agents use for predicting values 
% of the forward-looking variables. Only one model in the baseline
% specification, five-model set-up is below, commented

NUM_MODELS = 1;

% A model is a cell array named 'm1', 'm2', and so on. In case only one
% model is present, it should be 'm1'.
% Every model describes every forecasting equation. As in constant-gain
% learning, all forward-looking variables have to be depicted in every
% forecasting model, otherwise the program will abort.

% Within every model, an equation is described by a block of the form
%           'inve',{'inve','invel'},{};

% The first part of the block - 'inve' -  denotes which equation is being 
% described, in this case an equation for investment variable. This part
% could be thought of as one piece of the 'forwards' variable in the
% constant gain case.

% The second part - {'inve','invel'} - describes endogenous states used on 
% the RHS of this forecasting equation. The part is analoguous to the 
% 'states' variable in constant gain case.

% The third part - {} - contains shocks (more generally, any variables
% whose values at time t are assumed to be known to the agents for
% producing their forecasts of time t+1 variables. Analoguous to the
% 'shocks' variable in the constant gain case. Currently this is a reserved
% option, please don't use any variables there.

% As written, 'm1' forecasting model describes seven individual equations
% in which every forward-looking variable is explained using two lags of
% itself. A constant will be added automatically if kalman_algo=603 is
% used. For no constant case Kalman filter learning, kalman_algo=503 is
% reserved (but currently not implemented).

% AR(2) Model
m1 = {'inve',{'inve','invel'},{};'c',{'c','cl'},{};'lab',{'lab','labl'},{};'pinf',{'pinf','pinfl'},{};...
    'pk',{'pk','pkl'},{};'rk',{'rk','rkl'},{};'w',{'w','wl'},{};};

% A STANDARD SET-UP WITH 5 MODELS

% % AR(1) Model
% m1 = {'inve',{'inve'},{};'c',{'c'},{};'lab',{'lab'},{};'w',{'w'},{};'pinf',{'pinf'},{};...
%     'pk',{'pk'},{};'rk',{'rk'},{};};
% % AR(1) + 2 ('pinf' and 'r') Model
% m2 = {'inve',{'inve','pinf','r'},{};'c',{'c','pinf','r'},{};'lab',{'lab','pinf','r'},{};'pinf',{'pinf','pinf','r'},{};...
%     'pk',{'pk','pinf','r'},{};'rk',{'rk','pinf','r'},{};'w',{'w','pinf','r'},{};};
% % AR(2) Model
% m3 = {'inve',{'inve','invel'},{};'c',{'c','cl'},{};'lab',{'lab','labl'},{};'pinf',{'pinf','pinfl'},{};...
%     'pk',{'pk','pkl'},{};'rk',{'rk','rkl'},{};'w',{'w','wl'},{};};
% % AR(1) + 1 ('pinf') Model
% m4 = {'inve',{'inve','pinf'},{};'c',{'c','pinf'},{};'lab',{'lab','pinf'},{};'pinf',{'pinf','pinf'},{};...
%     'pk',{'pk','pinf'},{};'rk',{'rk','pinf'},{};'w',{'w','pinf'},{};};
% % AR(1) + 3 model
% m5 = {'inve',{'inve','pinf','y','r'},{};'c',{'c','pinf','y','r'},{};'lab',{'lab','pinf','y','r'},{};'pinf',{'pinf','pinf','y','r'},{};...
%     'pk',{'pk','pinf','y','r'},{};'rk',{'rk','pinf','y','r'},{};'w',{'w','pinf','y','r'},{};};

% OTHER MISC MODELS
% % AR(1) + 7 Model
% m1 = {'inve',{'c','inve','lab','pinf','y','r','w'},{};'c',{'c','inve','lab','pinf','y','r','w'},{};...
%     'lab',{'c','inve','lab','pinf','y','r','w'},{};'pinf',{'c','inve','lab','pinf','y','r','w'},{};...
%     'pk',{'pk','c','inve','lab','pinf','y','r','w'},{};'rk',{'rk','c','inve','lab','pinf','y','r','w'},{};...
%     'w',{'c','inve','lab','pinf','y','r','w'},{};};

% % VAR(7) Model
% m1 = {'inve',{'c','inve','lab','pinf','y','r','w'},{};'c',{'c','inve','lab','pinf','y','r','w'},{};...
%     'lab',{'c','inve','lab','pinf','y','r','w'},{};'pinf',{'c','inve','lab','pinf','y','r','w'},{};...
%     'pk',{'c','inve','lab','pinf','y','r','w'},{};'rk',{'c','inve','lab','pinf','y','r','w'},{};...
%     'w',{'c','inve','lab','pinf','y','r','w'},{};};

% % AR(1) + 2 ('pinf' and 'y') Model
% % m4 = {'inve',{'inve','pinf','y'},{};'c',{'c','pinf','y'},{};'lab',{'lab','pinf','y'},{};'pinf',{'pinf','pinf','y'},{};...
% %     'pk',{'pk','pinf','y'},{};'rk',{'rk','pinf','y'},{};'w',{'w','pinf','y'},{};};

% % AR(1) + 1 ('pinf') Model
% m4 = {'inve',{'inve','pinf'},{};'c',{'c','pinf'},{};'lab',{'lab','pinf'},{};'pinf',{'pinf','pinf'},{};...
%     'pk',{'pk','pinf'},{};'rk',{'rk','pinf'},{};'w',{'w','pinf'},{};};

% This variable points to the file name which includes a parameter vector
% that is used to calculate REE for the 2nd way of generating initial
% beliefs. The file name should end in '_mode' which isn't part of the
% variable use_model
% The variable should be set for both constant gain (kalman_algo=102 and
% 202) and Kalman filter learning (kalman_algo=502 and 602)
use_model='MSV_REE_consistent_presample_beliefs';

% This variable contains the initial beliefs (including second moments) for
% the 4th way of forming initial beliefs (based on pre-sample regression). 
% Order of variables: Endogenous first, then shocks, alphabetical within groups.
% File name should end in '_bet' which isn't part of the use_bet variable.
% Pre-sample regression initial beliefs are not currently working for the
% Kalman filter learning
use_bet='complete_sample_19_mode';

options_.varlist = varlist;

options_.lgyidx2varobs = zeros(size(lgy_,1),1);
for i = 1:size(lgy_,1)
  tmp = strmatch(deblank(lgy_(i,:)),options_.varobs,'exact');
  if ~isempty(tmp)
    options_.lgyidx2varobs(i,1) = tmp;
  end			
end
options_ = set_default_option(options_,'kalman_algo',1);
options_ = set_default_option(options_,'lik_penalty',LIK_penalty);
options_ = set_default_option(options_,'upd_start',UPD_start);
options_ = set_default_option(options_,'corr',Corr);
options_ = set_default_option(options_,'ridge',Ridge);
options_ = set_default_option(options_,'ridge_adj',Ridge_adj);
% options_ = set_default_option(options_,'gsg',GSG);
options_ = set_default_option(options_,'pj_threshold',PJ_threshold);
options_ = set_default_option(options_,'first_obs',1);
options_ = set_default_option(options_,'prefilter',0);
options_ = set_default_option(options_,'presample',0);
options_ = set_default_option(options_,'lik_algo',1);
options_ = set_default_option(options_,'lik_init',1);
options_ = set_default_option(options_,'nograph',0);
options_ = set_default_option(options_,'mh_conf_sig',0.90);
options_ = set_default_option(options_,'mh_replic',20000);
options_ = set_default_option(options_,'mh_drop',0.5);
options_ = set_default_option(options_,'mh_jscale',0.2);
options_ = set_default_option(options_,'mh_init_scale',2*options_.mh_jscale);
options_ = set_default_option(options_,'mode_file','');
options_ = set_default_option(options_,'mode_compute',4);
options_ = set_default_option(options_,'mode_check',0);
options_ = set_default_option(options_,'prior_trunc',1e-10);
options_ = set_default_option(options_,'mh_mode',1); 	
options_ = set_default_option(options_,'mh_nblck',2);	
options_ = set_default_option(options_,'load_mh_file',0);
options_ = set_default_option(options_,'nodiagnostic',0);
options_ = set_default_option(options_,'loglinear',0);
options_ = set_default_option(options_,'unit_root_vars',[]);
options_ = set_default_option(options_,'XTick',[]);
options_ = set_default_option(options_,'XTickLabel',[]);
options_ = set_default_option(options_,'bayesian_irf',0);
options_ = set_default_option(options_,'bayesian_th_moments',0);
options_ = set_default_option(options_,'TeX',0);
options_ = set_default_option(options_,'irf',40);
options_ = set_default_option(options_,'relative_irf',0);
options_ = set_default_option(options_,'order',1);
options_ = set_default_option(options_,'ar',5);
options_ = set_default_option(options_,'dr_algo',0);
options_ = set_default_option(options_,'linear',1);
options_ = set_default_option(options_,'drop',0);
options_ = set_default_option(options_,'replic',1);
options_ = set_default_option(options_,'hp_filter',0);
options_ = set_default_option(options_,'forecast',0);
options_ = set_default_option(options_,'smoother',0);
options_ = set_default_option(options_,'moments_varendo',0);
options_ = set_default_option(options_,'filtered_vars',0);
options_ = set_default_option(options_,'kalman_algo',1);
options_ = set_default_option(options_,'kalman_tol',10^(-12));
options_ = set_default_option(options_,'diffuse_d',[]);
options_ = set_default_option(options_,'nk',4);
% an experimental option for Great Moderation - change in std between the
% subperiods
options_ = set_default_option(options_,'gm',1);

optim_options = optimset('display','iter','LargeScale','off',...
			 'MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
         
options_.save_print = 1;

if isfield(options_,'optim_opt')
  eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
end

pnames=['     ';'beta ';'gamm ';'norm ';'invg ';'unif ';'invg2'];
n_varobs = size(options_.varobs,1);
[xparam1,estim_params_,bayestopt_,lb,ub]=set_prior(estim_params_);
if any(bayestopt_.pshape > 0)
%     plot_priors
%     pause
else
  options_.mh_replic = 0;
end
bounds = prior_bounds(bayestopt_);
bounds(:,1)=max(bounds(:,1),lb);
bounds(:,2)=min(bounds(:,2),ub);
if any(xparam1 < bounds(:,1)) || any(xparam1 > bounds(:,2))
  find(xparam1 < bounds(:,1))
  find(xparam1 > bounds(:,2))
  error('Initial parameter values are outside parameter bounds')
end
lb = bounds(:,1);
ub = bounds(:,2);
bayestopt_.lb = lb;
bayestopt_.ub = ub;

if isempty(trend_coeff_)
  bayestopt_.with_trend = 0;
else
  bayestopt_.with_trend = 1;
  bayestopt_.trend_coeff = {};
  for i=1:n_varobs
    if i > length(trend_coeff_) || isempty(trend_coeff_{i})
      bayestopt_.trend_coeff{i} = '0';
    else
      bayestopt_.trend_coeff{i} = trend_coeff_{i};
    end
  end
end

bayestopt_.penalty=1e8;             % penalty 

nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np ;
nx = nvx+nvn+ncx+ncn+np;

check_model;

%static solver
if exist([fname_ '_steadystate'],'file')
  bayestopt_.static_solve = [fname_ '_steadystate'];
else
  bayestopt_.static_solve = 'dynare_solve';
end

dr = set_state_space([]);

tmp = deblank(lgy_(dr.order_var,:));

% Creation of lists for learning step. kalman_algo=1 means RE estimation.
% The variables are pre-sorted, thus 'states' and 'shocks' don't need to be
% put in alpha order above. Similarly, variable names within the models 
% (m1, m2, etc.) and equations are sorted automatically.
% NOTE that DYNARE 3 sorts variables found in the .mod file in one register
% (lower, I believe). If some of your variables combine upper and lower 
% letters and especially an '_' sign, this might present problem for the AL
% estimation. The best approach is to use only lower case letters in 
% variable names in the .mod file.
if options_.kalman_algo > 500
    forwards = sort(forwards);
    yf_list = zeros(length(forwards),1);
    for k = 1:length(yf_list)
        yf_list(k) = strmatch(deblank(forwards(k)),tmp,'exact');
    end
    
    m(NUM_MODELS) = struct;    
    for i = 1:NUM_MODELS
        eval(['m_t = m',deblank(int2str(i)),';']);
        f_list = m_t(1,1);
        for j = 2:length(forwards) % assume every model describes all forecasting equations
            f_list = [f_list;m_t(j,1)];
        end
        [f_list_sort,ind_f_sort] = sort(f_list);
        m(i).yf_list = f_list_sort;
        y_st_all = sort(m_t{ind_f_sort(1),2})';
        shock_all = sort(m_t{ind_f_sort(1),3})';
        for j = 1:length(forwards)
            jp = ind_f_sort(j);
            y_st  = sort(m_t{jp,2})';
            y_st_list = zeros(size(y_st));
            shock = sort(m_t{jp,3})';
            y_sh_list = zeros(size(shock));
            if j > 1
                y_st_all = [y_st_all;y_st];
                shock_all = [shock_all;shock];
            end
            for k = 1:length(y_st)
                y_st_list(k) = strmatch(deblank(y_st(k)),tmp,'exact');
            end
            for k = 1:length(shock)
                y_sh_list(k) = strmatch(deblank(shock(k)),tmp,'exact');
            end
            m(i).(['y_st_list_',int2str(j)]) = y_st_list;
            m(i).(['y_sh_list_',int2str(j)]) = y_sh_list;
        end
        [m(i).state_unique,m(i).st_iu,m(i).st_ju] = unique(deblank(y_st_all));
        [m(i).shock_unique,m(i).sh_iu,m(i).sh_ju] = unique(deblank(shock_all));
        m(i).y_st_full = zeros(length(m(i).state_unique),length(forwards));
        m(i).y_sh_full = zeros(length(m(i).shock_unique),length(forwards));
        offset_y = 0;
        offset_s = 0;
        for j = 1:length(forwards)
            len_y = length(m(i).(['y_st_list_',int2str(j)]));
            len_s = length(m(i).(['y_sh_list_',int2str(j)]));
            m(i).y_st_full(m(i).st_ju(offset_y+1:offset_y+len_y),j) = m(i).(['y_st_list_',int2str(j)]);
            m(i).y_sh_full(m(i).sh_ju(offset_s+1:offset_s+len_s),j) = m(i).(['y_sh_list_',int2str(j)]);
            offset_y = offset_y + len_y;
            offset_s = offset_s + len_s;
        end
    end

    % Creating a joint to all models list of RHS variables
    unique_joint_st = m(1).state_unique;
    unique_joint_sh = m(1).shock_unique;
    if NUM_MODELS > 1
        for i = 2:NUM_MODELS
            unique_joint_st = union(unique_joint_st,m(i).state_unique);
            unique_joint_sh = union(unique_joint_sh,m(i).shock_unique);
        end
        for i = 1:NUM_MODELS
            pos = zeros(length(m(i).state_unique),1);
            for j = 1:length(m(i).state_unique)
                pos(j) = strmatch(m(i).state_unique(j),unique_joint_st,'exact');
            end
            y_st_full_temp = zeros(length(unique_joint_st),length(forwards));
            y_st_full_temp(pos,:) = m(i).y_st_full;
            m(i).y_st_full = y_st_full_temp;
            pos = zeros(length(m(i).shock_unique),1);
            for j = 1:length(m(i).shock_unique)
                pos(j) = strmatch(m(i).shock_unique(j),unique_joint_sh,'exact');
            end
            y_sh_full_temp = zeros(length(unique_joint_sh),length(forwards));
            y_sh_full_temp(pos,:) = m(i).y_sh_full;
            m(i).y_sh_full = y_sh_full_temp;
        end
    end
    
    % Transforming joint lists (cell arrays of chars) into lists of numbers
    ys_list = zeros(length(unique_joint_st),1);
    shock_list = zeros(length(unique_joint_sh),1);
    for i = 1:length(ys_list)
        ys_list(i) = strmatch(unique_joint_st(i),tmp,'exact');
    end
    for i = 1:length(shock_list)
        shock_list(i) = strmatch(unique_joint_sh(i),tmp,'exact');
    end
    % creating reorder vector
    [B,reorder] = sort([ys_list; shock_list]);
    
    ys_list_MSV = [];
    shock_list_MSV = [];    

    options_.m = m;
    options_.num_mod = NUM_MODELS;
    options_.y_st = length(unique_joint_st);
    options_.y_sh = length(unique_joint_sh);
elseif options_.kalman_algo > 100
    % Creation of lists for learning step
    states_MSV = sort(states_MSV);
    shocks_MSV = sort(shocks_MSV);    
    states = sort(states);
    shocks = sort(shocks);
    forwards = sort(forwards);
    tmp = lgy_(dr.order_var,:);
    ys_list_MSV = zeros(length(states_MSV),1);
    shock_list_MSV = zeros(length(shocks_MSV),1);
    ys_list = zeros(length(states),1);
    yf_list = zeros(length(forwards),1);
    shock_list = zeros(length(shocks),1);
    for i = 1:size(lgy_,1)
        j = strmatch(deblank(tmp(i,:)),states_MSV,'exact');
        if ~isempty(j)
            ys_list_MSV(j) = i;
        end
        j = strmatch(deblank(tmp(i,:)),shocks_MSV,'exact');
        if ~isempty(j)
            shock_list_MSV(j) = i;
        end
        j = strmatch(deblank(tmp(i,:)),states,'exact');
        if ~isempty(j)
            ys_list(j) = i;
        end
        j = strmatch(deblank(tmp(i,:)),shocks,'exact');
        if ~isempty(j)
            shock_list(j) = i;
        end
        j = strmatch(deblank(tmp(i,:)),forwards,'exact');
        if ~isempty(j)
            yf_list(j) = i;
        end
    end
    % creating reorder vector
    [B,reorder] = sort([ys_list; shock_list]);
else
    ys_list_MSV = [];
    shock_list_MSV = [];    
    ys_list = [];
    shock_list = [];
    yf_list = [];
    reorder = [];
end

% AL-related objects are mostly transferred through options_ structure
options_.ys_list_MSV = ys_list_MSV;
options_.shock_list_MSV = shock_list_MSV;    
options_.ys_list = ys_list;
options_.shock_list = shock_list;
options_.yf_list = yf_list;
options_.reorder = reorder;
    
% Initialization with unit-root variables
if ~isempty(options_.unit_root_vars)
  n_ur = length(options_.unit_root_vars);
  i_ur = zeros(n_ur,1);
  for i=1:n_ur
    i1 = strmatch(options_.unit_root_vars{i},lgy_(dr.order_var,:),'exact');
    if isempty(i1)
      error('Undeclared variable in unit_root_vars statement')
    end
    i_ur(i) = i1;
  end
  if ykmin_ > 1
    l1 = flipud([cumsum(iy_(1:ykmin_-1,dr.order_var),1);ones(1, ...
						  endo_nbr)]);
    n1 = nnz(l1);
    bayestopt_.Pinf = zeros(n1,n1);
    l2 = find(l1');
    l3 = zeros(endo_nbr,ykmin_);
    l3(i_ur,:) = l1(:,i_ur)';
    l3 = l3(:);
    i_ur1 = find(l3(l2));
    i_stable = ones(endo_nbr,1);
    i_stable(i_ur) = zeros(n_ur,1);
    i_stable = find(i_stable);
    bayestopt_.Pinf(i_ur1,i_ur1) = diag(ones(1,length(i_ur1)));
    bayestopt_.i_var_stable = i_stable;
    l3 = zeros(endo_nbr,ykmin_);
    l3(i_stable,:) = l1(:,i_stable)';
    l3 = l3(:);
    bayestopt_.i_T_var_stable = find(l3(l2));
  else
    n1 = endo_nbr;
    bayestopt_.Pinf = zeros(n1,n1);
    bayestopt_.Pinf(i_ur,i_ur) = diag(ones(1,length(i_ur)));
    l1 = ones(endo_nbr,1);
    l1(i_ur,:) = zeros(length(i_ur),1);
    bayestopt_.i_T_var_stable = find(l1);
  end
  options_.lik_init = 3;
end % if ~isempty(options_.unit_root_vars)

if isempty(options_.datafile)
  error('ESTIMATION: datafile option is missing')
end

if isempty(options_.varobs)
  error('ESTIMATION: VAROBS is missing')
end

% if jscale isn't specified for an estimated parameter, use
% global option options_.jscale, set to 0.2, by default
k = find(isnan(bayestopt_.jscale));
bayestopt_.jscale(k) = options_.mh_jscale;

% read and demean data 

if exist(options_.datafile)
  instr = options_.datafile;
else
  instr = ['load ' options_.datafile];
end

instr;

eval(instr);

rawdata = [];
k = [];
k1 = [];
for i=1:n_varobs
  rawdata = [rawdata eval(deblank(options_.varobs(i,:)))];
  k = [k strmatch(deblank(options_.varobs(i,:)),lgy_(dr.order_var,:),'exact')];
  k1 = [k1 strmatch(deblank(options_.varobs(i,:)),lgy_, 'exact')];
end

bayestopt_.mf = k;
bayestopt_.mfys = k1;

options_ = set_default_option(options_,'nobs',size(rawdata,1)-options_.first_obs+1);
gend = options_.nobs;

rawdata = rawdata(options_.first_obs:options_.first_obs+gend-1,:);
if options_.loglinear == 1
  rawdata = log(rawdata);
end
if options_.prefilter == 1
  bayestopt_.mean_varobs = mean(rawdata,1);
  data = (rawdata-ones(gend,1)*bayestopt_.mean_varobs)';
else
  data = rawdata';
end

if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
	 ' transformation'])
end

if ~isempty(options_.mode_file)
  eval(['load ' options_.mode_file ';']');
end
 
% Initial beliefs. If kalman_algo has three digits (AL estimation) and ends 
% in 1 (pre-sample regression based beliefs) or 2 (beliefs based on a 
% pre-specified REE), the initial beliefs are initialized below. For 
% kalman_algo ending in 3 (REE-consistent beliefs), the beliefs are always
% calculated within the DsgeLikelihood.m
if options_.kalman_algo > 100
    if mod(options_.kalman_algo,100) == 1
        R_beta = []; % in older versions of the file, R_beta might not exist
        eval(['load ' use_bet '_bet']);
    elseif mod(options_.kalman_algo,100) == 2
        xparam_save = xparam1;
        hh_save = hh;
        eval(['load ' [use_model,'_mode'] ' xparam1 hh;']');
        [info,betamat,SecondMoments,R_beta] = BetaFromTR_102(xparam1,gend,data);
        xparam1 = xparam_save;
        hh = hh_save;
    elseif mod(options_.kalman_algo,100) == 3
        betamat = [];
        SecondMoments = [];        
        R_beta = [];
    end
else
    betamat = [];
    SecondMoments = [];
    R_beta = [];
end

% AL-related objects are mostly transferred through options_ structure
options_.betamat = betamat;
options_.SecondMoments = SecondMoments;
options_.R_beta = R_beta;

tic
initial_estimation_checks(xparam1,gend,data);
toc

% This block has been modified to allow for one-sided Hessian derivations
% All the optimization options should work without problems.
if options_.mode_compute > 0
    fh=str2func('DsgeLikelihood');
%    info = 1;
%    while info
        if options_.mode_compute == 1 
          [xparam1,fval,exitflag,output,lamdba,grad,hessian_fmincon] = ...
              fmincon(fh,xparam1,[],[],[],[],lb,ub,[],optim_options,gend,data);
          hh = hessian_fmincon;
        eval(['save ' fname_ '_mode1 xparam1 hh hessian_fmincon;']);
        elseif options_.mode_compute == 2
            %	asamin('set','maximum_cost_repeat',0);
            [fval,xparam1,grad,hessian_asamin,exitflag] = ...
            asamin('minimize','DsgeLikelihood',xparam1,lb,ub,- ...
                   ones(size(xparam1)),gend,data);   
        elseif options_.mode_compute == 3
            [xparam1,fval,exitflag] = fminunc(fh,xparam1,optim_options,gend,data);
        elseif options_.mode_compute == 4
            H0 = 1e-4*eye(nx);
            crit = 1e-7;
            nit = 1000;
            disp('nit=250');
            %    nit=20;
            verbose = 2;
            [fval,xparam1,grad,hessian_csminwel,itct,fcount,retcodehat] = ...
            csminwel('DsgeLikelihood',xparam1,H0,[],crit,nit,gend,data);
            disp(sprintf('Objective function at mode: %f',fval))
            disp(sprintf('Objective function at mode: %f',DsgeLikelihood(xparam1,gend,data)))
            hh = hessian_csminwel;
            eval(['save ' fname_ '_mode1 xparam1 hh hessian_csminwel;']);
        elseif options_.mode_compute == 5
            flag = 0;
            [xparam1, hh, gg, fval] = newrat('DsgeLikelihood',xparam1,[],[],flag,gend,data);
            eval(['save ' fname_ '_mode xparam1 hh gg fval;']);
        end
        % Hessian derivations has been modified. First, a slightly different 
        % procedure that returns a matrix, not a vector, is used. Second, 
        % the derivation  is 'one-sided': if 'cliffs' are observed along
        % some direction, the derivation uses only function evaluations to
        % the 'good side'.
        if options_.mode_compute ~= 5
%            options_.lik_penalty = 0;
            %      options_.ridge = 0;
            scale = 1;
            % standard DYNARE hessian calculation
%           [hh, info, xparam1_min, f_min] = hessian_min('DsgeLikelihood',xparam1,scale,gend,data);
%            hh = reshape(hh,nx,nx);
            % one-sided Hessian; no need to remove the penalty
            [hh, info, xparam1_min, f_min] = hessian_min_1S('DsgeLikelihood',xparam1,scale,gend,data);
%            fprintf('Scale = %8.4f \n',scale);
            if info && ((fval - f_min)  > optim_options.TolFun)
                fprintf('A better point found: fval= %10.6f \n',f_min);
                % In the line below, the minimum found by the optimizer is
                % replaced by the point found during Hessian evaluation,
                % which amounts to a grid search around the minimum point.
                % This point might fail to respect the bounds imposed on
                % parameters. Using this option might be helpful if
                % optimization results in a corner point.
%                xparam1 = xparam1_min;
%             else
%                 info = 0;
            end
            eval(['save ' fname_ '_mode xparam1* hh* ;']);
            options_.lik_penalty = LIK_penalty;
            options_.ridge = Ridge;
        end
%   end
    eval(['save ' fname_ '_mode xparam1 hh ;']);
else
%   gradient_check('DsgeLikelihood',xparam1,gend,data,betamat,SecondMoments,ys_list,shock_list,yf_list,reorder);
end

if options_.mode_check == 1
  mode_check(xparam1,0,hh,gend,data,lb,ub);
end

%%%% Hessian manipulation for MH
% Under AL, very ofter the Hessian would have several negative eigenvalues,
% no matter how fine is the grid in the Hessian derivation. Generalized
% Cholesky decomposition, used by DYNARE, results in close to zero standard
% errors for many parameters, making such Hessian useless for
% Metropolis-Hastings step.

% The manipulation makes all negative eigenvalues equal to the smallest
% positive one and preserves the same eigenvectors. The resulting matrix is
% commonly useful for proposal density in the Metropolis step.

% If manipulation has been performed, Laplace approximation is not valid.
% The larger is the most negative eigenvalue of the original Hessian (printed 
% by the hessian_min_1S procedure), the larger is the distortion.
[V,D] = eig(hh);
DD = diag(D);
neg_eigs = length(find(DD<0));
if neg_eigs > 0
    disp('HESSIAN MANIPULATED TO MAKE NEGATIVE EIGENVALUES EQUAL TO THE SMALLEST POSITIVE ONE')
    DD(DD<0) = DD(neg_eigs+1);
    hh = V * diag(DD) / V;
end
%%%% End of Hessian manipulation
hh = generalized_cholesky(hh);
invhess = inv(hh);
stdh = sqrt(diag(invhess));
if any(bayestopt_.pshape > 0)
  disp(' ')
  disp('RESULTS FROM POSTERIOR MAXIMIZATION')
  tstath = zeros(nx,1);
  for i = 1:nx
    tstath(i) = abs(xparam1(i))/stdh(i);
  end
  tit1 = sprintf('%10s %7s %8s %7s %6s %4s %6s\n',' ','prior mean', ...
		 'mode','s.d.','t-stat','prior','pstdev');
  if np
    ip = nvx+nvn+ncx+ncn+1;
    disp('parameters')
    disp(tit1)
    for i=1:np
      disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   deblank(estim_params_.param_names(i,:)), ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip)));
		   eval(['oo_.posterior_mode.parameters.' deblank(estim_params_.param_names(i,:)) ' = xparam1(ip);']);
		   eval(['oo_.posterior_std.parameters.' deblank(estim_params_.param_names(i,:)) ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  if nvx
    ip = 1;
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:nvx
		k = estim_params_.var_exo(i,1);
      	disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   deblank(lgx_(k,:)),bayestopt_.pmean(ip),xparam1(ip), ...
		   stdh(ip),tstath(ip),pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip))); 
      	Sigma_e_(k,k) = xparam1(ip)*xparam1(ip);
	   	eval(['oo_.posterior_mode.shocks_std.' deblank(lgx_(k,:)) ' = xparam1(ip);']);
		eval(['oo_.posterior_std.shocks_std.' deblank(lgx_(k,:)) ' = stdh(ip);']); 
	    ip = ip+1;
    end
  end
  if nvn
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = nvx+1;
    for i=1:nvn
    	disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', ...
		   deblank(options_.varobs(estim_params_.var_endo(i,1),: ...
					   )),bayestopt_.pmean(ip), ...
		   xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), ...
		   bayestopt_.pstdev(ip)));
	   	eval(['oo_.posterior_mode.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = xparam1(ip);']);
		eval(['oo_.posterior_std.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = stdh(ip);']); 
      	ip = ip+1;
    end
  end
  if ncx
    disp('correlation of shocks')
    disp(tit1)
    ip = nvx+nvn+1;
    for i=1:ncx
		k1 = estim_params_.corrx(i,1);
      	k2 = estim_params_.corrx(i,2);
      	name = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
      	disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip),  ...
		   pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
      	Sigma_e_(k1,k2) = xparam1(ip)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
      	Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
	   	eval(['oo_.posterior_mode.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = xparam1(ip);']);
		eval(['oo_.posterior_std.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = stdh(ip);']); 
      	ip = ip+1;
    end
  end
  if ncn
    disp('correlation of measurement errors')
    disp(tit1)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
    	k1 = estim_params_.corrn(i,1);
      	k2 = estim_params_.corrn(i,2);
      	name = [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))];
      	disp(sprintf('%12s %7.3f %8.4f %7.4f %7.4f %4s %6.4f', name, ...
		   bayestopt_.pmean(ip),xparam1(ip),stdh(ip),tstath(ip), ...
		   pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.pstdev(ip)));
	   	eval(['oo_.posterior_mode.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = xparam1(ip);']);
		eval(['oo_.posterior_std.measurement_errors_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = stdh(ip);']); 
      ip = ip+1;
    end
  end
  %%% Laplace approximation to the marginal log density: 
  poweradj=0;
  if det(invhess) == 0
      poweradj=2;
  end;    
  
  oo_.PosteriorMode = DsgeLikelihood(xparam1,gend,data);
  md_Laplace = .5*size(xparam1,1)*log(2*pi) + .5* (log(det(invhess*10^poweradj))-log(10^(poweradj*size(invhess,1)))) ...
      - oo_.PosteriorMode;  
  oo_.MarginalDensity.LaplaceApproximation = md_Laplace;    
  
  disp(' ')
  disp(sprintf('Log data density [Laplace approximation] is %f.',md_Laplace))
  disp(sprintf('Posterior Mode                           is %f.',oo_.PosteriorMode))
  disp(' ')
else
  disp(' ')
  disp('RESULTS FROM MAXIMUM LIKELIHOOD')
  tstath = zeros(nx,1);
  for i = 1:nx
    tstath(i) = abs(xparam1(i))/stdh(i);
  end
  tit1 = sprintf('%10s %10s %7s %6s\n',' ', ...
		 'Estimate','s.d.','t-stat');
  if np
    ip = nvx+nvn+ncx+ncn+1;
    disp('parameters')
    disp(tit1)
    for i=1:np
    	disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
		   deblank(estim_params_.param_names(i,:)), ...
		   xparam1(ip),stdh(ip),tstath(ip)));
	   	eval(['oo_.mle_mode.parameters.' deblank(estim_params_.param_names(i,:)) ' = xparam1(ip);']);
		eval(['oo_.mle_std.parameters.' deblank(estim_params_.param_names(i,:)) ' = stdh(ip);']); 
      	ip = ip+1;
    end
  end
  if nvx
    ip = 1;
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:nvx
	    k = estim_params_.var_exo(i,1);
      	disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
		   deblank(lgx_(k,:)),xparam1(ip), ...
		   stdh(ip),tstath(ip)));
      	Sigma_e_(k,k) = xparam1(ip)*xparam1(ip);
	   	eval(['oo_.mle_mode.shocks_std.' deblank(lgx_(k,:)) ' = xparam1(ip);']);
		eval(['oo_.mle_std.shocks_std.' deblank(lgx_(k,:)) ' = stdh(ip);']); 
	    ip = ip+1;
    end
  end
  if nvn
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = nvx+1;
    for i=1:nvn
      disp(sprintf('%12s %8.4f %7.4f %7.4f', ...
		   deblank(options_.varobs(estim_params_.var_endo(i,1),: ...
					   )), ...
		   xparam1(ip),stdh(ip),tstath(ip)))
	   	eval(['oo_.mle_mode.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = xparam1(ip);']);
		eval(['oo_.mle_std.measurement_errors_std.' deblank(options_.varobs(estim_params_.var_endo(i,1),:)) ' = stdh(ip);']);      
      ip = ip+1;
    end
  end
  if ncx
    disp('correlation of shocks')
    disp(tit1)
    ip = nvx+nvn+1;
    for i=1:ncx
      	k1 = estim_params_.corrx(i,1);
      	k2 = estim_params_.corrx(i,2);
      	name = [deblank(lgx_(k1,:)) ',' deblank(lgx_(k2,:))];
      	disp(sprintf('%12s %8.4f %7.4f %7.4f', name, ...
		   xparam1(ip),stdh(ip),tstath(ip)));
      	Sigma_e_(k1,k2) = xparam1(ip)*sqrt(Sigma_e_(k1,k1)*Sigma_e_(k2,k2));
      	Sigma_e_(k2,k1) = Sigma_e_(k1,k2);
	   	eval(['oo_.mle_mode.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = xparam1(ip);']);
		eval(['oo_.mle_std.shocks_corr.' deblank(lgx_(k1,:)) '_' deblank(lgx_(k2,:)) ' = stdh(ip);']);      
      	ip = ip+1;
    end
  end
  if ncn
    disp('correlation of measurement errors')
    disp(tit1)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      name = [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))];
      disp(sprintf('%12s %8.4f %7.4f %7.4f', name, ...
		   xparam1(ip),stdh(ip),tstath(ip)));
	   	eval(['oo_.mle_mode.measurement_error_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = xparam1(ip);']);
		eval(['oo_.mle_std.measurement_error_corr.' deblank(lgy_(k1,:)) '_' deblank(lgy_(k2,:)) ' = stdh(ip);']);      
      ip = ip+1;
    end
  end
end

if any(bayestopt_.pshape > 0) & options_.TeX %%%% Bayesian estimation (posterior mode) Latex output
  if np
    filename = [fname_ '_Posterior_Mode_1.TeX'];
    fidTeX = fopen(filename,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (parameters)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'{\\tiny \n')
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+ncx+ncn+1;
    for i=1:np
      fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      deblank(estim_params_.tex(i,:)),...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
	      bayestopt_.pmean(ip),...
	      estim_params_.param_vals(i,6),...
	      xparam1(ip),...
	      stdh(ip));
      ip = ip + 1;    
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (parameters)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:1}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'} \n')
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if nvx
    TeXfile = [fname_ '_Posterior_Mode_2.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of structural shocks)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'{\\tiny \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = 1;
    for i=1:nvx
      k = estim_params_.var_exo(i,1);
      fprintf(fidTeX,[ '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
	      deblank(lgx_TeX_(k,:)),...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
	      bayestopt_.pmean(ip),...
	      estim_params_.var_exo(i,7),...
	      xparam1(ip), ...
	      stdh(ip)); 
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of structural shocks)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:2}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'} \n')
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if nvn
    TeXfile = [fname_ '_Posterior_Mode_3.TeX'];
    fidTeX  = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of measurement errors)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+1;
    for i=1:nvn
      fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      deblank(options_.varobs_TeX(estim_params_.var_endo(i,1),:)), ...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...        
	      bayestopt_.pmean(ip), ...
	      estim_params_.var_endo(i,7),...        
	      xparam1(ip),...
	      stdh(ip)); 
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (standard deviation of measurement errors)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:3}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if ncx
    TeXfile = [fname_ '_Posterior_Mode_4.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of structural shocks)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+1;
    for i=1:ncx
      k1 = estim_params_.corrx(i,1);
      k2 = estim_params_.corrx(i,2);
      fprintf(fidTeX,[ '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
	      [deblank(lgx_TeX_(k1,:)) ',' deblank(lgx_TeX_(k2,:))], ...
	      deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
	      bayestopt_.pmean(ip), ...
	      estim_params_.corrx(i,8), ...
	      xparam1(ip), ...
	      stdh(ip));
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of structural shocks)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:4}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end
  if ncn
    TeXfile = [fname_ '_Posterior_Mode_5.TeX'];
    fidTeX = fopen(TeXfile,'w');
    fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
    fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of measurement errors)\n');
    fprintf(fidTeX,['%% ' datestr(now,0)]);
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,' \n');
    fprintf(fidTeX,'\\begin{table}\n');
    fprintf(fidTeX,'\\centering\n');
    fprintf(fidTeX,'\\begin{tabular}{l|lcccc} \n');
    fprintf(fidTeX,'\\hline\\hline \\\\ \n');
    fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n')
    fprintf(fidTeX,'\\hline \\\\ \n');
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
      k1 = estim_params_.corrn(i,1);
      k2 = estim_params_.corrn(i,2);
      fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
	      [deblank(lgy_(k1,:)) ',' deblank(lgy_(k2,:))], ...
	      pnames(bayestopt_.pshape(ip)+1,:), ...
	      bayestopt_.pmean(ip), ...
	      estim_params_.corrn(i,8), ...
	      xparam1(ip), ...
	      stdh(ip));
      ip = ip+1;
    end
    fprintf(fidTeX,'\\hline\\hline \n');
    fprintf(fidTeX,'\\end{tabular}\n ');    
    fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of measurement errors)}\n ');
    fprintf(fidTeX,'\\label{Table:Posterior:5}\n');
    fprintf(fidTeX,'\\end{table}\n');
    fprintf(fidTeX,'%% End of TeX file.\n');
    fclose(fidTeX);
  end 
end    

if (any(bayestopt_.pshape  >0 ) & options_.mh_replic) | (any(bayestopt_.pshape >0 ) & options_.load_mh_file)  % not ML estimation
  bounds = prior_bounds(bayestopt_);
  bayestopt_.lb = bounds(:,1);
  bayestopt_.ub = bounds(:,2);
  if any(xparam1 < bounds(:,1)) | any(xparam1 > bounds(:,2))
    find(xparam1 < bounds(:,1))
    find(xparam1 > bounds(:,2))
    error('Mode values are outside prior bounds. Reduce prior_trunc.')
  end
  metropolis(xparam1,invhess,gend,data,rawdata,bounds);
end

if ~((any(bayestopt_.pshape > 0) & options_.mh_replic) | (any(bayestopt_.pshape > 0) & options_.load_mh_file)) | ~options_.smoother  % ML estimation, or posterior mode without metropolis-hastings or metropolis without bayesian smooth variables
	options_.lik_algo = 2;

    % 'filteratt' is a variable that is true filtered vector, a_{t|t}, while standard DYNARE variable
    % denoted 'filtered_state_vector' here is, in fact, filtered prediction
    % a_{t+1|t}.
    
%   	[atT,innov,measurement_error,filtered_state_vector,forecast_vector,ys,trend_coeff,filteratt,bet] = ...
%         DsgeSmoother(xparam1,gend,data);
  	[atT,innov,measurement_error,filtered_state_vector,forecast_vector,model_forecast_vector,ys,trend_coeff,filteratt,bet] = ...
        DsgeSmoother(xparam1,gend,data);
  	for i=1:endo_nbr
    	eval(['oo_.SmoothedVariables.' deblank(lgy_(dr.order_var(i),:)) ' = atT(i,:)'';']);
        % ForecastVariables is what standard DYNARE calls
        % FilteredVariables!
    	eval(['oo_.ForecastVariables.' deblank(lgy_(dr.order_var(i),:)) ' = squeeze(forecast_vector(:,i,:))'';']);
    	eval(['oo_.FilteredVariables.' deblank(lgy_(dr.order_var(i),:)) ' = filteratt(i,:)'';']);
    end
    % It is possible to iterate forward small forecasting models in case
    % they are ARs. This will produce PLM-based multi-period forecasts of
    % forward-looking variables, see DiffuseKalmanSmoother*.m file.
    % The procedure is very ad-hoc, so its usage is not recommended

    %     for i = 1:length(yf_list)
    %     	eval(['oo_.Model_1_ForecastVariables.' deblank(lgy_(dr.order_var(yf_list(i)),:)) ' = squeeze(model_forecast_vector(1,:,i,:))'';']);
    %     	eval(['oo_.Model_2_ForecastVariables.' deblank(lgy_(dr.order_var(yf_list(i)),:)) ' = squeeze(model_forecast_vector(2,:,i,:))'';']);
    %     end
  	[nbplt,nr,nc,lr,lc,nstar] = pltorg(exo_nbr);
  	if options_.TeX
    	fidTeX = fopen([fname_ '_SmoothedShocks.TeX'],'w');
    	fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    	fprintf(fidTeX,' \n');
  	end    
  	if nbplt == 1
    	hh = figure('Name','Smoothed shocks');
    	NAMES = [];
    	if options_.TeX, TeXNAMES = [], end
    	for i=1:exo_nbr
      		subplot(nr,nc,i);
      		plot(1:gend,innov(i,:),'-k','linewidth',1)
      		hold on
      		plot([1 gend],[0 0],'-r','linewidth',.5)
      		hold off
      		xlim([1 gend])
      		name    = deblank(lgx_(i,:));
      		NAMES   = strvcat(NAMES,name);
      		if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
      		end
      		if options_.TeX
				texname = lgx_TeX_(i,1);
				TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      		end
      		title(name,'Interpreter','none')
      		eval(['oo_.SmoothedShocks.' deblank(lgx_(i,:)) ' = innov(i,:)'';']);
    	end
    	eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(1)]);
    	eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(1)]);
    	saveas(hh,[fname_ '_SmoothedShocks' int2str(1) '.fig']);
    	if options_.nograph, close(hh), end
    	if options_.TeX
      		fprintf(fidTeX,'\\begin{figure}[H]\n');
      		for jj = 1:exo_nbr
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      		end    
      		fprintf(fidTeX,'\\centering \n');
      		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(1));
      		fprintf(fidTeX,'\\caption{Smoothed shocks.}');
      		fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(1));
      		fprintf(fidTeX,'\\end{figure}\n');
      		fprintf(fidTeX,'\n');
      		fprintf(fidTeX,'%% End of TeX file.\n');
      		fclose(fidTeX);
    	end    
  	else
    	for plt = 1:nbplt-1
    		hh = figure('Name','Smoothed shocks');
    		set(0,'CurrentFigure',hh)
    		NAMES = [];
    		if options_.TeX, TeXNAMES = [], end
    		for i=1:nstar
				k = (plt-1)*nstar+i;
				subplot(nr,nc,i);
				plot([1 gend],[0 0],'-r','linewidth',.5)
				hold on
				plot(1:gend,innov(k,:),'-k','linewidth',1)
				hold off
				name = deblank(lgx_(k,:));
				NAMES = strvcat(NAMES,name);
				if ~isempty(options_.XTick)
	 				set(gca,'XTick',options_.XTick)
	 				set(gca,'XTickLabel',options_.XTickLabel)
				end
				xlim([1 gend])
				if options_.TeX
	  				texname = lgx_TeX_(k,:);
	  				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
				end    
				title(name,'Interpreter','none')
				eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
      		end
      		eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(plt)]);
      		eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(plt)]);
      		saveas(hh,[fname_ '_SmoothedShocks' int2str(plt) '.fig']);
      		if options_.nograph, close(hh), end
      		if options_.TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
				for jj = 1:nstar
	  				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
				end    
				fprintf(fidTeX,'\\centering \n');
				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(plt));
				fprintf(fidTeX,'\\caption{Smoothed shocks.}');
				fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(plt));
				fprintf(fidTeX,'\\end{figure}\n');
				fprintf(fidTeX,'\n');
      		end    
    	end
    	hh = figure('Name','Smoothed shocks');
    	set(0,'CurrentFigure',hh)
    	NAMES = [];
    	if options_.TeX, TeXNAMES = [], end
    	for i=1:exo_nbr-(nbplt-1)*nstar
      		k = (nbplt-1)*nstar+i;
      		if lr ~= 0
				subplot(lr,lc,i);
      		else
				subplot(nr,nc,i);
      		end    
      		plot([1 gend],[0 0],'-r','linewidth',0.5)
      		hold on
      		plot(1:gend,innov(k,:),'-k','linewidth',1)
      		hold off
      		name     = deblank(lgx_(k,:));
      		NAMES    = strvcat(NAMES,name);
      		if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
      		end
      		xlim([1 gend])
      		if options_.TeX
				texname  = lgx_TeX_(k,:);
				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
      		end
      		title(name,'Interpreter','none')
      		eval(['oo_.SmoothedShocks.' deblank(name) ' = innov(k,:)'';']);
    	end
    	eval(['print -depsc2 ' fname_ '_SmoothedShocks' int2str(nbplt)]);
    	eval(['print -dpdf ' fname_ '_SmoothedShocks' int2str(nbplt)]);
    	saveas(hh,[fname_ '_SmoothedShocks' int2str(nbplt) '.fig']);
    	if options_.nograph, close(hh), end
    	if options_.TeX
      		fprintf(fidTeX,'\\begin{figure}[H]\n');
      		for jj = 1:size(NAMES,1);
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
      		end    
      		fprintf(fidTeX,'\\centering \n');
      		fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedShocks%s}\n',fname_,int2str(nbplt));
      		fprintf(fidTeX,'\\caption{Smoothed shocks.}');
      		fprintf(fidTeX,'\\label{Fig:SmoothedShocks:%s}\n',int2str(nbplt));
      		fprintf(fidTeX,'\\end{figure}\n');
      		fprintf(fidTeX,'\n');
      		fprintf(fidTeX,'%% End of TeX file.\n');
      		fclose(fidTeX);
    	end    
  	end
  %%
  %%	Smooth observational errors...
  %%
  	yf = zeros(gend,n_varobs);
  	if options_.prefilter == 1
    	yf = atT(bayestopt_.mf,:)+repmat(mean_varobs',1,gend);
    	yfi = filtered_state_vector(bayestopt_.mf,2:end)+repmat(mean_varobs',1,gend);
  	elseif options_.loglinear == 1
    	yf = atT(bayestopt_.mf,:)+repmat(log(ys(bayestopt_.mfys)),1,gend)+...
	 		trend_coeff*[1:gend];
  	else
        % In AL with a constant (kalman_algo>200 and kalman_algo > 600), 
        % variables are modeled together with a constant, not relative to 
        % the trend as in standard DYNARE. For the AL without constant, the
        % same trend substraction as in DYNARE is used, therefore we need 
        % to use two different formulae here.
        % NOTE that the case with observational errors under AL is not 
        % included!
    	yf = atT(bayestopt_.mf,:)+repmat(ys(bayestopt_.mfys),1,gend)+...
	 		trend_coeff*[1:gend];
        if (options_.kalman_algo > 200 && options_.kalman_algo < 300) || ...
           (options_.kalman_algo > 600 && options_.kalman_algo < 700) 
            yfi = filtered_state_vector(bayestopt_.mf,1:end-1)+...
                trend_coeff*[1:gend];
        else
            yfi = filtered_state_vector(bayestopt_.mf,1:end-1)+repmat(ys(bayestopt_.mfys),1,gend)+...
                trend_coeff*[1:gend];
        end
            
  	end
  	if nvn
    	number_of_plots_to_draw = 0;
    	index = [];
    	for i=1:n_varobs
      		if max(abs(measurement_error(10:end))) > 0.000000001
				number_of_plots_to_draw = number_of_plots_to_draw + 1;
				index = cat(1,index,i);
      		end
      		eval(['oo_.SmoothedMeasurementErrors.' deblank(options_.varobs(i,:)) ...
	    			' = measurement_error(i,:)'';']);
    	end
    	[nbplt,nr,nc,lr,lc,nstar] = pltorg(number_of_plots_to_draw);
    	if options_.TeX
      		fidTeX = fopen([fname_ '_SmoothedObservationErrors.TeX'],'w');
      		fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
      		fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
      		fprintf(fidTeX,' \n');
    	end    
    	if nbplt == 1
      		hh = figure('Name','Smoothed observation errors');
      		set(0,'CurrentFigure',hh)
      		NAMES = [];
      		if options_.TeX, TeXNAMES = [], end
      		for i=1:number_of_plots_to_draw
				subplot(nr,nc,i);
				plot(1:gend,measurement_error(index(i),:),'-k','linewidth',1)
				hold on
				plot([1 gend],[0 0],'-r','linewidth',.5)
				hold off
				name    = deblank(options_.varobs(index(i),:));
				NAMES   = strvcat(NAMES,name);
				if ~isempty(options_.XTick)
	  				set(gca,'XTick',options_.XTick)
	  				set(gca,'XTickLabel',options_.XTickLabel)
				end
				if options_.TeX
	  				texname = options_.varobs_TeX(index(i),:);
	  				TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
				end
				title(name,'Interpreter','none')
      		end
      		eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(1)]);
      		eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(1)]);
      		saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(1) '.fig']);
      		if options_.nograph, close(hh), end
      		if options_.TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
				for jj = 1:number_of_plots_to_draw
	  				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
				end    
				fprintf(fidTeX,'\\centering \n');
				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',fname_,int2str(1));
				fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
				fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s',int2str(1));
				fprintf(fidTeX,'\\end{figure}\n');
				fprintf(fidTeX,'\n');
				fprintf(fidTeX,'%% End of TeX file.\n');
				fclose(fidTeX);
      		end    
    	else
      		for plt = 1:nbplt-1
				hh = figure('Name','Smoothed observation errors');
				set(0,'CurrentFigure',hh)
				NAMES = [];
				if options_.TeX, TeXNAMES = [], end
				for i=1:nstar
	  				k = (plt-1)*nstar+i;
	  				subplot(nr,nc,i);
	  				plot([1 gend],[0 0],'-r','linewidth',.5)
	  				hold on
	  				plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
	  				hold off
	  				name = deblank(options_.varobs(index(k),:));
	  				NAMES = strvcat(NAMES,name);
	  				if ~isempty(options_.XTick)
	    				set(gca,'XTick',options_.XTick)
	    				set(gca,'XTickLabel',options_.XTickLabel)
	  				end
	  				if options_.TeX
	    				texname = options_.varobs_TeX(k,:);
	    				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
	  				end    
	  				title(name,'Interpreter','none')
				end
				eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(plt)]);
				eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(plt)]);
				saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(plt) '.fig']);
				if options_.nograph, close(hh), end
				if options_.TeX
	  				fprintf(fidTeX,'\\begin{figure}[H]\n');
	  				for jj = 1:nstar
	    				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
	  				end    
	  				fprintf(fidTeX,'\\centering \n');
	  				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservationErrors%s}\n',fname_,int2str(plt));
	  				fprintf(fidTeX,'\\caption{Smoothed observation errors.}');
	  				fprintf(fidTeX,'\\label{Fig:SmoothedObservationErrors:%s}\n',int2str(plt));
	  				fprintf(fidTeX,'\\end{figure}\n');
	  				fprintf(fidTeX,'\n');
				end    
      		end
      		hh = figure('Name','Smoothed observation errors');
      		set(0,'CurrentFigure',hh)
      		NAMES = [];
      		if options_.TeX, TeXNAMES = [], end
      		for i=1:number_of_plots_to_draw-(nbplt-1)*nstar
				k = (nbplt-1)*nstar+i;
				if lr ~= 0
	  				subplot(lr,lc,i);
				else
	  				subplot(nr,nc,i);
				end    
				plot([1 gend],[0 0],'-r','linewidth',0.5)
				hold on
				plot(1:gend,measurement_error(index(k),:),'-k','linewidth',1)
				hold off
				name     = deblank(options_.varobs(index(k),:));
				NAMES    = strvcat(NAMES,name);
				if ~isempty(options_.XTick)
	  				set(gca,'XTick',options_.XTick)
	  				set(gca,'XTickLabel',options_.XTickLabel)
				end
				if options_.TeX
	  				texname  = options_.varobs_TeX(index(k),:);
	  				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
				end
				title(name,'Interpreter','none');
      		end
      		eval(['print -depsc2 ' fname_ '_SmoothedObservationErrors' int2str(nbplt)]);
      		eval(['print -dpdf ' fname_ '_SmoothedObservationErrors' int2str(nbplt)]);
      		saveas(hh,[fname_ '_SmoothedObservationErrors' int2str(nbplt) '.fig']);
      		if options_.nograph, close(hh), end
      		if options_.TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
				for jj = 1:size(NAMES,1);
	  				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
				end    
				fprintf(fidTeX,'\\centering \n');
				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_SmoothedObservedErrors%s}\n',fname_,int2str(nbplt));
				fprintf(fidTeX,'\\caption{Smoothed observed errors.}');
				fprintf(fidTeX,'\\label{Fig:SmoothedObservedErrors:%s}\n',int2str(nbplt));
				fprintf(fidTeX,'\\end{figure}\n');
				fprintf(fidTeX,'\n');
				fprintf(fidTeX,'%% End of TeX file.\n');
				fclose(fidTeX);
      		end    
    	end
    end	
%     %%
%     %%	Historical and smoothed variabes
%     %%
%     [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
%     if options_.TeX
%     	fidTeX = fopen([fname_ '_HistoricalAndSmoothedVariables.TeX'],'w');
%     	fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
%     	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
%     	fprintf(fidTeX,' \n');
%     end    
%     if nbplt == 1
%     	hh = figure('Name','Historical and smoothed variables');
%     	NAMES = [];
%     	if options_.TeX, TeXNAMES = [], end
%     	for i=1:n_varobs
% 			subplot(nr,nc,i);
% 			plot(1:gend,yf(i,:),'-r','linewidth',1)
% 			hold on
% 			plot(1:gend,rawdata(:,i),'-k','linewidth',1)
% 			hold off
% 			name    = deblank(options_.varobs(i,:));
% 			NAMES   = strvcat(NAMES,name);
% 			if ~isempty(options_.XTick)
% 				set(gca,'XTick',options_.XTick)
% 				set(gca,'XTickLabel',options_.XTickLabel)
% 			end
% 			xlim([1 gend])
% 			if options_.TeX
% 				texname = options_.varobs_TeX(i,1);
% 				TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
% 			end
% 			title(name,'Interpreter','none')
%     	end
%     	eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(1)]);
%     	eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(1)]);
%     	saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(1) '.fig']);
%     	if options_.nograph, close(hh), end
%     	if options_.TeX
% 			fprintf(fidTeX,'\\begin{figure}[H]\n');
% 			for jj = 1:n_varobs
% 				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
% 			end    
% 			fprintf(fidTeX,'\\centering \n');
% 			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(1));
% 			fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
% 			fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(1));
% 			fprintf(fidTeX,'\\end{figure}\n');
% 			fprintf(fidTeX,'\n');
% 			fprintf(fidTeX,'%% End of TeX file.\n');
% 			fclose(fidTeX);
%     	end    
%     else
%     	for plt = 1:nbplt-1
% 			hh = figure('Name','Historical and smoothed variables');
% 			set(0,'CurrentFigure',hh)
% 			NAMES = [];
% 			if options_.TeX, TeXNAMES = [], end
% 			for i=1:nstar
% 				k = (plt-1)*nstar+i;
% 				subplot(nr,nc,i);
% 				plot(1:gend,yf(k,:),'-r','linewidth',1)
% 				hold on
% 				plot(1:gend,rawdata(:,k),'-k','linewidth',1)
% 				hold off
% 				name = deblank(options_.varobs(k,:));
% 				NAMES = strvcat(NAMES,name);
% 				if ~isempty(options_.XTick)
% 					set(gca,'XTick',options_.XTick)
% 					set(gca,'XTickLabel',options_.XTickLabel)
% 				end
% 				xlim([1 gend])
% 				if options_.TeX
% 					texname = options_.varobs_TeX(k,:);
% 					TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
% 				end    
% 				title(name,'Interpreter','none')
% 			end
% 			eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(plt)]);
% 			eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(plt)]);
% 			saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(plt) '.fig']);
% 			if options_.nograph, close(hh), end
% 			if options_.TeX
% 				fprintf(fidTeX,'\\begin{figure}[H]\n');
% 				for jj = 1:nstar
% 					fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
% 				end    
% 				fprintf(fidTeX,'\\centering \n');
% 				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(plt));
% 				fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
% 				fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(plt));
% 				fprintf(fidTeX,'\\end{figure}\n');
% 				fprintf(fidTeX,'\n');
% 			end    
%     	end
%     	hh = figure('Name','Historical and smoothed variables');
%     	set(0,'CurrentFigure',hh)
%     	NAMES = [];
%     	if options_.TeX, TeXNAMES = [], end
%     	for i=1:n_varobs-(nbplt-1)*nstar
% 			k = (nbplt-1)*nstar+i;
% 			if lr ~= 0
% 				subplot(lr,lc,i);
% 			else
% 				subplot(nr,nc,i);
% 			end    
% 			plot(1:gend,yf(k,:),'-r','linewidth',1)
% 			hold on
% 			plot(1:gend,rawdata(:,k),'-k','linewidth',1)
% 			hold off
% 			name = deblank(options_.varobs(k,:));
% 			NAMES    = strvcat(NAMES,name);
% 			if ~isempty(options_.XTick)
% 				set(gca,'XTick',options_.XTick)
% 				set(gca,'XTickLabel',options_.XTickLabel)
% 			end
% 			xlim([1 gend])
% 			if options_.TeX
% 				texname  = options_.varobs_TeX(k,:);
% 				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
% 			end
% 			title(name,'Interpreter','none');
%     	end
%     	eval(['print -depsc2 ' fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
%     	eval(['print -dpdf ' fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt)]);
%     	saveas(hh,[fname_ '_HistoricalAndSmoothedVariables' int2str(nbplt) '.fig']);
%     	if options_.nograph, close(hh), end
%     	if options_.TeX
% 			fprintf(fidTeX,'\\begin{figure}[H]\n');
% 			for jj = 1:size(NAMES,1);
% 				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
% 			end    
% 			fprintf(fidTeX,'\\centering \n');
% 			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndSmoothedVariables%s}\n',fname_,int2str(nbplt));
% 			fprintf(fidTeX,'\\caption{Historical and smoothed variables.}');
% 			fprintf(fidTeX,'\\label{Fig:HistoricalAndSmoothedVariables:%s}\n',int2str(nbplt));
% 			fprintf(fidTeX,'\\end{figure}\n');
% 			fprintf(fidTeX,'\n');
% 			fprintf(fidTeX,'%% End of TeX file.\n');
% 			fclose(fidTeX);
%     	end    
% 	end
% end %	<--	if ML estimation, posterior mode without metropolis-hastings or metropolis 
% %		without bayesian posterior forecasts.
% 	%stoch_simul(lgy_);
% %end%	<--	if ML estimation, posterior mode without metropolis-hastings or metropolis 
% %		without bayesian posterior forecasts.
% 
% 
% 

    %%
    %%	Historical and Filtered variabes
    %%
    [nbplt,nr,nc,lr,lc,nstar] = pltorg(n_varobs);
    if options_.TeX
    	fidTeX = fopen([fname_ '_HistoricalAndFilteredVariables.TeX'],'w');
    	fprintf(fidTeX,'%% TeX eps-loader file generated by dynare_estimation.m (Dynare).\n');
    	fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    	fprintf(fidTeX,' \n');
    end    
    if nbplt == 1
    	hh = figure('Name','Historical and Filtered variables');
    	NAMES = [];
    	if options_.TeX, TeXNAMES = [], end
    	for i=1:n_varobs
			subplot(nr,nc,i);
			plot(1:gend,yfi(i,:),'-r','linewidth',1)
			hold on
			plot(1:gend,rawdata(:,i),'-k','linewidth',1)
			hold off
			name    = deblank(options_.varobs(i,:));
			NAMES   = strvcat(NAMES,name);
			if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
			end
			xlim([1 gend])
			if options_.TeX
				texname = options_.varobs_TeX(i,1);
				TeXNAMES   = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
			end
			title(name,'Interpreter','none')
    	end
    	eval(['print -depsc2 ' fname_ '_HistoricalAndFilteredVariables' int2str(1)]);
    	eval(['print -dpdf ' fname_ '_HistoricalAndFilteredVariables' int2str(1)]);
    	saveas(hh,[fname_ '_HistoricalAndFilteredVariables' int2str(1) '.fig']);
    	if options_.nograph, close(hh), end
    	if options_.TeX
			fprintf(fidTeX,'\\begin{figure}[H]\n');
			for jj = 1:n_varobs
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
			end    
			fprintf(fidTeX,'\\centering \n');
			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndFilteredVariables%s}\n',fname_,int2str(1));
			fprintf(fidTeX,'\\caption{Historical and filtered variables.}');
			fprintf(fidTeX,'\\label{Fig:HistoricalAndFilteredVariables:%s}\n',int2str(1));
			fprintf(fidTeX,'\\end{figure}\n');
			fprintf(fidTeX,'\n');
			fprintf(fidTeX,'%% End of TeX file.\n');
			fclose(fidTeX);
    	end    
    else
    	for plt = 1:nbplt-1
			hh = figure('Name','Historical and Filtered variables');
			set(0,'CurrentFigure',hh)
			NAMES = [];
			if options_.TeX, TeXNAMES = [], end
			for i=1:nstar
				k = (plt-1)*nstar+i;
				subplot(nr,nc,i);
				plot(1:gend,yfi(k,:),'-r','linewidth',1)
				hold on
				plot(1:gend,rawdata(:,k),'-k','linewidth',1)
				hold off
				name = deblank(options_.varobs(k,:));
				NAMES = strvcat(NAMES,name);
				if ~isempty(options_.XTick)
					set(gca,'XTick',options_.XTick)
					set(gca,'XTickLabel',options_.XTickLabel)
				end
				xlim([1 gend])
				if options_.TeX
					texname = options_.varobs_TeX(k,:);
					TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
				end    
				title(name,'Interpreter','none')
			end
			eval(['print -depsc2 ' fname_ '_HistoricalAndFilteredVariables' int2str(plt)]);
			eval(['print -dpdf ' fname_ '_HistoricalAndFilteredVariables' int2str(plt)]);
			saveas(hh,[fname_ '_HistoricalAndFilteredVariables' int2str(plt) '.fig']);
			if options_.nograph, close(hh), end
			if options_.TeX
				fprintf(fidTeX,'\\begin{figure}[H]\n');
				for jj = 1:nstar
					fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
				end    
				fprintf(fidTeX,'\\centering \n');
				fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndFilteredVariables%s}\n',fname_,int2str(plt));
				fprintf(fidTeX,'\\caption{Historical and Filtered variables.}');
				fprintf(fidTeX,'\\label{Fig:HistoricalAndFilteredVariables:%s}\n',int2str(plt));
				fprintf(fidTeX,'\\end{figure}\n');
				fprintf(fidTeX,'\n');
			end    
    	end
    	hh = figure('Name','Historical and Filtered variables');
    	set(0,'CurrentFigure',hh)
    	NAMES = [];
    	if options_.TeX, TeXNAMES = [], end
    	for i=1:n_varobs-(nbplt-1)*nstar
			k = (nbplt-1)*nstar+i;
			if lr ~= 0
				subplot(lr,lc,i);
			else
				subplot(nr,nc,i);
			end    
			plot(1:gend,yfi(k,:),'-r','linewidth',1)
			hold on
			plot(1:gend,rawdata(:,k),'-k','linewidth',1)
			hold off
			name = deblank(options_.varobs(k,:));
			NAMES    = strvcat(NAMES,name);
			if ~isempty(options_.XTick)
				set(gca,'XTick',options_.XTick)
				set(gca,'XTickLabel',options_.XTickLabel)
			end
			xlim([1 gend])
			if options_.TeX
				texname  = options_.varobs_TeX(k,:);
				TeXNAMES = strvcat(TeXNAMES,['$ ' deblank(texname) ' $']);
			end
			title(name,'Interpreter','none');
    	end
    	eval(['print -depsc2 ' fname_ '_HistoricalAndFilteredVariables' int2str(nbplt)]);
    	eval(['print -dpdf ' fname_ '_HistoricalAndFilteredVariables' int2str(nbplt)]);
    	saveas(hh,[fname_ '_HistoricalAndFilteredVariables' int2str(nbplt) '.fig']);
    	if options_.nograph, close(hh), end
    	if options_.TeX
			fprintf(fidTeX,'\\begin{figure}[H]\n');
			for jj = 1:size(NAMES,1);
				fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(jj,:)),deblank(TeXNAMES(jj,:)));
			end    
			fprintf(fidTeX,'\\centering \n');
			fprintf(fidTeX,'\\includegraphics[scale=0.5]{%s_HistoricalAndFilteredVariables%s}\n',fname_,int2str(nbplt));
			fprintf(fidTeX,'\\caption{Historical and Filtered variables.}');
			fprintf(fidTeX,'\\label{Fig:HistoricalAndFilteredVariables:%s}\n',int2str(nbplt));
			fprintf(fidTeX,'\\end{figure}\n');
			fprintf(fidTeX,'\n');
			fprintf(fidTeX,'%% End of TeX file.\n');
			fclose(fidTeX);
    	end    
	end
end %	<--	if ML estimation, posterior mode without metropolis-hastings or metropolis 
%		without bayesian posterior forecasts.
	%stoch_simul(lgy_);
%end%	<--	if ML estimation, posterior mode without metropolis-hastings or metropolis 
%		without bayesian posterior forecasts.







      
% SA 07-31-2004		* Added TeX output.
%					* Prior plots are done by calling plot_priors.m.
%					* All the computations related to the metropolis-hastings are made
%					in a new version of metropolis.m.
%					* Corrected a bug related to prior's bounds.
%					* ...
%					* If you do not want to see all the figures generated by dynare, you can use the option
%					nograph. The figures will be done and saved in formats eps, pdf and fig (so that you
%					should be able to modify the plots within matlab) but each figure will be erased from the
%					workspace when completed.
% SA 08-04-2004		Corrected a bug related to the display of the Smooth shocks and variables plots,
%					for ML and posterior mode estimation. 
% SA 09-03-2004		Compilation of TeX appendix moved to dynare.m.	                                                                                                                                                                       h