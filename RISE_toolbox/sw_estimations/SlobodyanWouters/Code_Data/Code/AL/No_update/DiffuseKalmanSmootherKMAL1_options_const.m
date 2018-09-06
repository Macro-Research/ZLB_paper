function [alphahat,etahat,a, aK, aKM, filteratt,bet] = DiffuseKalmanSmootherKMAL1_options_const(T,R,Q,Pinf,Pstar,Y,SState,start,np,smpl,mf,endpars)

% modified by M. Ratto:
% new output argument aK (1-step to k-step predictions)
% new options_.nk: the max step ahed prediction in aK (default is 4)
% new crit1 value for rank of Pinf
% it is assured that P is symmetric 
%
% stephane.adjemian@cepremap.cnrs.fr [09-16-2004]
% 
%   See "Filtering and Smoothing of State Vector for Diffuse State Space
%   Models", S.J. Koopman and J. Durbin (2003, in Journal of Time Series 
%   Analysis, vol. 24(1), pp. 85-98).  
%   calculations of smoothed vars are probably still incorrect - Sergey

%  The program was modified by Sergey Slobodyan and Raf Wouters to allow for
%  adaptive learning
%  Sergey.Slobodyan@cerge-ei.cz, August 18, 2011

global bayestopt_ options_  dr_ estim_params_ lgy_ oo_
% gain = 0.030969052313799; ro = endpars; sigma = 0.003;
% gain = endpars(1); ro = 1; sigma = endpars(2);

% Beliefs updating is switched off
gain = 0.030969052313799; ro = 0; sigma = 0.003;

% As a rule, only features specific to DiffuseKalmanSmoother*.m are documented
% here, the rest is identical with DiffuseLikelihood*.m

% this should be set to 0 if I intend to use MCMC to generate posterior
% distributions of forecast variables etc.
save_print = options_.save_print;

ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;
reorder = options_.reorder;  
upd_start = options_.upd_start;
lik_penalty = options_.lik_penalty;
YF = length(yf_list);
ZZ = length(shock_list) + length(ys_list);
NUM_MODELS = options_.num_mod;
df = zeros(NUM_MODELS,1);
dR = zeros(NUM_MODELS,1);
BIC = ones(NUM_MODELS,1);
% AIC = ones(NUM_MODELS,1);
weights = ones(NUM_MODELS,1) / NUM_MODELS;
err = zeros(length(yf_list),NUM_MODELS);
leave = [1:4,7];
Sig_joint = zeros(length(yf_list));
weight = zeros(NUM_MODELS,smpl);
weigth(:,1) = weights;
errors = zeros(length(yf_list),NUM_MODELS+1,smpl);
TT  = zeros(size(T,1),size(T,2),smpl+1);
TT(:,:,1) = T;
RR  = zeros(size(R,1),size(R,2),smpl+1);
RR(:,:,1) = R;
mus = zeros(size(T,1),smpl+1);
exps = zeros(length(yf_list),smpl);
tS = zeros(length(leave),length(leave),NUM_MODELS);
means = SState(dr_.order_var);
tmp_ = zeros(size(options_.tmp,1)+1);
tmp_(1:end-1,1:end-1) = options_.tmp;
tmp_ = tmp_ + [means;1] * [means;1]';
const = size(T,1) + 1;
projection_facility = zeros(NUM_MODELS,1);

for i = 1:NUM_MODELS
    ys_lists{i} = [const * ones(1,YF);options_.m(i).y_st_full];
    betamat_New{i} = ys_lists{i};
    ys_list_all{i} = ys_lists{i};
    ys_list_all{i} = ys_list_all{i}(:);
    tb{i} = (ys_list_all{i} ~=0) & (ys_list_all{i} ~= const);
    tc{i} = (ys_list_all{i} == const);
    ys_list_all{i} = nonzeros(ys_list_all{i});
    ys_vars{i} = ys_list_all{i}(ys_list_all{i} ~= const);
    beta = options_.m(i).beta;
%    df(i) = length(beta) + YF;
    df(i) = length( nonzeros(ys_lists{i}(:,leave)));
    Mask_const{i} = options_.m(i).Mask_const;
    X{i} = Mask_const{i};
    tX{i} = (Mask_const{i} ~= 0); % logical array showing location of ANY variable in the data matrix
    tXc{i} = zeros(size(Mask_const{i},1),size(Mask_const{i},2)); 
    for j = 1:size(Mask_const{i},2)
        tXc{i}(find(Mask_const{i}(:,j),1),j) = 1;
    end
    temp = (tXc{i} ~= 0); % logical array showing location of constants in the data matrix
    tXc{i} = temp;
    tX{i} = tX{i} & ~tXc{i}; % now tX shows location of non-constants in the data matrix
    X{i}(tX{i}) = means(ys_vars{i});
    X_noconst = options_.m(i).Mask;
    X_noconst(X_noconst ~= 0) = means(ys_vars{i});
    consts{i} = means(yf_list) - X_noconst' * beta;
    betamat_all{i} = zeros(YF*(ZZ+1),1);
    betamat_all{i}(tb{i}) = beta;
    betamat_all{i}(tc{i}) = consts{i}; 
    tt{i} = tb{i} | tc{i};
    betamat{i} = betamat_all{i}(find(tt{i}));
    bet{i} = zeros(length(betamat{i}),smpl);
    bet{i}(:,1) = betamat{i};
    betamat_bar{i} = betamat{i};
    betamat_bar_all{i} = betamat_all{i};
    ksi{i} = zeros(length(ys_list_all{i}),1);
    ksi_all{i} = zeros(YF*(ZZ+1),1);
    size_beta{i} = length(betamat{i}(:)); 
    %%%% If we make R_beta diagonal, all the matrices will remain
    %%%% block-diagonal, effectively generating equation-by-equation estimation
    %%%% -------------
    R_beta{i} = options_.m(i).R_beta;
    %R_beta = diag(diag(options_.R_beta));
    dR(i) = log(det(R_beta{i}));
    Sigma{i} = zeros(size(R_beta{i},1));
    inv_V = Mask_const{i} * (R_beta{i} \ Mask_const{i}');
    Q_bet = inv(tmp_(ys_list_all{i},ys_list_all{i}) .* inv_V);
    P_beta{i} = gain * Q_bet;
    P_b{i} = zeros(size(P_beta{i},1),size(P_beta{i},2),smpl);
    Q_beta{i} = sigma * Q_bet;
    T_beta{i} = ro * eye(length(betamat{i}+YF));
end
bet{NUM_MODELS+1} = zeros(ZZ+1,YF,smpl);
for i = 1:NUM_MODELS
    bet{NUM_MODELS+1}(:,:,1) = bet{NUM_MODELS+1}(:,:,1) + ...
        weights(i) * reshape(betamat_all{i},ZZ+1,YF);
end

PF_total = 0;

pj_threshold = options_.pj_threshold;
pj_mult=(1-pj_threshold)*2/pi;
pj_threshold_pi = pj_threshold - 0.15;
pj_mult_pi = (1-pj_threshold_pi)*2/pi;

nk = options_.nk;
mm   = size(T,2);
pp   = size(Y,1);
mu_in = (eye(mm) - T) * means;
mu = mu_in;
mus(:,1) = mu;
v       	= zeros(pp,smpl);
a         = zeros(mm,smpl+1);
% a(:,1)    = mu;
a(:,1) = means;
aK        = zeros(nk,mm,smpl+1);
aKM       = zeros(2,nk,mm,smpl+1);
filteratt = mu * ones(1,smpl);
QQ   = R*Q*transpose(R);
crit1       = 1.e-8;
n_y = length(ys_list);
beta_lim = 1e+8;
warning off

LG = zeros(size(dr_.jacobian,1));
LG(:,sort(dr_.order_var(dr_.nstatic+1:dr_.nstatic+dr_.npred))) = ...
  dr_.jacobian(:,1:dr_.npred);
EXP = zeros(mm);
vars = dr_.order_var(sort([ys_list; shock_list]));
fw = size(dr_.jacobian,2) - estim_params_.nvx;
fwd_sort = dr_.order_var(yf_list);
t_vars = dr_.npred+1:fw-dr_.nsfwrd;
A1 = dr_.jacobian(:,t_vars);
A2 = dr_.jacobian(:,fw-dr_.nsfwrd+1:fw);
A2_full = zeros(mm);
A2_full(:,fwd_sort) = A2;
B  = dr_.jacobian(:,fw+1:end);

Sum_A_SState = (LG + A1 + A2_full) * SState;

iF      	= zeros(pp,pp,smpl);
K       	= zeros(mm,pp,smpl);
L       	= zeros(mm,mm,smpl);
P       	= zeros(mm,mm,smpl+1);
P(:,:,1)    = Pstar;
crit    	= options_.kalman_tol;
rr      	= size(Q,1);
QQ      	= R*Q*transpose(R);
alphahat   	= zeros(mm,smpl);
etahat	   	= zeros(rr,smpl);
r 		   	= zeros(mm,smpl);
n_y = length(ys_list);
N = size(dr_.jacobian,1);

% block for forecasting models' predictions: highly specific! In our 5 
% forecasting models setup, only models 1 and 3 are pure ARs - AR(1) and
% AR(2). If you introduce your own sets, this will need to be changed
iM = [1;3];
tempF = zeros(2*length(yf_list),1);

Z = zeros(pp,mm);
for i=1:pp;
	Z(i,mf(i)) = 1;
end

t = 0;
notsteady = 1;

while notsteady && t<smpl
    t = t+1;
    v(:,t)      = Y(:,t) - a(mf,t);
    P(:,:,t)    = tril(P(:,:,t)) + transpose(tril(P(:,:,t),-1));
    if rcond(P(mf,mf,t)) < crit
    	return		
    end    
    iF(:,:,t)   = inv(P(mf,mf,t));
    K(:,:,t)    = T * P(:,mf,t) * iF(:,:,t);
    L(:,:,t)    = T - K(:,:,t)*Z;
    a(:,t+1)    = mu + T*a(:,t) + K(:,:,t)*v(:,t);    
    %-------------------------------------------------- 
    filteratt(:,t) = a(:,t) + P(:,mf,t)*iF(:,:,t) * v(:,t);
    %-------------------------------------------------- 
    % ALM forecasting block. This generates forecasts based on ALM, when
    % beliefs are already incorporated into the DSGE model, replacing
    % expectations, and the backward model is solved to produce the new
    % time-varying transmission mechanism. Then, this transmission
    % mechanism is simply iterated forward, ignoring any beliefs updating
    % that could have happened within the forecast window.
    
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk
        aK(jnk,:,t+jnk) = (mu + T * aK(jnk-1,:,t+jnk-1)')';
    end
%     % Individual models forecasting block - highly specific, check
%     % dynare_estimation.m!
%     % Model 1 - AR(1), Model 3 - AR(2)
%     for j = 1:length(iM)
%         ii = iM(j);
%         X{ii}(tX{ii}) = filteratt(ys_vars{ii},t);
%         X{ii}(tXc{ii}) = 1;
%         % These are the forecasted variables than we actually need in Model
%         % 1
%         aKM(j,1,:,t+1) = X{ii}' * betamat{ii};
%         if ii == 1
%             X{ii}(tX{ii}) = aKM(j,1,:,t+1);
%         elseif ii == 3
%             tempF(1:2:end) = aKM(j,1,:,t+1);
%             tempF(2:2:end) = filteratt(yf_list,t);
%             X{ii}(tX{ii}) = tempF;
%         end
%         for jnk=2:nk
%             aKM(j,jnk,:,t+jnk) = X{ii}' * betamat{ii};
%             if ii == 1
%                 X{ii}(tX{ii}) = aKM(j,jnk,:,t+jnk);
%             elseif ii == 3
%                 tempF(2:2:end) = tempF(1:2:end);
%                 tempF(1:2:end) = aKM(j,jnk,:,t+jnk);
%                 X{ii}(tX{ii}) = tempF;
%             end
%         end
%     end

    P(:,:,t+1)  = T*P(:,:,t)*transpose(T)-T*P(:,mf,t)*transpose(K(:,:,t)) + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
    %%% updating step: deriving new T and R; begins in period t >= upd_start
    if t >= upd_start
        % Updating ALL models
        for i = 1:NUM_MODELS
            if t == 1
                X{i}(tX{i}) = a(ys_vars{i},1);
            else
                X{i}(tX{i}) = filteratt(ys_vars{i},t-1);
            end
            X{i}(tXc{i}) = 1;
            betamat_all{i}(tt{i}) = betamat{i};
            err(:,i) = filteratt(yf_list,t) - X{i}' * betamat{i};
            ksi{i} = betamat{i} - betamat_bar{i};
            PX = P_beta{i} * X{i};
            F_beta = X{i}' * PX + R_beta{i};
            PX_iF = PX / F_beta;
            ksi{i} = T_beta{i} * (ksi{i} + PX_iF * err(:,i) );            
            ksi_all{i}(tt{i}) = ksi{i};
            betamat_New{i} = reshape(betamat_bar_all{i} + ksi_all{i},(ZZ+1),YF);
            P_b{i}(:,:,t) = P_beta{i};
            P_beta{i} = P_beta{i} - PX_iF * PX';
            P_beta{i} = T_beta{i} * P_beta{i} * T_beta{i}' + Q_beta{i};
            Sigma{i} = Sigma{i} + err(:,i) * err(:,i)';
            tS(:,:,i) = Sigma{i}(leave,leave);
            if any(isnan(betamat_New{i}(:))) || norm(betamat_New{i}) > beta_lim % beliefs often explode under GSG
                return;
            end
            % The only thing which remains to be settled is whether this
            % model retains the new betamat!
            betama_ = [betamat_New{i}(2:n_y+1,:)' betamat_New{i}(n_y+2:end,:)' * T(shock_list,shock_list)];
            betama_ = betama_(:,reorder);
            EXP(:,vars) = A2 * betama_;
            T_New = - (A1 + EXP) \ LG;
            if any(isnan(T_New(:))) 
                return;
            end
            hit = any(abs(eig(T_New)) > pj_threshold);
            if hit > 0
                projection_facility(i) = projection_facility(i) +  1;
                 betamat_New{i} = reshape(betamat_all{i},(ZZ+1),YF);
            else
                betamat{i} = betamat_bar{i} + ksi{i};
            end
            
            if save_print % printing information about individual forecasting model performance
                fprintf('%8.2f %7.2f ',t * log(det(Sigma{i}/t)),t*log(det(squeeze(tS(:,:,i)/t))));
                fprintf('%5.2f %2i %2i',weights(i),projection_facility(i),i);
            end
            bet{i}(:,t) = betamat{i};
        end
        % Integrating predictions of all the models. 
        if (t >= start) && (t >= length(leave)+1)
            for i = 1:NUM_MODELS
%                 AIC(i) = t * log(det(squeeze(tS(:,:,i)/t))) + 2 * df(i);
               BIC(i) = t * log(det(squeeze(tS(:,:,i)/t))) + log(t) * df(i);
            end
        end
        exp_dBIC = exp(-0.5 * (BIC - min(BIC)));
        sum_exp_dBIC = sum(exp_dBIC);
        weights_New = exp_dBIC / sum_exp_dBIC; 
        if save_print % printing information about aggregate forecasting model performance
            Sig_joint = Sig_joint + (err*weights) * (err*weights)';
            fprintf(' %8.2f %7.2f %3i\n',t * log(det(Sig_joint/t)),t * log(det(Sig_joint(leave,leave)/t)),t);
        end
        errors(:,1:NUM_MODELS,t) = err;
        errors(:,NUM_MODELS+1,t) = err * weights;
        exps(:,t) = filteratt(yf_list,t) - errors(:,NUM_MODELS+1,t);
        weight(:,t) = weights_New;
%         % If equal weights are assumed, the next line should be commented
%          weights = weights_New;
        betamat_Int = weights(1) * betamat_New{1};
        for i = 2:NUM_MODELS
            betamat_Int = betamat_Int + weights(i) * betamat_New{i};
        end
        bet{NUM_MODELS+1}(:,:,t) = betamat_Int;
        if any(isnan(betamat_Int(:))) || norm(betamat_Int) > beta_lim % beliefs often explode under GSG
            return;
        end
        betama_ = [betamat_Int(2:n_y+1,:)' betamat_Int(n_y+2:end,:)' * T(shock_list,shock_list)];
        betama_ = betama_(:,reorder);
        EXP(:,vars) = A2 * betama_;
        T_New = - (A1 + EXP) \ LG;
        if any(isnan(T_New(:))) 
            return;
        end
        [V_,D_] = eig(T_New);
        DD_ = diag(D_);
        hit = any(abs(DD_) > pj_threshold);
        if options_.corr == 0 || hit == 0
            R = - (A1 + EXP) \ B;
            T = T_New(dr_.order_var,dr_.order_var);
            R = R(dr_.order_var,:);
            QQ = R * Q * R';
        elseif hit > 0
            PF_total = PF_total + 1;
            ii = find(DD_>pj_threshold);
            DD_(ii) = pj_threshold + pj_mult * atan((DD_(ii)-pj_threshold)/pj_mult);
            % Eigenvalue adjustment can lead to tiny imaginary components
            % in T_New, which results in imaginary likelihood and the point
            % being discarded
            T_New = real(V_ * diag(DD_) * pinv(V_));
            T = T_New(dr_.order_var,dr_.order_var);
            R = - (A1 + EXP) \ B;
            R = R(dr_.order_var,:);
            QQ = R * Q * R';
        end
        mu_new = (A1 + EXP)\ (Sum_A_SState - A2 * betamat_Int(1,:)');
        mu = mu_new(dr_.order_var);
    end
    mus(:,t+1) = mu;
    TT(:,:,t+1) = T;
    RR(:,:,t+1) = R;    
end

t = smpl+1;
while t>2
	t = t-1;
    r(:,t-1) = Z' * iF(:,:,t) * v(:,t) + L(:,:,t)'*r(:,t);
    alphahat(:,t)	= a(:,t) + P(:,:,t)*r(:,t-1);
	etahat(:,t)		= Q * RR(:,:,t)'*r(:,t);
end

% GM_simulation(TT,RR,Q,mus,lgy_,dr_.order_var)
% This calculates average variances of some variables which we used to see
% if the Great Moderation was reproduced, and variance decomposition. All
% time-varying. You will need to comment out this line if using any other
% model!
% The simulation is performed at the posterior mode beliefs only, while the
% paper typically reports values obtained using MCMC draws.
GM_simulation_var_dec(TT,RR,Q,mus,lgy_,dr_.order_var)

% Saving time-varying transmission mechanism, beliefs, expectations implied
% by the average model (average PLM expectations), and expectation errors
oo_.TT = TT;
oo_.RR = RR;
oo_.mus = mus;
oo_.Q = Q;
oo_.bet = bet;
oo_.exps = exps;
oo_.errors = errors;
oo_.P_b = P_b;

% The following block prints the number of projection facility hits for
% individual model and the aggregate model, and also plot beliefs graphs.
% Time units are highly specific, will have to be adjusted if any other
% sample is used

if save_print
    if any(projection_facility)
        fprintf('\n %4i violations \n',projection_facility);
        fprintf('\n %4i violations in aggregate\n',PF_total);    
    end
    t_in = 1965.;
    t_b = (t_in:0.25:2008.75)';
    t_b(smpl+1:end) = [];

    for m = 1:NUM_MODELS
        offset = 0;
        for i = 1:length(yf_list)
            vars_X = nonzeros(ys_lists{m}(2:end,i));
            leg = cell(length(vars_X)+1,1);
            leg{1} = 'const';
            for j = 2:length(leg)
                leg{j} = deblank(lgy_(dr_.order_var(vars_X(j-1)),:));
            end
            figure(i+20*m)
            plot(datenum(t_b,1,1),bet{m}(offset+1:offset+length(leg),:)');
            datetick('x',17,'keeplimits');
            offset = offset+length(leg);
            title(['Model ',int2str(m),', equation for ',lgy_(dr_.order_var(yf_list(i)),:)]);
            legend('v6',leg)
        end
    end

    % plotting aggregate beliefs
    leg = cell(length(ys_list)+1,1);
    leg{1} = 'const';
    for j = 2:length(leg)
        leg{j} = deblank(lgy_(dr_.order_var(ys_list(j-1)),:));
    end

    for j = 1:length(yf_list)
        figure(i+20*(m + 1)+j)
        plot(datenum(t_b,1,1),squeeze(bet{NUM_MODELS+1}(:,j,:))');
        datetick('x',17,'keeplimits');
        offset = offset+length(leg);
        title(['Aggregate Model, equation for ',lgy_(dr_.order_var(yf_list(j)),:)]);
        legend('v6',leg)
    end

    % plotting model weights
    leg_m = cell(NUM_MODELS,1);
    for j = 1:NUM_MODELS
        leg_m{j} = ['Model ',int2str(j)];
    end
    figure(i+20*m +1)
    plot(weight');
    title('Models weights');
    legend('v6',leg_m)
    
    leg_err = [leg_m;'Average'];    
    for i = 1:length(yf_list)
        figure(20*(NUM_MODELS+1) + i);
        plot(squeeze(errors(i,:,:))');
        title([' Errors for ',lgy_(dr_.order_var(yf_list(i)),:)]);
        legend('v6',leg_err);
    end
    
end

% Graphs of perceived volatility of the beliefs for inflation; assumes that
% Model 1 is AR2
var_pinf_const = squeeze(P_b{1}(10,10,:));
var_pinf_pers = squeeze(P_b{1}(11,11,:) + P_b{1}(12,12,:) + 2 * P_b{1}(11,12,:));
figure(100)
plot(datenum(t_b,1,1),var_pinf_const);
datetick('x',17,'keeplimits');
title('Perceived volatility of beliefs about constant: AR2 pinf model')

figure(101)
plot(datenum(t_b,1,1),var_pinf_pers);
datetick('x',17,'keeplimits');
title('Perceived volatility of beliefs about persistence: AR2 pinf model')

