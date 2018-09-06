function LIK = DiffuseLikelihoodKMAL1_options_const(T,R,Q,Pinf,Pstar,Y,SState,start,endpars)

% stephane.adjemian@cepremap.cnrs.fr [07-19-2004]
%
% Same as DiffuseLikelihoodH1 without measurement error.

%  The program was modified by Sergey Slobodyan and Raf Wouters to allow for
%  adaptive learning
%  Sergey.Slobodyan@cerge-ei.cz, August 18, 2011

% This program delivers LIK for a model with Kalman filter learning
% WITH a constant

global bayestopt_ options_  dr_ estim_params_ lgy_
crit        = options_.kalman_tol;

% In this version, ro is estimated. gain and sigma set as described in the
% paper. Commented out line fixes ro at 1 and estimated gain and sigma
% gain = 0.030969052313799; ro = endpars; sigma = 0.003;
% gain = endpars(1); ro = 1; sigma = endpars(2);

% Beliefs updating is switched off
gain = 0.030969052313799; ro = 0; sigma = 0.003;

ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;
reorder = options_.reorder;  
upd_start = options_.upd_start;
lik_penalty = options_.lik_penalty;
YF = length(yf_list);
Z = length(shock_list) + length(ys_list);
NUM_MODELS = options_.num_mod;
df = zeros(NUM_MODELS,1);
dR = zeros(NUM_MODELS,1);
BIC = ones(NUM_MODELS,1);
AIC = ones(NUM_MODELS,1);
exp_dBIC = BIC;
exp_dAIC = BIC;
smpl = size(Y,2);
weights = ones(NUM_MODELS,1) / NUM_MODELS;
err = zeros(length(yf_list),NUM_MODELS);
leave = [1:4,7];
Sig_joint = zeros(length(yf_list));
tS = zeros(length(leave),length(leave),NUM_MODELS);
means = SState(dr_.order_var);
tmp_ = zeros(size(options_.tmp,1)+1);
tmp_(1:end-1,1:end-1) = options_.tmp;
tmp_ = tmp_ + [means;1] * [means;1]';
const = size(T,1) + 1;
projection_facility = zeros(NUM_MODELS,1);

% Forecasting models could have different numbers of variables in
% individual equations, even within a single model. That is why all the
% relevant arrays are set as cell arrays
% Basic logic here: the beliefs are treated as a column vector for the
% purposes of beliefs Kalman filter. When the beliefs are plugged into the
% model to produce ALM, they are transformed into a rectangular array as in
% constant gain case. Mask variables are used to transform the vector into
% the array.
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
    betamat_all{i} = zeros(YF*(Z+1),1);
    betamat_all{i}(tb{i}) = beta;
    betamat_all{i}(tc{i}) = consts{i}; 
    tt{i} = tb{i} | tc{i};
    betamat{i} = betamat_all{i}(find(tt{i}));
    betamat_bar{i} = betamat{i};
    betamat_bar_all{i} = betamat_all{i};
    ksi{i} = zeros(length(ys_list_all{i}),1);
    ksi_all{i} = zeros(YF*(Z+1),1);
    size_beta{i} = length(betamat{i}(:)); 
    %%%% R_beta is the matrix \Sigma from the paper
    %%%% If we make R_beta diagonal, all the matrices will remain
    %%%% block-diagonal, effectively generating equation-by-equation estimation
    %%%% -------------
    R_beta{i} = options_.m(i).R_beta;
    %R_beta = diag(diag(options_.R_beta));
    Sigma{i} = zeros(size(R_beta{i},1));
    dR(i) = log(det(R_beta{i}));
    %%%% Q_bet is the matrix inv (X^T * inv(\Sigma) * X) from the paper
    inv_V = Mask_const{i} * (R_beta{i} \ Mask_const{i}');
    Q_bet = inv(tmp_(ys_list_all{i},ys_list_all{i}) .* inv_V);
    % P_beta is P_{1|0} from the paper
    P_beta{i} = gain * Q_bet;
    % Q_beta is V from the paper
    Q_beta{i} = sigma * Q_bet;
    T_beta{i} = ro * eye(length(betamat{i}+YF));
end

PF_total = 0;

% The next 4 variables could be used in alternative projection facility
% which smoothes out the standard PF's 1-0 behavior. Logic of the
% smoothing: the eigenvalues of the T matrix between pj_threshold and 
% infinity are mapped into [pj_threshold;1] interval, and the eigenvectors 
% are then used to reconstruct T in a manner similar to MH manipulation of 
% the Hessian 
pj_threshold = options_.pj_threshold;
pj_mult=(1-pj_threshold)*2/pi;
pj_threshold_pi = pj_threshold - 0.15;
pj_mult_pi = (1-pj_threshold_pi)*2/pi;

mf = bayestopt_.mf;
lik_beta = zeros(NUM_MODELS,smpl);
mm   = size(T,2);
pp   = size(Y,1);
mu_in = (eye(mm) - T) * means;
mu = mu_in;
a = mu;
f_a    = mu;
f_a_old    = mu;
v    = zeros(pp,1);
v_s    = zeros(pp,smpl);
dF = 1;
QQ   = R*Q*transpose(R);
lik  = zeros(smpl+1,1);
LIK  = Inf;
lik(smpl+1) = smpl*pp*log(2*pi);
crit1       = 1.e-8;
n_y = length(ys_list);
beta_lim = 1e+8;
LIK_lim = 1e+8;
warning off

% Matrix A_0 from the jacobia_ variable.
LG = zeros(size(dr_.jacobian,1));
% Need to re-sort the columns (originally in alpha order) into order_var.
LG(:,sort(dr_.order_var(dr_.nstatic+1:dr_.nstatic+dr_.npred))) = ...
  dr_.jacobian(:,1:dr_.npred);
% Preparations to construct matrix A_2 * betamat
EXP = zeros(mm);
% A list of ALL variables present on the RHS of forecasting equations, in
% order_var
vars = dr_.order_var(sort([ys_list; shock_list]));
fwd_sort = dr_.order_var(yf_list);
fw = size(dr_.jacobian,2) - estim_params_.nvx;
t_vars = dr_.npred+1:fw-dr_.nsfwrd;
% A_1 matrix of the jacobia_ variable. Used only in its original order
% (columns alpha)
A1 = dr_.jacobian(:,t_vars);
% matrix A_2 from jacobian, used only in its original (alpha) form.
A2 = dr_.jacobian(:,fw-dr_.nsfwrd+1:fw);
% matrix A_2_full is useful in updating the constant
A2_full = zeros(mm);
A2_full(:,fwd_sort) = A2;
% Matrix which multiplies innovations in jacobia_
B  = dr_.jacobian(:,fw+1:end);

% Vector that is used in updating the mu in VAR representation of the
% model, given beliefs
Sum_A_SState = (LG + A1 + A2_full) * SState;

t = 0;

P = Pstar;
notsteady = 1;

while notsteady && t<smpl
    t = t+1;
%    v = Y(:,t) - a(mf) - trend(:,t);
    % We are working with the original variables, not detrended ones in
    % this subroutine!
    v = Y(:,t) - a(mf);
    v_s(:,t) = v;
    P=tril(P)+tril(P,-1)';
    oldP  = P;
    F  = P(mf,mf);
    dF = det(F);
    if rcond(F) < crit
        if ~all(abs(F)<crit)
            return
        else
            f_a_old = f_a;
            f_a = a;
            a  = T*a;
            P     = T*P*transpose(T)+QQ;
        end
    else
        iF   = inv(F);
        K    = T*P(:,mf) * iF;
        KwithoutT  = P(:,mf) * iF;
        f_a_old = f_a;
        f_a  = a + KwithoutT*v;
        a          = mu + T*a + K*v;    
        lik(t) = log(dF) + v' * iF * v;
        P  = T*P*T' - T*P(:,mf)*K' + QQ;
    end    
    %%% updating step: deriving new T and R; begins in period t >= upd_start
    if t >= upd_start
        % Updating ALL models
        for i = 1:NUM_MODELS
            % Forming the X SURE-like matrix 
            X{i}(tX{i}) = f_a_old(ys_vars{i});
            X{i}(tXc{i}) = 1;
            % generating correct vector of beliefs
            betamat_all{i}(tt{i}) = betamat{i};
            % prediction error
            err(:,i) = f_a(yf_list) - X{i}' * betamat{i};
            % Kalman filter step for beliefs
            ksi{i} = betamat{i} - betamat_bar{i};
            PX = P_beta{i} * X{i};
            F_beta = X{i}' * PX + R_beta{i};
            PX_iF = PX / F_beta;
            ksi{i} = T_beta{i} * (ksi{i} + PX_iF * err(:,i) );            
            ksi_all{i}(tt{i}) = ksi{i};
            % transforming the updated belief vector into rectangular form
            betamat_New{i} = reshape(betamat_bar_all{i} + ksi_all{i},(Z+1),YF);
            % updating one step ahead var-covar matrix for beliefs
            P_beta{i} = P_beta{i} - PX_iF * PX';
            P_beta{i} = T_beta{i} * P_beta{i} * T_beta{i}' + Q_beta{i};
            % empirical var-covar matrix of forecast errors
            Sigma{i} = Sigma{i} + err(:,i) * err(:,i)';
            % we use only observable variables here
            tS(:,:,i) = Sigma{i}(leave,leave);
            if any(isnan(betamat_New{i}(:))) || norm(betamat_New{i}) > beta_lim % beliefs often explode under GSG
                LIK = LIK_lim;
                return;
            end
            % The only thing which remains to be settled is whether this
            % model retains the new betamat!
            % Per standard convention in AL literature, stochastic processes of
            % shocks are known to the agents exactly, therefore expectation of
            % future shocks are a simple matrix multiplication
            % NOTE that all 'shocks' are assumed to be members of the same DYNARE
            % group of purely predetermined variables. It doesn't make sense to
            % form expectations of y_{t+1} knowing w_t if w_t is i.i.d.,
            % anyway.
            % So far, shocks_list is always empty!
            betama_ = [betamat_New{i}(2:n_y+1,:)' betamat_New{i}(n_y+2:end,:)' * T(shock_list,shock_list)];
            % Columns of the beliefs matrix, with shocks columns translated one
            % period ahead, are put into purely alpha order and transposed to 
            % correspond to the was coluns of the matrix A_2 of jacobia_ are
            % ordered.
            betama_ = betama_(:,reorder);
            % Matrix EXP has the same dimensions as T, but only the columns
            % corresponding to the variables that are used on the RHS of
            % forecasting equation will be updated. This formulation allows ANY
            % model variable, including static or purely forward-looking, to be
            % used in forecasting equations
            EXP(:,vars) = A2 * betama_;
            % Solving the resulting purely backward-looking model
            T_New = - (A1 + EXP) \ LG;
            if any(isnan(T_New(:))) 
                LIK = LIK_lim;
                return;
            end
            % Are any eigenvalues of T_New unstable?
            hit = any(abs(eig(T_New)) > pj_threshold);
            %           OLD PF, 0-1
            % this is PF for individual model: whether a transmission
            % mechanism, formed if this were the only forecasting model the
            % agents would use, would be stable.
            if hit > 0 
                projection_facility(i) = projection_facility(i) +  1;
                 betamat_New{i} = reshape(betamat_all{i},(Z+1),YF);
            else
                betamat{i} = betamat_bar{i} + ksi{i};
            end
        end
        % INTEGRATING PREDICTIONS OF ALL MODELS
        % Forming BIC (or AIC) statistics
        if (t >= start) && (t >= length(leave)+1)
            for i = 1:NUM_MODELS
%                 AIC(i) = t * log(det(squeeze(tS(:,:,i)/t))) + 2 * df(i);
               BIC(i) = t * log(det(squeeze(tS(:,:,i)/t))) + log(t) * df(i);
            end
        end
        % Calculating new model weights based on past performance
        exp_dBIC = exp(-0.5 * (BIC - min(BIC)));
        sum_exp_dBIC = sum(exp_dBIC);
        weights_New = exp_dBIC / sum_exp_dBIC; 
%       If equal weights are assumed, the next line should be commented.
%          weights = weights_New;

        % Generating aggregate beliefs, using new model weights
        betamat_Int = weights(1) * betamat_New{1};
        for i = 2:NUM_MODELS
            betamat_Int = betamat_Int + weights(i) * betamat_New{i};
        end
        if any(isnan(betamat_Int(:))) || norm(betamat_Int) > beta_lim % beliefs often explode under GSG
            LIK = LIK_lim;
            return;
        end
        % Constructing new transmission mechanism - this time with
        % aggregate model
        betama_ = [betamat_Int(2:n_y+1,:)' betamat_Int(n_y+2:end,:)' * T(shock_list,shock_list)];
        betama_ = betama_(:,reorder);
        EXP(:,vars) = A2 * betama_;
        T_New = - (A1 + EXP) \ LG;
        if any(isnan(T_New(:))) 
            LIK = LIK_lim;
            return;
        end
        % Projection facility on aggregate model is done using smoothed PF
        % mechanism. Typically, after individual models' PFs are invoked,
        % the aggregate model produces a stable transmission mechanism
        [V_,D_] = eig(T_New);
        DD_ = diag(D_);
        hit = any(abs(DD_) > pj_threshold);
        if options_.corr == 0 || hit == 0
            R = - (A1 + EXP) \ B;
            T = T_New(dr_.order_var,dr_.order_var);
            R = R(dr_.order_var,:);
            QQ = R * Q * R';
        elseif hit > 0
            lik(t) = lik(t)+lik_penalty;
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
        % adjusting the constant in the transmission mechanism
        mu_new = (A1 + EXP)\ (Sum_A_SState - A2 * betamat_Int(1,:)');
        % putting the constant into proper order
        mu = mu_new(dr_.order_var);
    end
    notsteady   = ~(max(max(abs(P - oldP)))<crit);
end

reste = smpl-t;
lik(t) = lik(t) + reste*log(dF);

LIK    = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);% Minus the log-likelihood.