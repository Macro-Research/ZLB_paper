function LIK = DiffuseLikelihoodL1_corr_(T,R,Q,Pinf,Pstar,Y,trend,start,gain);

% stephane.adjemian@cepremap.cnrs.fr [07-19-2004]
%
% Same as DiffuseLikelihoodH1 without measurement error.

%  The program was modified by Sergey Slobodyan and Raf Wouters to allow for
%  adaptive learning
%  Sergey.Slobodyan@cerge-ei.cz, June 27, 2010

% This program delivers LIK for a model with simple constant gain learning
% WITHOUT constant
  global bayestopt_ options_ jacobia_ dr_ estim_params_ lgy_

betamat = options_.betamat;
SecMom = options_.SecondMoments;
SecMom_New = SecMom;
ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;
reorder = options_.reorder;  

upd_start = options_.upd_start;
lik_penalty = options_.lik_penalty;
projection_facility = 0;
ridge = options_.ridge;
ridge_adj = options_.ridge_adj;
smsize = size(SecMom,1);
mf = bayestopt_.mf;
smpl = size(Y,2);
mm   = size(T,2);
pp   = size(Y,1);
a    = zeros(mm,1);
f_a    = zeros(mm,1);
f_a_old    = zeros(mm,1);
v    = zeros(pp,1);
v_s    = zeros(pp,smpl);
dF = 1;
QQ   = R*Q*transpose(R);
t    = 0;
lik  = zeros(smpl+1,1);
LIK  = Inf;
lik(smpl+1) = smpl*pp*log(2*pi);
notsteady   = 1;
crit        = options_.kalman_tol;
crit1       = 1.e-8;
reste       = 0;
n_y = length(ys_list);
warning off

% Matrix A_0 from the jacobia_ variable.
LG = zeros(size(dr_.jacobian,1));
% Need to re-sort the columns (originally in alpha order) into order_var.
LG(:,sort(dr_.order_var(dr_.nstatic+1:dr_.nstatic+dr_.npred))) = ...
  dr_.jacobian(:,1:dr_.npred);
% Preparations to construct matrix A_2 * betamat
EXP = zeros(size(dr_.jacobian,1));
% A list of ALL variables present on the RHS of forecasting equations, in
% order_var
vars = dr_.order_var(sort([ys_list; shock_list]));
fw = size(dr_.jacobian,2) - estim_params_.nvx;
t_vars = dr_.npred+1:fw-dr_.nsfwrd;
% A_1 matrix of the jacobia_ variable. Used only in its original order
% (columns alpha)
A1 = dr_.jacobian(:,t_vars);
% matrix A_2 from jacobian, used only in its original (alpha) form.
A2 = dr_.jacobian(:,fw-dr_.nsfwrd+1:fw);
% Matrix which multiplies innovations in jacobia_
B  = dr_.jacobian(:,fw+1:end);

RO = T(shock_list,shock_list);

% Diffuse priors block DOES NOT work in AL, it shouldn't be used in this
% subroutine
t = 0;
while rank(Pinf,crit) && t<smpl
    t = t+1;
    v 		 	= Y(:,t) - a(mf) - trend(:,t);
    v_s(:,t) = v;
    Finf = Pinf(mf,mf);
    if rcond(Finf) < crit
        if ~all(abs(Finf(:))<crit)
            return		
        else
            iFstar	= inv(Pstar(mf,mf));
            dFstar	= det(Pstar(mf,mf));
            Kstar	= Pstar(:,mf)*iFstar;
            lik(t)	= log(dFstar) + transpose(v)*iFstar*v;
            Pinf	= T*Pinf*transpose(T);
            Pstar	= T*(Pstar-Pstar(:,mf)*transpose(Kstar))*transpose(T)+QQ;
            f_a  = a + Kstar*v;
            a		= T*(a+Kstar*v);
        end
    else
        lik(t)	= log(det(Finf));
        iFinf 	= inv(Finf);
        Kinf	= Pinf(:,mf)*iFinf;					%%	premultiplication by the transition matrix T is removed (stephane) 
        Fstar	= Pstar(mf,mf);
        Kstar	= (Pstar(:,mf)-Kinf*Fstar)*iFinf; 	%%	premultiplication by the transition matrix T is removed (stephane)
        Pstar	= T*(Pstar-Pstar(:,mf)*Kinf'-Pinf(:,mf)*Kstar')*T'+QQ;
        Pinf	= T*(Pinf-Pinf(:,mf)*Kinf')*T';
        f_a  = a + Kinf*v;
        a		= T*(a+Kinf*v);					
    end
end

if t == smpl
    error(['There isn''t enough information to estimate the initial' ... 
	   ' conditions of the nonstationary variables']);                   
end                                                                    

P = Pstar;
notsteady = 1;
F_singular = 1;

while notsteady && t<smpl
    t = t+1;
    v = Y(:,t) - a(mf) - trend(:,t);
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
            if options_.kalman_algo > 200
                f_a(1) = 1;
            end
            a  = T*a;
            P     = T*P*transpose(T)+QQ;
        end
    else
        F_singular = 0;
        iF   = inv(F);
        K    = T*P(:,mf)*iF;
        KwithoutT  = P(:,mf)*iF;
        f_a_old = f_a;
        f_a  = a + KwithoutT*v;
        if options_.kalman_algo > 200
            f_a(1) = 1;
        end
        a          = T*a + K*v;    
        lik(t) = log(dF) + transpose(v)*iF*v;
        P  = T*P*T' - T*P(:,mf)*K' + QQ;
    end    
    %%% updating step: deriving new T and R; begins in period t >= upd_start
    if t >= upd_start
        % RHS variables used in forecasting equations, endogenous first,
        % shocks last, both ordered alpha withing groups
        z = [f_a_old(ys_list); f_a(shock_list)];
        % Forecast errors
        err = f_a(yf_list) - betamat' * z;
        % Updating of the 2nd moments matrix
        SecMom_New = SecMom + gain * (z * z' - SecMom);
        if min(abs(eig(SecMom_New))) < ridge
            % beliefs updating with Ridge correction factor
            betamat_New = betamat + gain * ((SecMom_New + ridge / ridge_adj * eye(smsize)) \ (z * err'));
        else
            betamat_New = betamat + gain * (SecMom_New \ (z * err'));
        end
        % Per standard convention in AL literature, stochastic processes of
        % shocks are known to the agents exactly, therefore expectation of
        % future shocks are a simple matrix multiplication
        % NOTE that all 'shocks' are assumed to be members of the same DYNARE
        % group of purely predetermined variables. It doesn't make sense to
        % form expectations of y_{t+1} knowing w_t if w_t is i.i.d.,
        % anyway.
        betama_ = [betamat_New(1:n_y,:)' betamat_New(n_y+1:end,:)' * RO];
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
        % Are any eigenvalues of T_New unstable?
        hit = any(abs(eig(T_New)) > options_.qz_criterium);
        % No PF is desired or no need to invoke it
        if options_.corr == 0 || hit == 0
            R = - (A1 + EXP) \ B;
            % T and R thus obtained are ordered alpha, need to be put into
            % order_var to be used in the main Kalman filter step
            T = T_New(dr_.order_var,dr_.order_var);
            R = R(dr_.order_var,:);
            QQ = R * Q * R';
            % Finally, we can accept the update in beliefs
            SecMom = SecMom_New;
            betamat = betamat_New;
        end
        % Projection facility is needed and allowed. No updating happens,
        % just PF counter and likelihood penalty are incremented
        if hit > 0
            projection_facility = projection_facility +  1;
            lik(t) = lik(t)+lik_penalty;
        end
    end
    notsteady   = ~(max(max(abs(P - oldP)))<crit);
end

if F_singular == 1
error(['The variance of the forecast error remains singular until the' ...
  'end of the sample'])
end

% This piece of the code is not expected to run, as under constant gain
% learning, matrix P should be time-varying anyway
reste = smpl-t;
while t<smpl
    t=t+1;
    v = Y(:,t) - a(mf) - trend(:,t);
    v_s(:,t) = v;
    a = T*a + K*v;
end
lik(t) = lik(t) + reste*log(dF);

LIK    = .5*(sum(lik(start:end))-(start-1)*lik(smpl+1)/smpl);% Minus the log-likelihood.