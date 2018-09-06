function [alphahat,etahat,a, aK, filteratt,bet] = DiffuseKalmanSmootherL1_options_const(T,R,Q,Pinf1,Pstar1,Y,SState,pp,mm,smpl,mf,gain)
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
%  Sergey.Slobodyan@cerge-ei.cz, June 27, 2010

% This program calculates smmothed values and saves time-varying AL objects 
% for a model with simple constant gain learning WITH a constant

global options_ dr_ estim_params_ lgy_ oo_

betamat = options_.betamat;
SecMom = options_.SecondMoments;
if options_.kalman_algo > 200 % betamat is just the slopes, need to add 
    % proper intercepts to the beliefs. SecMom is constructed as a
    % variance-covariance matrix, if there's a constant we have to turn it
    % into 2nd moments by adding outer product of means vector
    means = SState(dr_.order_var);
    joint_list = [options_.ys_list; options_.shock_list];
    const = means(options_.yf_list) - betamat' * means(joint_list);
    betamat = [const'; betamat];
    SM = zeros(size(SecMom,1)+1);
    SM(2:end,2:end) = SecMom + means(joint_list) * means(joint_list)';
    SM(1,1) = 1;
    SM(2:end,1) = means(joint_list);
    SM(1,2:end) = means(joint_list)';  
    SecMom = SM;
end

SecMom_New = SecMom;
ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;
reorder = options_.reorder;  

upd_start = options_.upd_start;
ridge = options_.ridge;
ridge_adj = options_.ridge_adj;
smsize = size(SecMom,1);
projection_facility = 0;

nk = options_.nk;
spinf   	= size(Pinf1);
spstar  	= size(Pstar1);
v       	= zeros(pp,smpl);
N = length(SState);
mu_in      = (eye(N) - T) * SState(dr_.order_var);
mu = mu_in;
Sstates      = zeros(mm,smpl+1); 
mus(:,1)     = mu_in;
a         = mu * ones(1,smpl+1);
aK          = zeros(nk,mm,smpl+1);
filteratt = mu * ones(1,smpl+1);
iF      	= zeros(pp,pp,smpl);
Fstar   	= zeros(pp,pp,smpl);
iFinf   	= zeros(pp,pp,smpl);
K       	= zeros(mm,pp,smpl);
L       	= zeros(mm,mm,smpl);
Linf    	= zeros(mm,mm,smpl);
Kstar   	= zeros(mm,pp,smpl);
P       	= zeros(mm,mm,smpl+1);
Pstar   	= zeros(spstar(1),spstar(2),smpl+1); Pstar(:,:,1) = Pstar1;
Pinf    	= zeros(spinf(1),spinf(2),smpl+1); Pinf(:,:,1) = Pinf1;
crit    	= options_.kalman_tol;
crit1       = 1.e-8;
steady  	= smpl;
rr      	= size(Q,1);
QQ      	= R*Q*transpose(R);
alphahat   	= zeros(mm,smpl);
etahat	   	= zeros(rr,smpl);
r 		   	= zeros(mm,smpl);
n_y = length(ys_list);

bet = zeros(length(ys_list)+length(shock_list)+1,length(yf_list),smpl);
bet(:,:,1) = betamat;

% Matrix A_0 from the jacobia_ variable.
LG = zeros(size(dr_.jacobian,1));
% Need to re-sort the columns (originally in alpha order) into order_var.
LG(:,sort(dr_.order_var(dr_.nstatic+1:dr_.nstatic+dr_.npred))) = ...
  dr_.jacobian(:,1:dr_.npred);
% Preparations to construct matrix A_2 * betamat
EXP = zeros(size(dr_.jacobian,1));
% A list of ALL variables present on the RHS of forecasting equations, in
% order_var
fwd_sort = sort(dr_.order_var(mm-dr_.nsfwrd+1:end)); 
% the same as dr_.order_var(yf_list), if list of forwards is correct
fwd_sort = dr_.order_var(yf_list);
vars = dr_.order_var(sort([ys_list; shock_list]));
fw = size(dr_.jacobian,2) - estim_params_.nvx;
t_vars = dr_.npred+1:fw-dr_.nsfwrd;
% A_1 matrix of the jacobia_ variable. Used only in its original order
% (columns alpha)
A1 = dr_.jacobian(:,t_vars);
% matrix A_2 from jacobian, used only in its original (alpha) form.
A2 = dr_.jacobian(:,fw-dr_.nsfwrd+1:fw);
% matrix A_2_full is useful in updating the constant
A2_full = zeros(N);
A2_full(:,fwd_sort) = A2;
% Matrix which multiplies innovations in jacobia_
B  = dr_.jacobian(:,fw+1:end);

RO = T(shock_list,shock_list);

Z = zeros(pp,mm);
for i=1:pp;
	Z(i,mf(i)) = 1;
end

t = 0;
TT(:,:,1) = T;
RR(:,:,1) = R;

d = t;
P(:,:,d+1) = Pstar(:,:,d+1);
iFinf = iFinf(:,:,1:d);
Linf  = Linf(:,:,1:d);
Fstar = Fstar(:,:,1:d);
Kstar = Kstar(:,:,1:d);
Pstar = Pstar(:,:,1:d);
Pinf  = Pinf(:,:,1:d);
notsteady = 1;
while notsteady && t<smpl
    t = t+1;
    v(:,t)      = Y(:,t) - a(mf,t);
    P(:,:,t)=tril(P(:,:,t))+transpose(tril(P(:,:,t),-1));
    if rcond(P(mf,mf,t)) < crit
    	return		
    end    
    iF(:,:,t)   = inv(P(mf,mf,t));
    K(:,:,t)    = T*P(:,mf,t)*iF(:,:,t);
    L(:,:,t)    = T-K(:,:,t)*Z;
    a(:,t+1)    = mu + T*a(:,t) + K(:,:,t)*v(:,t);    
    %-------------------------------------------------- 
    KwithoutT(:,:,t)	 	= P(:,mf,t)*iF(:,:,t);
    filteratt(:,t)=a(:,t) + KwithoutT(:,:,t)*v(:,t);
    %-------------------------------------------------- 
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk
        aK(jnk,:,t+jnk) = (mu + T * aK(jnk-1,:,t+jnk-1)')';
    end
    P(:,:,t+1)  = T*P(:,:,t)*transpose(T)-T*P(:,mf,t)*transpose(K(:,:,t)) + QQ;
    notsteady   = ~(max(max(abs(P(:,:,t+1)-P(:,:,t))))<crit);
    if t >= upd_start
        if t > 1
            z = [filteratt(ys_list,t-1); filteratt(shock_list,t)];
        else
            z = [zeros(length(ys_list),1); filteratt(shock_list,t)];
        end
        if options_.kalman_algo > 200
            z = [1; z];
        end
        err = filteratt(yf_list,t) - betamat' * z;
        bet(:,:,t+1) = betamat;
        SecMom_New = SecMom + gain * (z * z' - SecMom);
        if min(abs(eig(SecMom_New))) < ridge
            betamat_New = betamat + gain * ((SecMom_New + ridge / ridge_adj * eye(smsize)) \ (z * err'));
        else
            betamat_New = betamat + gain * (SecMom_New \ (z * err'));
        end
        betama_ = [betamat_New(2:n_y+1,:)' betamat_New(n_y+2:end,:)' * RO];
        betama_ = betama_(:,reorder);
        EXP(:,vars) = A2 * betama_;
        T_New = - (A1 + EXP) \ LG;
        hit = any(abs(eig(T_New)) > options_.qz_criterium);
        if options_.corr == 0 || hit == 0
            R = - (A1 + EXP) \ dr_.jacobian(:,fw+1:end);
            T = T_New(dr_.order_var,dr_.order_var);
            R = R(dr_.order_var,:);
            TT(:,:,t+1) = T;
            RR(:,:,t+1) = R;
            QQ = R * Q * R';
            SecMom = SecMom_New;
            betamat = betamat_New;
        else
            TT(:,:,t+1) = T;
            RR(:,:,t+1) = R;
        end
        mu_new_1 = (A1 + EXP) \ ( (LG + A1 + A2_full) * SState - A2 * betamat(1,:)');
        mu =  mu_new_1(dr_.order_var);
        mus(:,t+1) = mu;
        if hit > 0
            fprintf('%4i',t);
            projection_facility = projection_facility + 1;
        end
    else
        TT(:,:,t+1) = T;
        RR(:,:,t+1) = R;
    end
end

% This piece works only when P becomes stationary, which can't happen
% under learning anyway
K_s = K(:,:,t);
%-------------------------------------------------- 
KwithoutT_s = KwithoutT(:,:,t);
%-------------------------------------------------- 
iF_s = iF(:,:,t);
P_s = P(:,:,t+1);
if t<smpl
	t_steady = t+1;
	P  = cat(3,P(:,:,1:t),repmat(P(:,:,t),[1 1 smpl-t_steady+1]));
	iF = cat(3,iF(:,:,1:t),repmat(inv(P_s(mf,mf)),[1 1 smpl-t_steady+1]));
	L  = cat(3,L(:,:,1:t),repmat(T-K_s*Z,[1 1 smpl-t_steady+1]));
	K  = cat(3,K(:,:,1:t),repmat(T*P_s(:,mf)*iF_s,[1 1 smpl-t_steady+1]));
end
while t<smpl
    t=t+1;
    v(:,t) = Y(:,t) - a(mf,t) - trend(:,t);
    a(:,t+1) = T*a(:,t) + K_s*v(:,t);
%-------------------------------------------------- 
filteratt(:,t+1)=a(:,t) + KwithoutT_s*v(:,t);
if options_.kalman_algo > 200
    filteratt(1,t+1) = 1;
end
%-------------------------------------------------- 
    aK(1,:,t+1) 	 	= a(:,t+1);
    for jnk=2:nk,
        aK(jnk,:,t+jnk) 	 	= T^(jnk-1)*a(:,t+1);
    end
end

% Smoothed variables should work OK
t = smpl+1;
while t>d+1 && t>2
	t = t-1;
    r(:,t-1) = transpose(Z)*iF(:,:,t)*v(:,t) + transpose(L(:,:,t))*r(:,t);
    alphahat(:,t)	= a(:,t) + P(:,:,t)*r(:,t-1);
	etahat(:,t)		= Q * RR(:,:,t)'*r(:,t);
end
if d
	r0 = zeros(mm,d); r0(:,d) = r(:,d);
	r1 = zeros(mm,d);
	for t = d:-1:2
    	r0(:,t-1) = transpose(Linf(:,:,t))*r0(:,t);
		r1(:,t-1) = transpose(Z)*(iFinf(:,:,t)*v(:,t)-transpose(Kstar(:,:,t))*r0(:,t)) + transpose(Linf(:,:,t))*r1(:,t);
		alphahat(:,t)	= a(:,t) + Pstar(:,:,t)*r0(:,t-1) + Pinf(:,:,t)*r1(:,t-1);
		etahat(:,t)		= Q * RR(:,:,t)'*r0(:,t);
	end
	r0_0 = transpose(Linf(:,:,1))*r0(:,1);
	r1_0 = transpose(Z)*(iFinf(:,:,1)*v(:,1)-transpose(Kstar(:,:,1))*r0(:,1)) + transpose(Linf(:,:,1))*r1(:,1);
	alphahat(:,1)  	= a(:,1) + Pstar(:,:,1)*r0_0 + Pinf(:,:,1)*r1_0;
	etahat(:,1)		= Q * RR(:,:,1)'*r0(:,1);
else
    r0 = transpose(Z)*iF(:,:,1)*v(:,1) + transpose(L(:,:,1))*r(:,1);
    alphahat(:,1)	= a(:,1) + P(:,:,1)*r0;
    etahat(:,1)	= Q * RR(:,:,1)'*r(:,1);
end

if projection_facility
    fprintf('\n %4i violations \n',projection_facility);
end

% Saving the beliefs and time-varying transmission mechanism

oo_.bet = bet;
oo_.TT = TT;
oo_.RR = RR;
oo_.mus = mus;

% Plotting the beliefs. Time axis has sample-specific start, therefore it
% is not included here

leg = cell(size(bet,1),1);
vars_X = [ys_list' shock_list'];
if options_.kalman_algo > 200
    leg{1} = 'const';
    for i = 2:length(leg)
        leg{i} = lgy_(dr_.order_var(vars_X(i-1)),:);
    end
else
    for i = 2:length(leg)
        leg{i} = lgy_(dr_.order_var(vars_X(i)),:);
    end
end   

for i = 1:size(bet,2)
    figure(i+20)
    temp = squeeze(bet(:,i,:));
    plot(temp')
    title(lgy_(dr_.order_var(yf_list(i)),:));
    legend('v6',leg)
end


