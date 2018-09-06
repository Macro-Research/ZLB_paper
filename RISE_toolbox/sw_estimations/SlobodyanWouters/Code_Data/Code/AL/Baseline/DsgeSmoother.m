function [alphahat,etahat,epsilonhat,ahat,aK,aKM,SteadyState,trend_coeff,filteratt,bet] = ...
    DsgeSmoother(xparam1,gend,Y)
% stephane.adjemian@cepremap.cnrs.fr [09-07-2004]
%
% Adapted from mj_optmumlik.m
global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_ xparam1_test
global trend_coeff_

betamat = options_.betamat;
SecondMoments = options_.SecondMoments;
ys_list = options_.ys_list;
shock_list = options_.shock_list;
yf_list = options_.yf_list;
reorder = options_.reorder;  

alphahat 	= [];
epsilonhat	= [];
etahat		= [];
bet         = [];
nobs 		= size(options_.varobs,1);
smpl        = size(Y,2);

Q = Sigma_e_;
for i=1:estim_params_.nvx
	k =estim_params_.var_exo(i,1);
	Q(k,k) = xparam1(i)*xparam1(i);
end
offset = estim_params_.nvx;
if estim_params_.nvn
	H = zeros(nobs,nobs);
	for i=1:estim_params_.nvn
		k = estim_params_.var_endo(i,1);
		H(k,k) = xparam1(i+offset)*xparam1(i+offset);
	end
end	
offset = offset+estim_params_.nvn;
for i=1:estim_params_.ncx
	k1 =estim_params_.corrx(i,1);
	k2 =estim_params_.corrx(i,2);
	Q(k1,k2) = xparam1(i+offset)*sqrt(Q(k1,k1)*Q(k2,k2));
	Q(k2,k1) = Q(k1,k2);
end
offset = offset+estim_params_.ncx;

if estim_params_.nvn && estim_params_.ncn 
	for i=1:estim_params_.ncn
		k1 = options_.lgyidx2varobs(estim_params_.corrn(i,1));
		k2 = options_.lgyidx2varobs(estim_params_.corrn(i,2));
		H(k1,k2) = xparam1(i+offset)*sqrt(H(k1,k1)*H(k2,k2));
		H(k2,k1) = H(k1,k2);
	end
	offset = offset+estim_params_.ncn;
end	
for i=1:estim_params_.np
	assignin('base',deblank(estim_params_.param_names(i,:)),xparam1(i+offset));
end
Sigma_e_ = Q;
%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
[T,R,SteadyState] = dynare_resolve;
if options_.kalman_algo > 100
    if mod(options_.kalman_algo,100) == 3
      [betamat,SecondMoments,R_beta] = BetaFromTR(T,R,Q);
      options_.betamat = betamat;
      options_.SecondMoments = SecondMoments;
      options_.R_beta = R_beta;
    end
   [T,R] = TRFromBeta(betamat,T);
end
if options_.loglinear == 1
	constant = log(SteadyState(bayestopt_.mfys));
else
	constant = SteadyState(bayestopt_.mfys);
end
trend_coeff = zeros(nobs,1);
% betamat is just the slopes, need to add proper intercepts
if bayestopt_.with_trend == 1
	trend_coeff = zeros(nobs,1);
	nx1 = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn;
	for i=1:nobs
    	trend_coeff(i) = evalin('base',bayestopt_.trend_coeff{i});
	end
	trend = constant*ones(1,gend)+trend_coeff*(1:gend);
else
	trend = constant*ones(1,gend);
end
start = options_.presample+1;
np    = size(T,1);
mf    = bayestopt_.mf;
% ------------------------------------------------------------------------------
%  3. Initial condition of the Kalman filter
% ------------------------------------------------------------------------------
% 
%  C'est ici qu'il faut déterminer Pinf et Pstar. Si le modèle est stationnaire,
%  alors il suffit de poser Pstar comme la solution de l'éuation de Lyapounov et
%  Pinf=[].
% 
if options_.lik_init == 1		% Kalman filter
  Pstar = lyapunov_symm(T,R*Q*transpose(R));
  Pinf	= [];
elseif options_.lik_init == 2 % Old Diffuse Kalman filter
  Pstar = 10*eye(np);
  Pinf	= [];
elseif options_.lik_init == 3 % Diffuse Kalman filter
  Pstar = zeros(np,np);
  ivs = bayestopt_.i_T_var_stable;
  Pstar(ivs,ivs) = lyapunov_symm(T(ivs,ivs),R(ivs,:)*Q* ...
			transpose(R(ivs,:)));
  Pinf  = bayestopt_.Pinf;
end
% -----------------------------------------------------------------------------
%  4. Kalman smoother
% -----------------------------------------------------------------------------
filteratt=0;
if estim_params_.nvn
	if options_.kalman_algo == 1
		[alphahat,epsilonhat,etahat,ahat] = DiffuseKalmanSmootherH1(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
		if all(alphahat(:)==0)
			[alphahat,epsilonhat,etahat,ahat] = DiffuseKalmanSmootherH3(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
		end
	elseif options_.kalman_algo == 3
		[alphahat,epsilonhat,etahat,ahat] = DiffuseKalmanSmootherH3(T,R,Q,H,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
	end	
else
	if options_.kalman_algo == 1
		[alphahat,etahat,ahat,aK,filteratt] = DiffuseKalmanSmoother1(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
		if all(alphahat(:)==0)
			[alphahat,etahat,ahat,aK,filteratt] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
        end
        aKM = [];
	elseif options_.kalman_algo == 3
		[alphahat,etahat,ahat,aK,filteratt] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
        aKM = [];
    elseif options_.kalman_algo > 600
            [alphahat,etahat,ahat,aK,aKM,filteratt,bet] = ...
                DiffuseKalmanSmootherKMAL1_options_const(T,R,Q,Pinf,Pstar,Y,SteadyState,nobs,np,smpl,mf,xparam1(end));
%             [alphahat,etahat,ahat,aK,aKM,filteratt,bet] = ...
%                 DiffuseKalmanSmootherKMAL1_options_const(T,R,Q,Pinf,Pstar,Y,SteadyState,nobs,np,smpl,mf,xparam1(end-1:end));
    elseif options_.kalman_algo > 300
%         [alphahat,etahat,ahat,xxx,filteratt,bet] = ...
%             DiffuseKalmanSmootherKL1_options(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf,xparam1(end-2:end));
%         [alphahat,etahat,ahat,xxx,filteratt,bet] = ...
%             DiffuseKalmanSmootherKL1_options(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf,xparam1(end));
        [alphahat,etahat,ahat,aK,aKM,filteratt,bet] = ...
            DiffuseKalmanSmootherKL1_options(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf,xparam1(end));
%             [alphahat,etahat,ahat,xxx,filteratt,bet] = ...
%                 DiffuseKalmanSmootherKL1_options(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf,xparam1(end-1:end));
    elseif options_.kalman_algo > 200
      [alphahat,etahat,ahat,aK,filteratt,bet] = ...
          DiffuseKalmanSmootherL1_options_const(T,R,Q,Pinf,Pstar,Y,SteadyState,nobs,np,smpl,mf,xparam1(end));
      aKM = [];
    elseif options_.kalman_algo > 100
        [alphahat,etahat,ahat,aK,filteratt,bet] = ...
            DiffuseKalmanSmootherL1_options(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf,xparam1(end));
        aKM = [];
      if all(alphahat(:)==0)
          [alphahat,etahat,ahat] = DiffuseKalmanSmoother3(T,R,Q,Pinf,Pstar,Y,trend,nobs,np,smpl,mf);
      end
	end		
end