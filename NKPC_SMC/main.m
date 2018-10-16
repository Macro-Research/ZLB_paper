clear;clc;%close all;
tic
load('full_dataset.mat');
% load('simulated_dataset.mat');
load('bounds.mat');
% first_obs=44;last_obs=length(dy);dataset=[gap_hp pinfobs robs];
first_obs=39;last_obs=length(gap_hp);dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);

tune.npara = 13;    
tune.npart = 5000;   
tune.nphi  = 50;    
tune.lam   = 2;    

tune.c      = 0.25;             
tune.acpt   = 0.25;             
tune.trgt   = 0.25;            
tune.alp    = 0.9;               

tune.phi = (1:1:tune.nphi)';
tune.phi = ((tune.phi-1)/(tune.nphi-1)).^tune.lam;

f = @(x1,x2) NKPC_objFcn(x1,x2,bounds,dataset);

parasim = zeros(tune.nphi, tune.npart, tune.npara); % parameter draws
wtsim   = zeros(tune.npart, tune.nphi);        % weights
zhat    = zeros(tune.nphi,1);                  % normalization constant
nresamp = 0; % record # of iteration resampled

csim    = zeros(tune.nphi,1); % scale parameter
ESSsim  = zeros(tune.nphi,1); % ESS
acptsim = zeros(tune.nphi,1); % average acceptance rate
rsmpsim = zeros(tune.nphi,1); % 1 if re-sampled

% % load('MH_Candidate');
% % c=0.4;
% % priorsim=mvnrnd(mode,c*Sigma,tune.npart); %initial proposal MH_Candidate instead of prior
priorsim       = prior_draws(tune.npart,bounds);
parasim(1,:,:) = priorsim;        % from prior

wtsim(:, 1)    = 1/tune.npart;    % initial weight is equal weights
zhat(1)        = sum(wtsim(:,1));

% Posterior values at prior draws
loglh   = zeros(tune.npart,1); %log-likelihood
logpost = zeros(tune.npart,1); %log-posterior

parfor i=1:1:tune.npart
    p0                     = priorsim(i,:)';
    [logpost(i), loglh(i)] = f(p0,tune.phi(1)); % likelihood
end

for i=2:1:tune.nphi
    disp(i)
    toc
    %-----------------------------------
    % (a) Correction
    %-----------------------------------
    % incremental weights
    incwt = exp((tune.phi(i)-tune.phi(i-1))*loglh);
   
    % update weights
    wtsim(:, i) = wtsim(:, i-1).*incwt;
    zhat(i)     = sum(wtsim(:, i));
    
    % normalize weights
    wtsim(:, i) = wtsim(:, i)/zhat(i);
    
    %-----------------------------------
    % (b) Selection
    %-----------------------------------
    ESS = 1/sum(wtsim(:, i).^2); % Effective sample size
    
    if (ESS < tune.npart/2)
        
        [id, m] = systematic_resampling(wtsim(:,i)'); %systematic resampling
        
        parasim(i-1, :, :) = squeeze(parasim(i-1, id, :));
        loglh              = loglh(id);
        logpost            = logpost(id);
        wtsim(:, i)        = 1/tune.npart; % resampled weights are equal weights
        nresamp            = nresamp + 1;
        rsmpsim(i)         = 1;
        
    end
    
    %--------------------------------------------------------
    % (c) Mutuation
    %--------------------------------------------------------
    % Adapting the transition kernel
    
    tune.c = tune.c*(0.95 + 0.10*exp(16*(tune.acpt-tune.trgt))/(1 + ...
             exp(16*(tune.acpt-tune.trgt))));
    
    % Calculate estimates of mean and variance
    para      = squeeze(parasim(i-1, :, :));
    wght      = repmat(wtsim(:, i), 1, tune.npara);
   
    tune.mu      = sum(para.*wght); % mean
    z            = (para - repmat(tune.mu, tune.npart, 1));
    tune.R       = (z.*wght)'*z;       % covariance
%     tune.R=nearestSPD(tune.R);%%%DOES THIS MAKE SENSE?
    tune.Rdiag   = diag(diag(tune.R)); % covariance with diag elements

    tune.Rchol   = chol(tune.R, 'lower');
    tune.Rchol2  = sqrt(tune.Rdiag);
 
     
    % Particle mutation (RWMH 2)
    temp_acpt = zeros(tune.npart,1); %initialize accpetance indicator
    
    parfor j = 1:tune.npart %iteration over particles
 
[ind_para, ind_loglh, ind_post, ind_acpt] = mutation_RWMH(para(j,:)', ...
                                         loglh(j), logpost(j), tune, i, f); 
        parasim(i,j,:) = ind_para;
        loglh(j)       = ind_loglh;
        logpost(j)     = ind_post;
        temp_acpt(j,1) = ind_acpt;
        
    end
    
    tune.acpt = mean(temp_acpt); % update average acceptance rate
    
    % store
    csim(i,:)    = tune.c; % scale parameter
    ESSsim(i,:)  = ESS; % ESS
    acptsim(i,:) = tune.acpt; % average acceptance rate
    
    % print some information
    if mod(i, 1) == 0
        
        para = squeeze(parasim(i, :, :));
        wght = repmat(wtsim(:, i), 1, 13);

mu  = sum(para.*wght);
sig = sum((para - repmat(mu, tune.npart, 1)).^2 .*wght);
sig = (sqrt(sig));

        % time calculation


      
        smctime = tic; % re-start clock
    end
end


figure;
for jj=1:tune.nphi 
    subplot(5,4,jj);
    auxiliary=sort(wtsim(:,jj));
    lower=auxiliary(round(0.05*length(auxiliary)));
    upper=auxiliary(round(0.95*length(auxiliary)));
    hist(wtsim(:,jj),100);
    if upper>lower
    xlim([lower upper]);
    end
end;
