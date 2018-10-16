alpha_old=alpha_tt;
beta_old=beta_tt;
cc_old=cc_tt;    
%---------
thetaOld=[alpha_tt(forward_indices)...
    beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];

%auxiliary=detrend(S_filtered(1:tt,:));
%regressor=[1;auxiliary(tt-1,backward_indices)';auxiliary(tt,18:24)'];
%regressand=auxiliary(tt,forward_indices)';
 regressor=[1;S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'];
 regressand=S_filtered(tt,forward_indices)';
%--------
[theta,rr_tt ,largest_eig1(tt),largest_eig2(tt),pr_flag(tt)] =...
    msv_learning4(regressand,regressor,thetaOld,rr_tt,gain,...
  forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,ergodic_states);
%---
% [theta,rr_tt ,largest_eig(tt),pr_flag(tt)] =...
%     msv_learning2(S_filtered(tt,forward_indices)',...
%     [1;S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,numBackward,backward_indices);

%--------
alpha_tt(forward_indices)=theta(1,:)';
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
% alpha_tt=alpha_old;
% beta_tt=beta_old;
% beta_tt=beta_old;
 
 
%    beta_old=beta_tt;cc_old=cc_tt;    
%    thetaOld=[beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];
%   [theta,rr_tt,largestEig(tt),pr_flag(tt)]=msv_learning2(S_filtered(tt,forward_indices)',[S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,numBackward,backward_indices);
% 
%  beta_tt(forward_indices,backward_indices)=theta(1:numBackward,:)';
% 
%  cc_tt(forward_indices,:)=theta(numBackward+1:end,:)';
%------------------------------------------