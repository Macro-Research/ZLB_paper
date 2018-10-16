alpha_old=alpha_tt;
beta_old=beta_tt;
cc_old=cc_tt;    
%---------
thetaOld=[alpha_tt(forward_indices)...
    beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];
%--------
[theta,rr_tt ,largestEig(tt),pr_flag(tt)] =...
    msv_learning2(S_filtered(tt,forward_indices)',...
    [1;S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,numBackward,backward_indices);
%--------
alpha_tt(forward_indices)=theta(1,:)';
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';

%alpha_tt(forward_indices)=zeros(7,1);
 
 
%    beta_old=beta_tt;cc_old=cc_tt;    
%    thetaOld=[beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];
%   [theta,rr_tt,largestEig(tt),pr_flag(tt)]=msv_learning2(S_filtered(tt,forward_indices)',[S_filtered(tt-1,backward_indices)';S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,numBackward,backward_indices);
% 
%  beta_tt(forward_indices,backward_indices)=theta(1:numBackward,:)';
% 
%  cc_tt(forward_indices,:)=theta(numBackward+1:end,:)';
%------------------------------------------