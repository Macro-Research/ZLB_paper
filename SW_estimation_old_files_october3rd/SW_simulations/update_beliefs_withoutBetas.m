alpha_old=alpha_tt;
cc_old=cc_tt;    
%---------
thetaOld=[alpha_tt(forward_indices)  cc_tt(forward_indices,:)];
%--------
[theta,rr_tt ,largest_eig1(tt),largest_eig2(tt),pr_flag(tt)] =...
    msv_learning_withoutBetas(S_filtered(tt,forward_indices)',...
    [1;S_filtered(tt,18:24)'],thetaOld,rr_tt,gain,...
  forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2);

alpha_tt(forward_indices)=theta(1,:)';
cc_tt(forward_indices,:)=theta(2:end,:)';
