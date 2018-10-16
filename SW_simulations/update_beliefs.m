alpha_old=alpha_tt;
beta_old=beta_tt;
cc_old=cc_tt;
rr_old=rr_tt;

thetaOld=[alpha_tt(forward_indices)...
    beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];

[alpha_tt,beta_tt,cc_tt,rr_tt] =...
    msv_learning4(XX(forward_indices,tt),...
    [1;XX(backward_indices,tt-1);XX(numEndo+1:end,tt)],...
    thetaOld,rr_old,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_crit,eig_eps);



eig_max(tt)=max(largest_eig1(tt),largest_eig2(tt));