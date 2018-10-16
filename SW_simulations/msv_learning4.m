function [alpha_tt beta_tt cc_tt rr_tt] =msv_learning4(x,regressor,thetaOld,rOld,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_crit,eig_eps)
numVar=17;
numBackward=length(backward_indices);
lambda=10e-5;
lastwarn('Success');
 yy=regressor;
 
% try 
rr_tt=rOld+gain*(yy*yy'-rOld);
theta=thetaOld'+gain*rr_tt^(-1)*yy*(x-thetaOld*yy)';

% catch
%rr_tt=increase_eigenvalue(rr_tt,lambda,lambda);
%rr_tt=nearestSPD(rr_tt);
%theta=thetaOld'+gain*(pinv(rr_tt)*yy)*(x-thetaOld*yy)';

% rr_tt=rOld;
% theta=thetaOld';
% end



alpha_tt=zeros(numVar,1);
alpha_tt(forward_indices)=theta(1,:)';

beta_tt=zeros(numVar,numVar);
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
beta_tt=reduce_eigenvalue(beta_tt,eig_crit,eig_eps);

cc_tt=zeros(numVar,7);
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';


% [msgstr, msgid] = lastwarn;

% if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
%     theta=thetaOld';
%     alpha_tt(forward_indices)=theta(1,:)';
%     beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
%     cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
%     rr_tt=rOld;
% elseif abs(eigs(beta_tt,1))>eig_crit
%      theta=thetaOld';
%     alpha_tt(forward_indices)=theta(1,:)';
%     beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
%     cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
%     rr_tt=rOld;
% end

end



