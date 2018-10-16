function [theta r largest_eig1 largest_eig2 projectionFlag] =msv_learning3(x,regressor,thetaOld,rOld,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_thr,ergodic_states)
numEndo=17;
numShocks=7;
largest_eig1=0;
largest_eig2=0;
projectionFlag=0;
numBackward=length(backward_indices);

lastwarn('Success');
lambda=10e-5;
projectionFlag=0;

 yy=regressor;


 
r=rOld+gain*(yy*yy'-rOld);

theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';    

 
beta_tt=zeros(numEndo,numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=zeros(numEndo,numShocks);


beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';


gamma1_1_tilde=AA1^(-1)*(BB1+CC1*beta_tt^2);
gamma1_2_tilde=AA2^(-1)*(BB2+CC2*beta_tt^2);

try
largest_eig1=abs(eigs(beta_tt,1));
catch
    largest_eig1=nan;
end

try
ev1=abs(eigs(gamma1_1_tilde,1));
ev2=abs(eigs(gamma1_2_tilde,1));
 largest_eig2=ergodic_states(1)*ev1+ergodic_states(2)*ev2;
catch

    largest_eig2=nan;
end
% 
% [msgstr, msgid] = lastwarn;

% if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
% r=rOld;
% [VV DD]=eig(r);
% DD=diag(DD);
% DD(DD<lambda)=min(DD(DD>lambda));
% r=VV* diag(DD) / VV;
% r=real(r);
% 
% theta=thetaOld';
%  projectionFlag=1;
% 
% elseif largest_eig1>eig_thr
%     r=rOld;
% theta=thetaOld';
% projectionFlag=1;

% elseif largest_eig2>eig_thr
%     
%         r=rOld;
% theta=thetaOld';
% projectionFlag=1;

% end



%beta_tt=zeros(numEndo,numEndo);
% beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
% gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);gamma1_1_tilde=gamma1_1_tilde(1:17,1:17);
% gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);gamma1_2_tilde=gamma1_2_tilde(1:17,1:17);
% try
% largest_eig1=abs(eigs(beta_tt,1));
% ev1=abs(eigs(gamma1_1_tilde,1));
% ev2=abs(eigs(gamma1_2_tilde,1));
% largest_eig2=max(ev1,ev2);
% catch
%    largest_eig1=1.01;
%    largest_eig2=1.01;
% end
    



end



