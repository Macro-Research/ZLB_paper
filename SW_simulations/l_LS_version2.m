function [theta r] =l_LS_version2(x,regressor,thetaOld,rOld,gain)

%thetaOld=[alphaOld,betaOld]';
 yy=regressor;

 
 
lambda=10e-3;


r=rOld+gain*(yy*yy'-rOld);
r=triu(r);
r=(r+r')/2;
%eig_min=abs(eigs(r,1));
% 
% [VV DD]=eig(r);
% DD=diag(DD);
% DD(0<DD<lambda)=min(DD> lambda);
% DD(-lambda<DD<0)=max(DD(DD<-lambda));
% r=VV *diag(DD) / VV;
% r=real(r);
% r=1*eye(length(rOld));


%if eig_min>lambda
theta=thetaOld'+gain*pinv(r)*yy*(x-thetaOld*yy)';
%else
%theta=thetaOld'+gain*pinv(r+lambda)*yy*(x-thetaOld*yy)';
%end



end



