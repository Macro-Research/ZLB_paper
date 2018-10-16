function [alpha,beta,r,projectionFac_flag] =msv_learning(x,regressor,alphaOld,betaOld,rOld,gain)



thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);

lambda=10e-5;
% 
try
eig_min=min(abs(eigs(r)));
catch
eig_min=10e-15;
end



% 
if eig_min>lambda
theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';
else
    theta=thetaOld+gain*((r+lambda*eye(length(r))))^(-1)*yy*(x-thetaOld'*yy)';


end



theta=theta';
alpha=theta(:,1);
beta=theta(:,2:end);

projectionFac_flag=0;
 

end
