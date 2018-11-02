function [alpha,beta_1,beta_2,r,projectionFac_flag] =...
    ar2_learning(x,regressor,alphaOld,beta1_old,beta2_old,rOld,gain)



thetaOld=[alphaOld,beta1_old,beta2_old]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);

lambda=10e-5;
% 
% try
% eig_min=min(abs(eigs(r)));
% catch
% eig_min=10e-15;
% end



% 
%if eig_min>lambda
theta=thetaOld+gain*pinv(r)*yy*(x-thetaOld'*yy)';    
%else
 %   theta=thetaOld+gain*((r+lambda*eye(length(r))))^(-1)*yy*(x-thetaOld'*yy)';


%end



theta=theta';
alpha=theta(1);
beta_1=theta(2);
beta_2=theta(3);

projectionFac_flag=0;
 

end
