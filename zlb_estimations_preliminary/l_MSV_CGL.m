function [alpha beta r] =l_MSV_CGL(x,regressor,alphaOld,betaOld,rOld,gain)
%function [alpha beta r] =msv_learning(x,regressor,alphaOld,betaOld,rOld,gain)
thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';


theta=theta';
alpha=theta(:,1);
beta=theta(:,2:end);




end



