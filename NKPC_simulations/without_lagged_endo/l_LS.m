function [alpha beta r,largestEig,projectionFac_flag] =l_LS(x,regressor,alphaOld,betaOld,rOld,gain)

thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
try
 r=nearestSPD(r);
catch
    r=10*eye(size(r,1));
end

theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';


theta=theta';
alpha=theta(:,1);
beta=theta(:,2:end);


%largestEig=abs(eigs(A^(-1)*(B+C*beta^2),1));

    try
largestEig=abs(eigs(beta,1));   
    catch
largestEig=999;
    end

if largestEig>.9999
    alpha=alphaOld;

    beta=betaOld;
    r=rOld;
    projectionFac_flag=1;
else 
     projectionFac_flag=0;
end

 % largestEig=0;projectionFac_flag=0;
end



