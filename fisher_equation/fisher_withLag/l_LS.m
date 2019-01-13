function [aa,bb,r,largestEig,projectionFac_flag] =l_LS(x,regressor,aOld,bOld,rOld,gain)
%function for least squares updating.
thetaOld=[aOld,bOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';


theta=theta';
aa=theta(:,1);
bb=theta(:,2:end);


projectionFac_flag=0;
if bb>1
    aa=aOld;
    bb=bOld;
    r=rOld;
projectionFac_flag=1;
end
largestEig=bb;

%largestEig=abs(eigs(A^(-1)*(B+C*beta^2),1));
% largestEig=abs(eigs(beta,1));

% if largestEig>1
%     alpha=alphaOld;
%     beta=betaOld;
%     r=rOld;
%     projectionFac_flag=1;
% else 
%      projectionFac_flag=0;
% end

%largestEig=0;projectionFac_flag=0;
end



