function [theta r,largestEig,projectionFlag] =l_LS_version2(x,regressor,thetaOld,rOld,gain,size)

%thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';


% theta=theta';
% alpha=theta(:,1);
%beta=theta(size(1)+1:size(1)+size(2),:);


projectionFlag=0;
largestEig=0;
%largestEig=abs(eigs(beta,1));
%if largestEig>1
 %   theta=thetaOld';
 %   projectionFlag=1;
%end
% 
% if largestEig>1
%     alpha=alphaOld;
%     beta=betaOld;
%     r=rOld;
%     projectionFac_flag=1;
% else 
%      projectionFac_flag=0;
% end

% largestEig=0;projectionFac_flag=0;
end



