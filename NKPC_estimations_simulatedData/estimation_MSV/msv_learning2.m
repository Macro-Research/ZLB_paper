function [theta r largestEig projectionFlag] =msv_learning2(x,regressor,thetaOld,rOld,gain)

%thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
try
r=nearestSPD(r);
catch
r=5*eye(size(r,1));
end

theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';
msgstr=[];
msgid=[];
[msgstr,msgid]=lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
r=5*eye(size(r,1));
theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';
end
% theta=theta';
% alpha=theta(:,1);
beta=theta(2,3)';
projectionFlag=0;
largestEig=0;

% %largestEig=abs(eigs(A^(-1)*(B+C*beta^2),1));
    try
largestEig=beta;
    catch
    largestEig=999;
    end
%     
if largestEig>1
    theta=thetaOld';
 projectionFlag=1;
end

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



