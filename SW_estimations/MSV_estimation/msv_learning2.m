function [theta r largestEig projectionFlag] =msv_learning2(x,regressor,thetaOld,rOld,gain,numBackward,backward_indices)

%thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
%  try
%  r=nearestSPD(r);
%  catch
%      r=100*eye(length(r));
%  end

theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';

msgstr=[];
msgid=[];
[msgstr, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
    r=1*eye(size(r,1));
    
theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';
end


% theta=theta';
% alpha=theta(:,1);
% beta=theta(2,3)';
% beta=theta(2:numBackward+1,backward_indices)';
% projectionFlag=0;
% largestEig=0;

% %largestEig=abs(eigs(A^(-1)*(B+C*beta^2),1));
% %     try
% largestEig=abs(eigs(beta,1));
%     catch
%     largestEig=999;
%     end
%     
% if largestEig>1
%     theta=thetaOld';
% % theta=zeros(14,13);
%  projectionFlag=1;
% end

% if largestEig>1
%     alpha=alphaOld;
%     beta=betaOld;
%     r=rOld;
%     projectionFac_flag=1;
% else 
%      projectionFac_flag=0;
% end

largestEig=0;projectionFlag=0;
end



