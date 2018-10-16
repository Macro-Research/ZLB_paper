function [alpha beta r,largestEig,projectionFac_flag] =msv_learning(x,regressor,alphaOld,betaOld,rOld,gain)
msgstr=[];
msgid=[];


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
    %theta=thetaOld+gain*(nearestSPD(r+lambda*eye(length(r))))^(-1)*yy*(x-thetaOld'*yy)';

        try
theta=thetaOld+gain*((r+lambda*eye(length(r))))^(-1)*yy*(x-thetaOld'*yy)';
        catch
    theta=thetaOld;
        end

end


[msgstr, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
    r=rOld;
    theta=thetaOld;
% theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';
end
% 
theta=theta';
alpha=theta(:,1);
beta=theta(:,2:end);


%largestEig=abs(eigs(A^(-1)*(B+C*beta^2),1));

%     try
% largestEig=abs(eigs(beta,1));   
%     catch
% largestEig=999;
%     end

% if largestEig>1
%     alpha=alphaOld;
%     beta=betaOld;
%     r=rOld;
%     projectionFac_flag=1;
% else 
%      projectionFac_flag=0;
% end

%  largestEig=0;projectionFac_flag=0;
end



