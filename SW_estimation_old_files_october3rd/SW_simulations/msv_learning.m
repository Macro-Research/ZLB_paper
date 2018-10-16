function [alpha beta r,largestEig,projectionFac_flag] =msv_learning(x,regressor,alphaOld,betaOld,rOld,gain,eig_crit)
largestEig=0;projectionFac_flag=0;

thetaOld=[alphaOld,betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);
lambda=10e-5;

r=triu(r);
r=(r+r')/2;

% [VV DD]=eig(r);
% DD=diag(DD);
% DD(DD<lambda)=min(DD(DD>lambda));
% r=VV* diag(DD) / VV;
% try
% eig_min=min(abs(eigs(r)));
% catch
% eig_min=10e-15;
% end

%theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';

    %if eig_min>lambda
     theta=thetaOld+gain*(pinv(r)*yy)*(x-thetaOld'*yy)';
   % else
   % theta=thetaOld+gain*(((r+lambda*eye(length(r))))\yy)*(x-thetaOld'*yy)';
   % end

msgstr=[];
msgid=[];
[msgstr, msgid] = lastwarn;

%if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
%    r=99*eye(size(r,1));
    
%theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';
%end

alpha=theta(1);
beta=theta(2);

% if abs(beta)>eig_crit
%     beta=betaOld;
%     alpha=alphaOld;
%     r=rOld;
%     projectionFac_flag=1;
% end


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

 
end



