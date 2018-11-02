function [ beta r,largestEig,projectionFac_flag] =msv_learningwithoutIntercept(x,regressor,betaOld,rOld,gain)

thetaOld=[betaOld]';
 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);


theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';



msgstr=[];
msgid=[];
[msgstr, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
    r=5*eye(size(r,1));
    
theta=thetaOld+gain*r^(-1)*yy*(x-thetaOld'*yy)';
end

% theta=theta';
% alpha=theta(:,1);
% beta=theta(:,2:end);
beta=theta;


 largestEig=0;projectionFac_flag=0;
end



