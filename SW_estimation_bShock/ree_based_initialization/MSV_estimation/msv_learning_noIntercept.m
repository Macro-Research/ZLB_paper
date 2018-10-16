function [theta r largestEig projectionFlag] =msv_learning2(x,regressor,thetaOld,rOld,gain,numBackward,backward_indices)
msgstr=[];
msgid=[];
lastwarn('Success');

projectionFlag=0;largestEig=0;

 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);

lambda=10e-5;

try
eig_min=min(abs(eigs(r)));
catch
eig_min=10e-15;
end







    if eig_min>lambda
    theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';
    else
    theta=thetaOld'+gain*((r+lambda*eye(length(r))))^(-1)*yy*(x-thetaOld*yy)';
    end
    
    
%alpha_tt=theta(1,:)';
beta_tt=theta(1:numBackward,:)';


%cc_tt=theta(numBackward+2:end,:)';

try
largestEig=abs(eigs(beta_tt,1));
catch
    largestEig=1.01;
end

[msgstr, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
r=rOld;
theta=thetaOld';
projectionFlag=1;
elseif largestEig>1
    r=rOld;
theta=thetaOld';
projectionFlag=1;

end


    



end



