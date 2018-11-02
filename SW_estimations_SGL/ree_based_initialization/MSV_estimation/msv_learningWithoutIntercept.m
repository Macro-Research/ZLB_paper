function [theta r largestEig projectionFlag] =msv_learningWithoutIntercept(x,regressor,thetaOld,rOld,gain,numBackward,backward_indices)

 yy=regressor;
 
r=rOld+gain*(yy*yy'-rOld);

theta=thetaOld'+gain*r^(-1)*yy*(x-thetaOld*yy)';

msgstr=[];
msgid=[];
[msgstr, msgid] = lastwarn;


largestEig=0;projectionFlag=0;
end



