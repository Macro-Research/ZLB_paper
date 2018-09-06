function [alphaNew] = msv_learningUpdate(Y,X,alphaOld,rOld)
% gain=1/length(Y);
gain=0.05;



r=rOld+gain*(X(end)^2-rOld);
alphaNew=alphaOld+gain*r^(-1)*X(end)*(Y(end)-alphaOld*X(end));



end