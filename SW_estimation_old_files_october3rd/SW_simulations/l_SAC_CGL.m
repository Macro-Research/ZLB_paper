function [alpha beta r] =l_SAC_CGL(x,xOld,alphaOld,betaOld,rOld,gain)
%constant gain sample autocorrelation learning with AR(1) rule, recursive
%form.
% alpha= alphaOld+gain*(x-alphaOld);
% 
% r = rOld+gain*((x-alphaOld)^2-rOld);
% 
% beta= betaOld+ gain*r^(-1)*((x-alphaOld)*(xOld-alphaOld)-betaOld*(x-alphaOld)^2);

alpha= alphaOld+gain*(x-alphaOld);

r = rOld+gain*((xOld-alphaOld)^2-rOld);

beta= betaOld+ gain*r^(-1)*((xOld-alphaOld))*((x-alpha)-betaOld*(xOld-alphaOld));



end