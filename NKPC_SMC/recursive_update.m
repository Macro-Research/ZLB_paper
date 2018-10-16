function [ alpha,beta,r] = recursive_update( X,t,alphaOld,betaOld,rOld)

x=X(t);
x_1=X(t-1);
x0=X(1);


alpha = alphaOld + 1/t * (x-alphaOld);

r = rOld + 1/t * (  (t-1)/t * (x-alphaOld)^2-rOld);

beta = betaOld + 1/(r*t) * ( (x- alphaOld )*...
       (x_1 + x0/t - ((t-1)^2+3*(t-1)+1)/t^2 * alphaOld - 1/t^2 *x )...
       -(t-1)/t*betaOld*(x-alphaOld)^2  );






end

