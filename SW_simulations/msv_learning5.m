function [alpha beta r,largestEig,projectionFac_flag] =msv_learning(x,regressor,alphaOld,betaOld,rOld,gain,eig_crit,tt)


regressand=x;

regression=(regressor'*regressor)\regressor'*regressand;

alpha=regression(1);
beta=regression(2);
r=0;
largestEig=0;
projectionFac_flag=0;

 
end

