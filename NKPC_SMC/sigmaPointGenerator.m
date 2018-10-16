function [sigmaPoints,weights]=sigmaPointGenerator(mu_uncond,p_uncond,gamma1,gamma2,gamma3,sigma)

numVar=length(gamma1);
numSigmaPoint=numVar*2+1;
lambda=3-numVar;
weights=nan(numSigmaPoint,1);

% [mu_uncond,p_uncond]=kalman_init(gamma1,gamma2,gamma3,sigma);

sigmaPoints=nan(numSigmaPoint,numVar);
sigmaPoints(1,:)=mu_uncond;
auxiliary=(sqrt((numVar+lambda)*p_uncond));
for mm=2:numVar+1   
    
    sigmaPoints(mm,:)=mu_uncond+auxiliary(mm-1,mm-1);
end

for mm=numVar+2:numSigmaPoint
    sigmaPoints(mm,:)=mu_uncond-auxiliary(mm-numVar-1,mm-numVar-1);
end

weights(1)=lambda/(lambda+numVar);
for i=2:numSigmaPoint
    weights(i)=1/(2*(numVar+lambda));
end

end