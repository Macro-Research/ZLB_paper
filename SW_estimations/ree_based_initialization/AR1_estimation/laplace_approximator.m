function[laplace]=laplace_approximator(likl,mode,sigma)

N=length(mode);

laplace=likl-N*log(2*pi)/2-log(det(sigma)^(0.5));

%laplace=exp(likl)/((1/((2*pi)^(N/2)))*det(sigma)^(0.5));
end



