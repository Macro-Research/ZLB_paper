function[alpha,beta]=cgl_alpha_beta(y,alphaOld)

N=length(y);
num=0;denum=0;
newObs=y(end);
gain=0.05;
alpha=gain*newObs+(1-gain)*alphaOld;
for i=1:N-1
    num=num+gain*(1-gain)^(N-i)*(y(i)-alpha)*(y(i+1)-alpha);
end

for i=1:N
    denum=denum+gain*(1-gain)^(N-i+1)*(y(i)-alpha)*(y(i)-alpha);
end


beta=sqrt(1-gain)*num/denum;


end