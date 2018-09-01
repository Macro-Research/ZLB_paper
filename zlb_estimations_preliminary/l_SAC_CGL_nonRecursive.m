function[alpha,beta]=cgl_learning_nonRecursive(y,alphaOld,gain)
%constant gain sample autocorrelation learning with AR(1) rule (univariate)
%non-recursive form.
N=length(y);
num=0;denum=0;
newObs=y(end);

alpha=gain*newObs+(1-gain)*alphaOld;

for tt=1:N
    if tt==1
    denum=denum+gain*(1-gain)^(N+1-tt)*(y(tt)-alpha)*(y(tt)-alpha);
    else 
    denum=denum+gain*(1-gain)^(N+1-tt)*(y(tt)-alpha)*(y(tt)-alpha);
    num=num+gain*(1-gain)^(N+1-tt)*(y(tt)-alpha)*(y(tt-1)-alpha);
    end
end




beta=sqrt(1-gain)*num/denum;


end