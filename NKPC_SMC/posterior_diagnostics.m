close all;

histGrid=50;

for i=1:tune.npara
    index=1;
    figure;
    for j=5:5:50;
    subplot(5,2,index);
    index=index+1;
    
    histogram(parasim(j,:,i),histGrid,'facealpha',0.9);
    hold on;
    histogram(parasim(j-1,:,i),histGrid,'facealpha',0.6);
    if j==2
        legend('target','proposal');
    end
   
    end
    
end

figure
for i=1:tune.npara
    subplot(3,5,i);
    hist(parasim(end,:,i),histGrid);
end

%extract final stage distributions

posteriorDist=reshape(parasim(end,:,:),[tune.npart tune.npara]);

%Compute mean, HPD interval and quantiles
mh_conf_sig=0.9;
 disp('POSTERIOR MEAN')
 mean(posteriorDist)
 posteriorDist=sort(posteriorDist);
hpdDraws=round((1-mh_conf_sig)*length(posteriorDist));
kk=zeros(hpdDraws,13);

for paraInd=1:tune.npara
    jj=length(posteriorDist)-hpdDraws-2;

for ii=1:hpdDraws
    kk(ii,paraInd)=posteriorDist(jj,paraInd)-posteriorDist(ii,paraInd);
    jj=jj+1;
end
[kmin,idx]=min(kk(:,paraInd));
hpd_interval(paraInd,:)=[posteriorDist(idx,paraInd) posteriorDist(idx,paraInd)+kmin];
post_deciles(paraInd,:)=posteriorDist([round(0.05*length(posteriorDist(:,paraInd)))...
    round(0.2*length(posteriorDist(:,paraInd)))...
    round(0.3*length(posteriorDist(:,paraInd)))...
    round(0.4*length(posteriorDist(:,paraInd)))...
    round(0.5*length(posteriorDist(:,paraInd)))...
    round(0.6*length(posteriorDist(:,paraInd)))...
    round(0.7*length(posteriorDist(:,paraInd)))...
    round(0.8*length(posteriorDist(:,paraInd)))...
    round(0.95*length(posteriorDist(:,paraInd)))],paraInd);
end