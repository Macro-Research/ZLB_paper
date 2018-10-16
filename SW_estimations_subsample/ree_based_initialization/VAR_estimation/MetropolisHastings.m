clear;clc;%close all;
rng(1);
tic
load('MH_Candidate');
Ndraws=20000;
numVar=length(mode);
c=0.1; %
recursiveAverages=nan(Ndraws,numVar);
Nburn=round(Ndraws/2);
crit=0.01;
% for jj=1:length(mode)
%     if Sigma(jj,jj)<crit
%         Sigma(jj,jj)=crit;
%     end
% end
        
% posteriorDraws=nan(Ndraws,numVar);
 currentDraw=mvnrnd(mode,c*Sigma);
%currentDraw=mvnrnd(mode,0.01*eye(length(mode)));
posteriorDraws(1,:)=currentDraw;
accept=0;
objective=-likelihood(posteriorDraws(1,:));
counter=0;
logposterior=objective*ones(Ndraws,1);

for i=1:Ndraws
%toc
% disp('REMAINING:');
% disp(Ndraws-i);

    currentDraw=mvnrnd(posteriorDraws(1,:),c*Sigma);
    objectiveNew=-likelihood(currentDraw);
    alpha=min(1,exp(objectiveNew-objective));
    u=rand(1);

    if u<=alpha
        posteriorDraws(i+1,:)=currentDraw;
        accept=accept+1;
        objective=objectiveNew;
        logposterior(i+1)=objectiveNew;
    else
        posteriorDraws(i+1,:)=posteriorDraws(i,:);
        logposterior(i+1)=objective;
    end

acceptanceRate=accept/i;
if acceptanceRate>0.3
    c=c+crit;
elseif c>0.05 
    c=c-crit;
end


disp(['ACCEPTANCE RATE:',num2str(acceptanceRate)]);
disp(['Scale Coefficient:',num2str(c)]);
counter=counter+1;
if counter==100
     disp(['Acceptance Rate: ', num2str(acceptanceRate)]);
    disp(['Remaining Draws: ', num2str(Ndraws-i)]);

    counter=0;
end


end


figure('Name','Posterior Distributions','units','normalized','outerposition',[0 0 1 1]);

title('Posterior Distributions');
for i=1:numVar;
subplot(6,7,i);hist(posteriorDraws(Nburn:Ndraws,i)); 

end
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'nkpc_mc_posteriors','-dpdf'); 

figure;
title('Recursive Averages');
 for i=1:numVar;
for j=1:Ndraws
    recursiveAverages(j,i)=mean(posteriorDraws(1:j,i));
    
    
end
   subplot(6,7,i);
   plot(recursiveAverages(:,i),'lineWidth',2);
   
 end

  disp('POSTERIOR MEAN')
  Nburn=1;
  mh_conf_sig=0.95;
 mean(posteriorDraws(Nburn+1:end,:))
posteriorDist=posteriorDraws(Nburn+1:end,:);
posteriorDist=sort(posteriorDist);
hpdDraws=round((1-mh_conf_sig)*length(posteriorDist));
kk=zeros(hpdDraws,13);


for paraInd=1:length(mode)
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
 
 
toc