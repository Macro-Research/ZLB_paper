clear;clc;close all;
rng(1);
load('MH_Candidate');
load('bounds.mat');
Ndraws=200000;
numVar=length(mode);
c=0.4; %
load('full_dataset.mat');
first_obs=200;last_obs=220;dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);
first_obs=200;last_obs=220;dataset=[gap_hp pinfobs robs];
dataset=dataset(first_obs:last_obs,:);

recursiveAverages=nan(Ndraws,numVar);
Nburn=round(Ndraws/2);
posteriorDraws=nan(Ndraws,numVar);
currentDraw=mvnrnd(mode,c*Sigma);
posteriorDraws(1,:)=currentDraw;
accept=0;
% objective=-likelihoodKalman(posteriorDraws(1,:),dataset);
[loglh,objective]=NKPC_objFcn(posteriorDraws(1,:)',1,bounds,dataset);

counter=0;
logposterior=objective*ones(Ndraws,1); 

for i=1:Ndraws
    
    currentDraw=mvnrnd(posteriorDraws(1,:),c*Sigma);
%     objectiveNew=-likelihoodKalman(currentDraw,dataset);
[loglh,objectiveNew]=NKPC_objFcn(currentDraw',1,bounds,dataset);
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
counter=counter+1;
if counter==500
     disp(['Acceptance Rate: ', num2str(acceptanceRate)]);
    disp(['Remaining Draws: ', num2str(Ndraws-i)]);
    counter=0;
end


end


figure;
title('Posterior Distributions');
for i=1:numVar;
subplot(3,5,i);hist(posteriorDraws(Nburn:Ndraws,i));hold on;plot([mode(i) mode(i)],[ylim],'r'); 
end

figure;
title('Recursive Averages');
 for i=1:numVar;
for j=1:Ndraws
    recursiveAverages(j,i)=mean(posteriorDraws(1:j,i));
    
end
   subplot(3,5,i);
   plot(recursiveAverages(:,i));
  
end
