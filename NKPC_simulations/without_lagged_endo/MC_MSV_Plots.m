clear;clc;close all;
dynare NKPC_REE_Simulation;
policyMatrix=zeros(5,2);
policyMatrix(4:5,1:2)=oo_.dr.ghu(4:5,1:2)';
clearvars -except policyMatrix;
load('MC_MSV_results.mat');

J=5;
I=2;

figure;
index=0;
for jj=1:J
    for ii=1:I
        index=index+1;
        subplot(J,I,index);
        hist(learning(:,jj,ii));
              hold on;
line([policyMatrix(jj,ii) policyMatrix(jj,ii)],[0 400],'lineStyle','--','lineWidth',3,'color','red');
    end
end

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MC_MSV_withoutLags','-dpdf');
      
