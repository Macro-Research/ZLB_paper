clear;clc;close all;
load('MC_MS_MSV_results.mat');
d_REE=d_REE(1:2,1:2);
d_regime1=d_regime1(1:2,1:2);
size=[1 3 3];
J=6;
I=3;

figure;
index=0;
for jj=1:J
    for ii=1:I
        index=index+1;
        subplot(J,I,index);
        hist(learning(:,jj,ii));
        
%         if jj>1
%             hold on;
%            
% %             line([d_REE(jj-1,ii) d_REE(jj-1,ii)],[0 200],'lineStyle','--','lineWidth',3,'color','red');
% %             hold on;
% %              line([d_regime1(jj-1,ii) d_regime1(jj-1,ii)],[0 200],'lineStyle','--','lineWidth',3,'color','red');
%         else 
%             hold on;
%              line([0 0],[0 200],'lineStyle','--','lineWidth',3,'color','red');
%         end
        
    end
end
      

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MC_MS_MSV_withLags','-dpdf');
