clear;clc;%close all;
load('kf_output_workspace.mat');
periods=40;
numPeriods=50;
numVar=5;numShocks=3;
shocks={'\eta_y','\eta_{pi}','\eta_r'};
vars={'y','pi','r'};
beta1=zeros(numVar,numVar);
index=0;
for jj=T-numPeriods:T-1
    index=index+1;
    beta1(1:2,1:2)=diag(beta_tt(:,jj));
gamma1_1=A1_inv*(B1+C1*beta1^2);gamma3_1=A1_inv*D1;
gamma1_2=A2_inv*(B2+C2*beta1^2);gamma3_2=A1_inv*D2;

gamma1_avg=pp_filtered(jj)*gamma1_1+(1-pp_filtered(jj))*gamma1_2;
gamma3_avg=pp_filtered(jj)*gamma3_1+(1-pp_filtered(jj))*gamma3_2;

imp11(index,:,:,:)=impulse_response(parameters,gamma1_1,gamma3_1,periods);
imp22(index,:,:,:)=impulse_response(parameters,gamma1_2,gamma3_2,periods);
imp33(index,:,:,:)=impulse_response(parameters,gamma1_avg,gamma3_avg,periods);
end

refLine1=linspace(1,periods,periods);refLine1=repmat(refLine1,[numPeriods 1]);
refLine2=linspace(1,numPeriods,numPeriods);refLine2=repmat(refLine2,[periods 1])';




index=0;
 figure('Name','Impulse Responses-AR(1) Model');
for ii=1:2
    for jj=1:2
      
        index=index+1;
        subplot(2,2,index);
%         subplot(numVar-2,numShocks,index);
  plot3(refLine1',refLine2',imp33(:,:,ii,jj)','color','red','lineWidth',0.5);
%   plot3(refLine1',refLine2',imp11(:,:,ii,jj)','color','red','lineWidth',0.5);
%  hold on;
% plot3(refLine1',refLine2',imp22(:,:,ii,jj)','color','black','lineWidth',0.5);% plot(imp11(10,:,ii,jj),'lineWidth',3,'color','blue');
% hold on;
% plot(imp22(10,:,ii,jj),'lineWidth',3,'color','green');

       title([vars(ii) shocks(jj)]);
    end
end
% legend('normal regime','zlb regime');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_ree_init_AR1_IR_timeVarying','-dpdf');

% index=0;
% for ii=1:numVar-2
%     for jj=1:numShocks
%         figure;
%         index=index+1;
%         subplot(numVar-2,numShocks,index);
% plot3(refLine1',refLine2',imp33(:,:,ii,jj)','k','lineWidth',1);
% legend('weighted average');
%        title([vars(ii) shocks(jj)]);
%     end
% end

