clear;clc;close all;
load('nkpc_estimation_results.mat');
myirfs0=irf(mm,'irf_periods',40);

myirfs1=irf(mm,'irf_periods',40,'irf_type','girf');


startDate=datenum('01-01-1966');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,203);

regime_prob=double(mm.filtering.filtered_regime_probabilities.regime_2);

figure;
area(Date,regime_prob(2:end));
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_sigmaPoint_regimeProb','-dpdf');


figure;
plot(Date,regime_prob(2:end),'lineWidth',3);
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
hold on;
plot(Date,robs(end-202:end),'--');
xlim([startDate endDate])
datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_sigmaPoint_regime','-dpdf');


figure;
plot(myirfs0.eps_y.y(:,1));
hold on;
plot(myirfs0.eps_y.y(:,2));
legend('normal regime','zlb regime');

std_y=mm.estimation.posterior_maximization.mode(14);
std_pinf=mm.estimation.posterior_maximization.mode(15);

imp_coll(:,:,1)=double(myirfs0.eps_y.y);
imp_coll(:,:,2)=double(myirfs0.eps_pi.y);
imp_coll(:,:,3)=double(myirfs0.eps_y.pinf);
imp_coll(:,:,4)=double(myirfs0.eps_pi.pinf);
shocks={'\eta_y','\eta_{pi}','\eta_r'};
vars={'y','pi','r','u_y','u_{\pi}'};
figure('Name','Impulse Responses-REE with switching');
index=0;


for ii=1:2
    for jj=1:2
        index=index+1;
        subplot(2,2,index);
    plot(imp_coll(2:end,1,index),'lineWidth',3,'color','blue');
    hold on;
    plot(imp_coll(2:end,2,index),'lineWidth',3,'color','green');
         title([vars(ii) shocks(jj)]);
    end
end
legend('normal regime','zlb regime');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_ree_init_REE_MS_IR','-dpdf');
