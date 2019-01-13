clear;clc;close all;
dynare NKPC_REE_Estimation;

imp_coll(:,:,1)=y_eps_y;
imp_coll(:,:,2)=y_eps_pi;
imp_coll(:,:,3)=pi_eps_y;
imp_coll(:,:,4)=pi_eps_pi;


shocks={'\eta_y','\eta_{pi}','\eta_r'};
vars={'y','pi','r','u_y','u_{\pi}'};
figure('Name','Impulse Responses-REE without switching');
index=0;


for ii=1:2
    for jj=1:2
        index=index+1;
        subplot(2,2,index);
    plot(imp_coll(2:end,index),'lineWidth',3,'color','blue');
%     hold on;
%     plot(imp_coll(2:end,2,index),'lineWidth',3,'color','green');
         title([vars(ii) shocks(jj)]);
    end
end
% legend('normal regime','zlb regime');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'NKPC_ree_init_REE_IR','-dpdf');
