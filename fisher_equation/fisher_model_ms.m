clear;clc;%close all;
%Fisher equation with two (Markov) regimes, msv learning. 
%rng(2);
N=1000;
phi_pinf1=.899;
phi_pinf2=2;
rho=0.9;
eta_sigma=0.1;
eta=normrnd(0,eta_sigma,[N 1]);
p_11=0.95;p_22=0.95;
Q=[p_11,1-p_11;1-p_22,p_22];
regime=zeros(N,1);
regime(1)=1;
r=zeros(N,1);
pinf=zeros(N,1);


learning=nan(N,2);
learningCovariance=nan(N,2,2);

ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
beta_ree=1/(ergodic_states(1)*phi_pinf1+ergodic_states(2)*phi_pinf2-rho);
alpha_ree=0;
E_stabilityMarkov=rho/(phi_pinf1*ergodic_states(1)+phi_pinf2*ergodic_states(2));
E_stabilityRegime1=rho/(phi_pinf1);
E_stabilityRegime2=rho/(phi_pinf2);
if E_stabilityMarkov<1
    disp('MARKOV CHAIN E-STABILITY SATISFIED');
else disp('MARKOV CHAIN E-STABILITY NOT SATISFIED');
end

if E_stabilityRegime1<1
    disp('REGIME 1 E-STABILITY SATISFIED');
else disp('REGIME 1 E-STABILITY NOT SATISFIED');
end

if E_stabilityRegime2<1
    disp('REGIME 2 E-STABILITY SATISFIED');
else disp('REGIME 2 E-STABILITY NOT SATISFIED');
end

beta_tt=beta_ree;
alpha_tt=alpha_ree;
rr_tt=eye(2);

forecast_errors=zeros(N,1);
forecast=zeros(N,1);
for jj=2:N
    gain=0.05;
    regime(jj)=findRegime(regime(jj-1),p_11,p_22);
    r(jj)=rho*r(jj-1)+eta(jj);
     learning(jj,:)=[alpha_tt,beta_tt];%forecast coef for t
     learningCovariance(jj,:,:)=rr_tt;
    
     forecast_errors(jj-1)=pinf(jj-1)-alpha_tt+beta_tt*r(jj-1);
%      if jj>10
   [alpha_tt beta_tt rr_tt] =...
       l_LS(pinf(jj-1),[1 r(jj-1)]',alpha_tt,beta_tt,rr_tt,gain);%forecast for t+1
%      end
    pinf(jj)=(regime(jj)*(1/phi_pinf1)+(1-regime(jj))*(1/phi_pinf2))...
    *((alpha_tt+beta_tt*rho*r(jj))+r(jj));


    
   
end

figure('Name','Fisher eqn-learning coef','units','normalized','outerposition',[0 0 1 1]);

% subplot(3,1,1);
% plot(learning(:,1),'lineWidth',3);
% hold on;
% plot(ones(N,1)*alpha_ree,'--');
% title('alpha');
subplot(2,1,1);
plot(learning(:,2),'lineWidth',3);
hold on;
plot(ones(N,1)*beta_ree,'--');
%hold on;
%plot(ones(N,1)*(1/(phi_pinf1-rho)),'--');
hold on;
plot(ones(N,1)*(1/(phi_pinf2-rho)),'--');
legend('learning coef','unconditional','regime 1','regime 2');
% hold on;
% area(regime);
title('beta');
subplot(2,1,2);
area(regime);
% hold on;
% plot(forecast_errors,'lineWidth',3);
title('regime');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fisher_simulation2_learningCoef','-dpdf'); 

figure('Name','simulated inflation','units','normalized','outerposition',[0 0 1 1]);
plot(pinf,'lineWidth',3);


fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fisher_simulation2_pinf','-dpdf'); 


% figure;
% subplot(2,1,1);
% plot(forecast_errors,'lineWidth',3);
% title('forecast errors');
% subplot(2,1,2);
% area(regime);
% title('regime');
% 
% figure;
% subplot(2,1,1);
% plot(pinf,'lineWidth',3);
% title('Inflation');
% subplot(2,1,2);
% area(regime);
% title('regime');
% figure;
% plotyy(linspace(0,N,N),regime,...
%     linspace(0,N,N),forecast_errors);
% title('forecast errors');

% figure;
% subplot(2,2,1);
% plot(learningCovariance(:,1,1),'lineWidth',3);
% subplot(2,2,2);
% plot(learningCovariance(:,1,2),'lineWidth',3);
% subplot(2,2,3);
% plot(learningCovariance(:,2,1),'lineWidth',3);
% subplot(2,2,4);
% plot(learningCovariance(:,2,2),'lineWidth',3);
