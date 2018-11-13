clear;clc;%close all;
load('param_init.mat');
param=param_init;

% msgstr=[];
% msgid=[];
% param=x;
% load('estimation_results.mat');
% names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
%     'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
%            'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
%            'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','mu_p','mu_w','rho_ga',...
%            'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r1','eta_r2', 'eta_p', 'eta_w',...
%            'gain','1-p_11','1-p_22','rbar_zlb'} ] ; 

%variable order: mc zcap rk k1   q c inve y lab pinf w r kp 
numVar=24;numShocks=7;numEndo=17;numExo=7;numBackward=7;numForward=7;
backward_indices=[6 7 8 10 11 12 13];
forward_indices=[3 5 6 7 9 10 11];
forward_names=[{'rk','q','c','i','l','pi','w'}];
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [param(1),param(1)] ;    %phi
parameters(7,:)  =[param(2),param(2)] ;      %sigma_c
parameters(8,:)  =  [param(3),param(3)];   %lambda 
parameters(9,:)  =  [param(4),param(4)] ; %xi_w
parameters(10,:) = [param(5),param(5)];        %sigma_l 
parameters(11,:) =  [param(6),param(6)];      %xi_p 
parameters(12,:) =  [param(7),param(7)] ;     %iota_w
parameters(13,:) =  [param(8),param(8)];         %iota_p
parameters(14,:) = [param(9),param(9)] ;   %psi
parameters(15,:) = [param(10),param(10)]; %phi_p

%policy related parameters

parameters(16,:)    =   [param(11),0];   %r_pi
parameters(17,:)    =  [param(12),0]; %rho
parameters(18,:)    =   [param(13),0];%0.0746;    %r_y
parameters(19,:)    =  [param(14),0];     %r_dy

%SS related parameters
parameters(20,:)    = [param(15),param(15)] ;    %pi_bar
parameters(21,:)    = [param(16),param(16)];       %beta_const
parameters(22,:)    = [param(17),param(17)];     %l_bar
parameters(23,:)    = [param(18),param(18)];         %gamma_bar
parameters(24,:)    = [param(19),param(19)] ;    %alpha

% %shock persistence
parameters(25,:) = [param(20),param(20)];      %rho_a
parameters(26,:) =[param(21),param(21)];  %rho_b
parameters(27,:) =[param(22),param(22)];   %rho_g
parameters(28,:) =[param(23),param(23)];   %rho_i
parameters(29,:) =  [param(24),0];   %rho_r
parameters(30,:) =  [param(25),param(25)];  %rho_p
parameters(31,:) =[param(26),param(26)];%rho_w 


% parameters(32,:) = [param(27),param(27)] ;    %mu_p 
% parameters(33,:) =[param(28),param(28)];    %mu_w
parameters(32,:) =  [param(27),param(27)]      ;  %rho_ga

%shock standard deviations
parameters(33,:)=[param(28),param(28)];  %sigma_a
parameters(34,:)=  [param(29),param(29)];  %sigma_b
parameters(35,:)= [param(30),param(30)]; %sigma_g
parameters(36,:)=  [param(31),param(31)] ;   %sigma_i
parameters(37,:)=[param(32),param(33)];   %sigma_r1 & sigma_r2
parameters(38,:)= [param(34),param(34)];  %sigma_p
parameters(39,:)=   [param(35),param(35)];   %sigma_w
gain=param(36);
% gain=0.025;
p_11=1-param(37);p_22=1-param(38); 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

[AA1, BB1, CC1, DD1, EE1 ,RHO1 ,FF1, GG1 ,E1 ,F1] =SW_sysmat_MSV_filter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2, RHO2, FF2 ,GG2, E2, F2]=SW_sysmat_MSV_filter(parameters(:,2));

beta_tt=0*eye(numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=0*ones(numEndo,numShocks);
 % rr_tt=5*eye(numBackward+numExo+1);
%  rr_tt=triu(rr_tt);
% rr_tt=(rr_tt+rr_tt')/2;
 load('initial_beliefs_msv.mat');
 
  beta_tt(forward_indices,backward_indices)=beta_init;
   cc_tt(forward_indices,:)=cc_init;
rr_tt=diag(diag(rr_init));
   
%rr_tt=diag(diag(rr_init));%rr_tt(rr_tt<0.01)=0.01;
 % rr_tt=triu(rr_tt);
 % rr_tt=(rr_tt+rr_tt')/2;
 %rr_tt=diag(diag(rr_init));

%  load('rr_init1.mat');
%  rr_tt=rr_init1;


E2(6)=param(39);%ss level of interest rate
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);
% Sigma1=diag(parameters(end-numShocks+1:end,1))^2;
% Sigma2=diag(parameters(end-numShocks+1:end,2))^2;
Sigma1=diag(parameters(33:39,1))^2;
Sigma2=diag(parameters(33:39,2))^2;

load('raf_dataset.mat');first_obs=71;last_obs=length(dy);
dataset=[dy dc dinve dw pinfobs robs labobs];
dataset=dataset(first_obs:last_obs,:);l=7;N=length(dataset);numVar=24;burnIn=06;
T=size(dataset,1);numObs=7;
startDate=datenum('01-01-1966');
endDate = datenum('01-12-2016');
Date=linspace(startDate,endDate,T);

obs_indices=[14 15 16 17 10 12 9];
obs_intercepts=[parameters(23,1),parameters(23,1),parameters(23,1),parameters(23,1),...
parameters(20,1),E1(6),parameters(22,1)]';    
obs_names=[{'Output Growth','Consumption Growth','Real Wage Growth','Real Investment Growth',...
    'Inflation','Fed Funds Rate','Hours Worked'}];


H1=0*eye(numObs);H2=0*eye(numObs);


S_fore11=zeros(numVar,1);S_fore12=zeros(numVar,1);S_fore21=zeros(numVar,1);S_fore22=zeros(numVar,1);
P_fore11=zeros(numVar,numVar);P_fore12=zeros(numVar,numVar);P_fore21=zeros(numVar,numVar);P_fore22=zeros(numVar,numVar);
v11=zeros(numObs,1);v12=zeros(numObs,1);v21=zeros(numObs,1);v22=zeros(numObs,1);
Fe11=zeros(numObs,numObs);Fe12=zeros(numObs,numObs);Fe21=zeros(numObs,numObs);Fe22=zeros(numObs,numObs);
S_upd11=zeros(numVar,1);S_upd12=zeros(numVar,1);S_upd21=zeros(numVar,1);S_upd22=zeros(numVar,1);
P_upd11=zeros(numVar,numVar);P_upd12=zeros(numVar,numVar);P_upd21=zeros(numVar,numVar);P_upd22=zeros(numVar,numVar);
ml11=0;ml12=0;ml21=0;ml22=0;
likl=zeros(T,1);
%pp_fore11=0;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_fore11=1;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_upd11=1;pp_upd12=0;pp_upd21=0;pp_upd22=0;
%pp_collapse1=regime(1);pp_collapse2=1-regime(1);
pp_collapse1=1;pp_collapse2=0;
S_collapse1=zeros(numVar,1);S_collapse2=zeros(numVar,1);
P_collapse1=eye(numVar);P_collapse2=eye(numVar);
q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
S_filtered=zeros(T,numVar);

pp_filtered=ones(T,1);


%----------------------------------------------------------------------------

update_matrices;

for tt=2:T
    
  
    x_tt=dataset(tt,:)';

%kalman block
S_fore11=gamma1_1*S_collapse1+gamma2_1;%
S_fore12=gamma1_2*S_collapse1+gamma2_2;%
%S_fore12(19)=S_fore12(18)-1.38;
S_fore12(12)=S_fore12(12)-E1(6)+E2(6);
S_fore21=gamma1_1*S_collapse2+gamma2_1;%
S_fore22=gamma1_2*S_collapse2+gamma2_2;%
%S_fore22(19)=S_fore22(18)-1.38;
S_fore22(12)=S_fore22(12)-E1(6)+E2(6);
%
P_fore11=gamma1_1*P_collapse1*gamma1_1'+gamma3_1*Sigma1*gamma3_1';%+gamma4_1*Sigma1*gamma4_1';
P_fore12=gamma1_2*P_collapse1*gamma1_2'+gamma3_2*Sigma2*gamma3_2';%+gamma4_2*Sigma2*gamma4_2';;  
P_fore21=gamma1_1*P_collapse2*gamma1_1'+gamma3_1*Sigma1*gamma3_1';%+gamma4_1*Sigma1*gamma4_1';;
P_fore22=gamma1_2*P_collapse2*gamma1_2'+gamma3_2*Sigma2*gamma3_2';%+gamma4_2*Sigma2*gamma4_2';;
%    
v11=x_tt-E1-F1*S_fore11;
v12=x_tt-E1-F2*S_fore12;
v21=x_tt-E1-F1*S_fore21;
v22=x_tt-E1-F2*S_fore22;
% %

% v11=x_tt-E2-F1*S_fore11;
% v12=x_tt-E2-F2*S_fore12;
% v21=x_tt-E2-F1*S_fore21;
% v22=x_tt-E2-F2*S_fore22;


Fe11=F1*P_fore11*F1'+H1;
Fe12=F2*P_fore12*F2'+H2;
Fe21=F1*P_fore21*F1'+H1;
Fe22=F2*P_fore22*F2'+H2;
%
S_upd11=S_fore11+P_fore11*F1'*(Fe11^(-1))*v11;
S_upd12=S_fore12+P_fore12*F2'*(Fe12^(-1))*v12;
S_upd21=S_fore21+P_fore21*F1'*(Fe21^(-1))*v21;
S_upd22=S_fore22+P_fore22*F2'*(Fe22^(-1))*v22;
%
P_upd11=(eye(numVar)-P_fore11*F1'*Fe11^(-1)*F1)*P_fore11;    
P_upd12=(eye(numVar)-P_fore12*F2'*Fe12^(-1)*F2)*P_fore12;  
P_upd21=(eye(numVar)-P_fore21*F1'*Fe21^(-1)*F1)*P_fore21;  
P_upd22=(eye(numVar)-P_fore22*F2'*Fe22^(-1)*F2)*P_fore22;  
%
pp_fore11=q_11*pp_collapse1;
pp_fore12=q_12*pp_collapse1;
pp_fore21=q_21*pp_collapse2;
pp_fore22=q_22*pp_collapse2;
%
ml11=(2*pi)^(-l/2)*det(Fe11)^(-0.5)*exp(-0.5*v11'*Fe11^(-1)*v11);
ml12=(2*pi)^(-l/2)*det(Fe12)^(-0.5)*exp(-0.5*v12'*Fe12^(-1)*v12);
ml21=(2*pi)^(-l/2)*det(Fe21)^(-0.5)*exp(-0.5*v21'*Fe21^(-1)*v21);
ml22=(2*pi)^(-l/2)*det(Fe22)^(-0.5)*exp(-0.5*v22'*Fe22^(-1)*v22);
% 
likl(tt)=ml11*pp_fore11+ml12*pp_fore12+ml21*pp_fore21+ml22*pp_fore22;
%
pp_upd11=(ml11*pp_fore11)/likl(tt);
pp_upd12=(ml12*pp_fore12)/likl(tt);
pp_upd21=(ml21*pp_fore21)/likl(tt);
pp_upd22=(ml22*pp_fore22)/likl(tt);
%
pp_collapse1=pp_upd11+pp_upd21;
pp_collapse2=pp_upd12+pp_upd22;
%
if pp_collapse1>10e-10;
    S_collapse1=(pp_upd11*S_upd11+pp_upd21*S_upd21)/pp_collapse1;
P_collapse1=(pp_upd11*(P_upd11+(S_collapse1-S_upd11)*(S_collapse1-S_upd11)')+...
    pp_upd21*(P_upd21+(S_collapse1-S_upd21)*(S_collapse1-S_upd21)'))/pp_collapse1;

else
    S_collapse1=(S_upd11+S_upd21)/2;
    P_collapse1=(P_upd11+P_upd21)/2;
%     S_collapse1=zeros(numVar,1);
%     P_collapse1=1*eye(numVar);
 end

 if pp_collapse2>10e-10;
P_collapse2=(pp_upd12*(P_upd12+(S_collapse2-S_upd12)*(S_collapse2-S_upd12)')+...
    pp_upd22*(P_upd22+(S_collapse2-S_upd22)*(S_collapse2-S_upd22)'))/pp_collapse2;
S_collapse2=(pp_upd12*S_upd12+pp_upd22*S_upd22)/pp_collapse2;
%=nearestSPD(P_collapse1);
%P_collapse2=nearestSPD(P_collapse2);
  else
     S_collapse2=(S_upd12+S_upd22)/2;
     P_collapse2=(P_upd12+P_upd22)/2;
     
%     S_collapse2=zeros(numVar,1);
%     P_collapse2=1*eye(numVar);
 end

%
S_filtered(tt,:)=pp_collapse1*S_collapse1+pp_collapse2*S_collapse2;
pp_filtered(tt)=pp_collapse1;
% %------------------------------------------
update_beliefs;
update_matrices;
% check_eigenvalue;
% revert_if_explosive;
% update_matrices;
% check_eigenvalue;  
%-----------------------------------------    
learning_filtered(tt,:,:)=theta;

obs_forecasts(:,tt)=obs_intercepts + pp_fore11*S_fore11(obs_indices)+...
    pp_fore12*S_fore12(obs_indices)+...
    pp_fore21*S_fore21(obs_indices)+pp_fore22*S_fore22(obs_indices);

forecast_errors(:,tt)=pp_fore11*v11+pp_fore12*v12+pp_fore21*v21+pp_fore22*v22;

largest_eig2(tt)=abs(eigs(gamma1_2,1));

rr_aggregate(tt,:,:)=rr_tt;
end

likl=-sum(log(likl(burnIn+1:end)))
yylim=[min(dataset)',max(dataset)'];

figure('Name','MSV Forecast Errors','units','normalized','outerposition',[0 0 1 1]);
for jj=1:length(obs_indices)

subplot(7,1,jj);
plot(Date,obs_forecasts(jj,:),'-','color','black','lineWidth',1);
hold on;
plot(Date,dataset(:,jj),'-','color','red','lineWidth',1);
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
title(obs_names(jj));
grid on;
grid minor;
ylim([yylim(jj,:)]);
if jj==length(obs_indices)
    legend('forecast','observable');
end

% subplot(7,2,2*jj)
% plot(forecast_errors(jj,2:end),'color','black','lineWidth',3);
% % hold on;
% % plot(dataset(:,jj),'--','color','red');
% %title(obs_names(jj));
% grid on;
% grid minor;
% xlim([0 length(dataset(:,jj))]);
% if jj==length(obs_indices)
%     legend('forecast error');
% end   
end

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_forecast_errors','-dpdf');  
%--------------------------------------------------
figure('Name','MSV Eigenvalues','units','normalized','outerposition',[0 0 1 1]);
subplot(4,1,1);
plot(Date(2:end),largest_eig1(2:end),'color','black','lineWidth',3);
hold on;plot(Date(2:end),ones(length(largest_eig1(2:end)),1),'color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
title('Largest eigenvalue: normal regime');
subplot(4,1,2);
plot(Date(2:end),largest_eig2(2:end),'color','black','lineWidth',3);
hold on;
plot(Date(2:end),ones(length(largest_eig2(2:end)),1),'color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
title('Largest eigenvalue: zlb regime','lineWidth',3);
subplot(4,1,3);
plot(Date(2:end),average_eig(2:end),'color','black','lineWidth',3);
title('Largest eigenvalue: weighted average','lineWidth',3);
hold on;
plot(Date(2:end),ones(length(average_eig(2:end)),1),'color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
subplot(4,1,4);

plot(Date(2:end),pr_flag(2:end),'color','blac');title('projection facility activity','lineWidth',3);
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_eigenvalues','-dpdf');  
%------------------------------------------------------

figure('Name','MSV Regime Probability','units','normalized','outerposition',[0 0 1 1]);
plot(Date,1-pp_filtered,'lineWidth',3,'color','blue');
hold on;
plot(Date,S_filtered(:,12),'lineWidth',3);
legend('Regime Probability','Interest Rate');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_regimeProb','-dpdf'); 
%----------------------------------------------------
figure('Name','MSV Learning Coef-Intercept','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,1,jj)'),'color','black');

text(Date(end),learning_filtered(end,1,jj),forward_names(index));
hold on;
end
  xlim([startDate endDate])
  ylim([-0.5, 0.3]);%same scale as AR(1) model
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_alphas','-dpdf');  



figure('Name','Implied variances of learning terms','units','normalized','outerposition',[0 0 1 1]);
for jj=1:size(diag(rr_tt));
    subplot(5,3,jj);
    plot(rr_aggregate(2:end,jj,jj),'lineWidth',3);
    
    
    
end
%------------------------------------------------
%------------------------------------------------
figure('Name','MSV Learning Coef-lagged consumption','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,2,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,2,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_laggedCons','-dpdf'); 
%------------------------------------------
figure('Name','MSV Learning Coef-lagged investment','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,3,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,3,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_laggedInve','-dpdf'); 
%----------------
figure('Name','MSV Learning Coef-lagged output','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,4,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,4,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_laggedOutput','-dpdf'); 
%--------------------------------
figure('Name','MSV Learning Coef-lagged inflation','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,5,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,5,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_laggedInfl','-dpdf'); 
%=----------------------------------------------
figure('Name','MSV Learning Coef-lagged wage','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,6,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,6,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_laggedCons','-dpdf'); 


%--------------------------------------------------
figure('Name','MSV Learning Coef-b shock','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,10,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,10,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_bShock','-dpdf');

figure('Name','MSV Learning Coef-g shock','units','normalized','outerposition',[0 0 1 1]);
index=0;
for jj=1:length(forward_indices)
    index=index+1;
plot(Date(2:end),squeeze(learning_filtered(2:end,11,jj)'),'color','black');
hold on;
text(Date(end),learning_filtered(end,11,jj),forward_names(index));
hold on;
end
plot(Date,ones(T,1),'--','color','red');
  xlim([startDate endDate])
  datetick('x','yyyy','keeplimits');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_msv_learning_gShock','-dpdf');
%---------------------
% output_gap=S_filtered(:,8)-param(10)*S_filtered(:,18);
% gap_growth=diff(output_gap(2:end));
% figure('Name','Output gap and growth rate');
% subplot(2,1,1);
% plot(output_gap);
% subplot(2,1,2);
% plot(gap_growth);
% 
% 
% figure('Name','inflation intercept');
% subplot(2,1,1);
% plot(S_filtered(:,10));
% subplot(2,1,2);
% plot(learning_filtered(2:end,10,1));