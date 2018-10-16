
clear;clc;%close all;
load('param_init.mat');
param=param_init;
forward_indices=[3 5 6 7 9 10 11];
numVar=24;numShocks=7;numEndo=13;numExo=7;
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

% parameters(12,:) =  [0,0] ;     %iota_w
% parameters(13,:) =  [0,0];         %iota_p
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
parameters(29,:) =  [param(24),param(24)];   %rho_r
parameters(30,:) =  [param(25),param(25)];  %rho_p
parameters(31,:) =[param(26),param(26)];%rho_w 

% parameters(32,:) = [0,0] ;    %mu_p 
% parameters(33,:) =[0,0];    %mu_w
parameters(32,:) = [param(27),param(27)] ;    %mu_p 
parameters(33,:) =[param(28),param(28)];    %mu_w
parameters(34,:) =  [param(29),param(29)]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[param(30),param(30)];  %sigma_a
parameters(36,:)=  [param(31),param(31)];  %sigma_b
parameters(37,:)= [param(32),param(32)]; %sigma_g
parameters(38,:)=  [param(33),param(33)] ;   %sigma_i
parameters(39,:)=[param(34),param(35)];   %sigma_r1 & sigma_r2
parameters(40,:)= [param(36),param(36)];  %sigma_p
parameters(41,:)=   [param(37),param(37)];   %sigma_w
gain=param(38);
p_11=1-param(39);p_22=1-param(40); 
%p_11=1;p_22=0;
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

[AA1,BB1,CC1,DD1,EE1,E1,F1]=SW_sysmat_VAR_filter(parameters(1:end,1));
%E1(6)=1;E2(6)=0;%ss level of interest rate
[AA2,BB2,CC2,DD2,EE2,E2,F2]=SW_sysmat_VAR_filter(parameters(1:end,2));
E2(6)=param(41);%ss level of interest rate
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);
Sigma1=diag(parameters(end-numShocks+1:end,1))^2;
Sigma2=diag(parameters(end-numShocks+1:end,2))^2;
load('raf_dataset.mat');first_obs=147;last_obs=length(dy);
dataset=[dy dc dinve dw pinfobs robs labobs];
dataset=dataset(first_obs:last_obs,:);l=7;N=length(dataset);numVar=24;burnIn=6;
T=size(dataset,1);numObs=7;
alpha1=0*ones(numVar,1);
beta1_tt=0*eye(numVar);
beta2_tt=0*eye(numVar);
rr=100*repmat(eye(3),[1 1 numVar]);
%rr=100*ones(1,1,numVar);
%  beta1_tt(3,3)=0.98;beta1_tt(5,5)=0.8 ;beta1_tt(6,6)=0.98;beta1_tt(7,7)=0.98;beta1_tt(9,9)=0.98;
%  beta1_tt(10,10)=0.6;beta1_tt(11,11)=0.98;
load('AR1_initial_beliefs.mat');
for jj=1:length(beta1_init)
    %rr(1,1,forward_indices(jj))=rr_init(2,2,jj);
    beta1_tt(forward_indices(jj),forward_indices(jj))=beta1_init(jj);
    beta2_tt(forward_indices(jj),forward_indices(jj))=beta2_init(jj);
    rr(:,:,forward_indices(jj))=rr_init(:,:,jj);
end
 





H1=0*eye(numObs);H2=0*eye(numObs);



S_fore11=zeros(numVar,T);
S_fore12=zeros(numVar,T);
S_fore21=zeros(numVar,T);
S_fore22=zeros(numVar,T);

P_fore11=zeros(numVar,numVar,T);
P_fore12=zeros(numVar,numVar,T);
P_fore21=zeros(numVar,numVar,T);
P_fore22=zeros(numVar,numVar,T);

v11=zeros(numObs,1);v12=zeros(numObs,1);v21=zeros(numObs,1);v22=zeros(numObs,1);

Fe11=zeros(numObs,numObs,1);
Fe12=zeros(numObs,numObs,1);
Fe21=zeros(numObs,numObs,1);
Fe22=zeros(numObs,numObs,1);

S_upd11=zeros(numVar,T);
S_upd12=zeros(numVar,T);
S_upd21=zeros(numVar,T);
S_upd22=zeros(numVar,T);

P_upd11=zeros(numVar,numVar,T);
P_upd12=zeros(numVar,numVar,T);
P_upd21=zeros(numVar,numVar,T);
P_upd22=zeros(numVar,numVar,T);

S_collapse1=zeros(numVar,T);S_collapse2=zeros(numVar,T);

P_collapse1=repmat(eye(numVar),[1 1 T]);
P_collapse2=repmat(eye(numVar),[1 1 T]);

ml11=0;ml12=0;ml21=0;ml22=0;
likl=zeros(T,1);

pp_fore11=1;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_upd11=1;pp_upd12=0;pp_upd21=0;pp_upd22=0;

pp_collapse1=1;pp_collapse2=0;

q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
S_filtered=zeros(T,numVar);
pp_filtered=ones(T,1);


%----------------------------------------------------------------------------
gamma1_1=AA1_inv*(BB1+CC1*(beta1_tt^2+beta2_tt));
gamma2_1=AA1_inv*CC1*(eye(numVar)+beta1_tt)*alpha1;
gamma3_1=AA1_inv*DD1;
gamma4_1=AA1_inv*CC1*beta1_tt*beta2_tt;

gamma1_2=AA2_inv*(BB2+CC2*(beta1_tt^2+beta2_tt));
gamma2_2=AA2_inv*CC2*(eye(numVar)+beta1_tt)*alpha1;
gamma3_2=AA2_inv*DD2;
gamma4_2=AA2_inv*CC2*beta1_tt*beta2_tt;



%---------------------------------------

pr_flag=zeros(T,1);
for tt=3:T
    


    x_tt=dataset(tt,:)';

%kalman block
S_fore11(:,tt)=gamma1_1*S_collapse1(:,tt-1)+gamma4_1*S_collapse1(:,tt-2)+gamma2_1;%
S_fore12(:,tt)=gamma1_2*S_collapse1(:,tt-1)+gamma4_2*S_collapse1(:,tt-2)+gamma2_2;%
S_fore12(12,tt)=S_fore12(12,tt)-E1(6)+E2(6);
S_fore21(:,tt)=gamma1_1*S_collapse2(:,tt-1)+gamma4_1*S_collapse2(:,tt-2)+gamma2_1;%
S_fore22(:,tt)=gamma1_2*S_collapse2(:,tt-1)+gamma4_1*S_collapse2(:,tt-2)+gamma2_2;%
S_fore22(12,tt)=S_fore22(12,tt)-E1(6)+E2(6);
%
P_fore11(:,:,tt)=gamma1_1*P_collapse1(:,:,tt-1)*gamma1_1'+gamma4_1*P_collapse1(:,:,tt-2)*gamma4_1'+gamma3_1*Sigma1*gamma3_1';%+gamma4_1*Sigma1*gamma4_1';
P_fore12(:,:,tt)=gamma1_2*P_collapse1(:,:,tt-1)*gamma1_2'+gamma4_2*P_collapse1(:,:,tt-2)*gamma4_2'+gamma3_2*Sigma2*gamma3_2';%+gamma4_2*Sigma2*gamma4_2';  
P_fore21(:,:,tt)=gamma1_1*P_collapse2(:,:,tt-1)*gamma1_1'+gamma4_1*P_collapse2(:,:,tt-2)*gamma4_1'+gamma3_1*Sigma1*gamma3_1';%+gamma4_1*Sigma1*gamma4_1';
P_fore22(:,:,tt)=gamma1_2*P_collapse2(:,:,tt-1)*gamma1_2'+gamma4_2*P_collapse2(:,:,tt-2)*gamma4_2'+gamma3_2*Sigma2*gamma3_2';%+gamma4_2*Sigma2*gamma4_2';
%    
v11=x_tt-E1-F1*S_fore11(:,tt);
v12=x_tt-E1-F2*S_fore12(:,tt);
v21=x_tt-E1-F1*S_fore21(:,tt);
v22=x_tt-E1-F2*S_fore22(:,tt);

%
Fe11=F1*P_fore11(:,:,tt)*F1'+H1;
Fe12=F2*P_fore12(:,:,tt)*F2'+H2;
Fe21=F1*P_fore21(:,:,tt)*F1'+H1;
Fe22=F2*P_fore22(:,:,tt)*F2'+H2;
%
S_upd11(:,tt)=S_fore11(:,tt)+P_fore11(:,:,tt)*F1'*(Fe11^(-1))*v11;
S_upd12(:,tt)=S_fore12(:,tt)+P_fore12(:,:,tt)*F2'*(Fe12^(-1))*v12;
S_upd21(:,tt)=S_fore21(:,tt)+P_fore21(:,:,tt)*F1'*(Fe21^(-1))*v21;
S_upd22(:,tt)=S_fore22(:,tt)+P_fore22(:,:,tt)*F2'*(Fe22^(-1))*v22;
%
P_upd11(:,:,tt)=(eye(numVar)-P_fore11(:,:,tt)*F1'*Fe11^(-1)*F1)*P_fore11(:,:,tt);    
P_upd12(:,:,tt)=(eye(numVar)-P_fore12(:,:,tt)*F2'*Fe12^(-1)*F2)*P_fore12(:,:,tt);  
P_upd21(:,:,tt)=(eye(numVar)-P_fore21(:,:,tt)*F1'*Fe21^(-1)*F1)*P_fore21(:,:,tt);  
P_upd22(:,:,tt)=(eye(numVar)-P_fore22(:,:,tt)*F2'*Fe22^(-1)*F2)*P_fore22(:,:,tt);  
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

pp_upd11=(ml11*pp_fore11)/likl(tt);
pp_upd12=(ml12*pp_fore12)/likl(tt);
pp_upd21=(ml21*pp_fore21)/likl(tt);
pp_upd22=(ml22*pp_fore22)/likl(tt);
%
pp_collapse1=pp_upd11+pp_upd21;
pp_collapse2=pp_upd12+pp_upd22;
%
if pp_collapse1>10e-10;
    S_collapse1(:,tt)=(pp_upd11*S_upd11(:,tt)+pp_upd21*S_upd21(:,tt))/pp_collapse1;
P_collapse1(:,:,tt)=(pp_upd11*(P_upd11(:,:,tt)+(S_collapse1(:,tt)-S_upd11(:,tt))*...
    (S_collapse1(:,tt)-S_upd11(:,tt))')+...
    pp_upd21*(P_upd21(:,:,tt)+(S_collapse1(:,tt)-S_upd21(:,tt))*...
    (S_collapse1(:,tt)-S_upd21(:,tt))'))/pp_collapse1;

else
    S_collapse1(:,tt)=(S_upd11(:,tt)+S_upd21(:,tt))/2;
    P_collapse1(:,:,tt)=(P_upd11(:,:,tt)+P_upd21(:,:,tt))/2;

 end

 if pp_collapse2>10e-10;
P_collapse2(:,:,tt)=(pp_upd12*(P_upd12(:,:,tt)+(S_collapse2(:,tt)-S_upd12(:,tt))*...
    (S_collapse2(:,tt)-S_upd12(:,tt))')+...
    pp_upd22*(P_upd22(:,:,tt)+(S_collapse2(:,tt)-S_upd22(:,tt))*...
    (S_collapse2(:,tt)-S_upd22(:,tt))'))/pp_collapse2;
S_collapse2(:,tt)=(pp_upd12*S_upd12(:,tt)+pp_upd22*S_upd22(:,tt))/pp_collapse2;

  else
     S_collapse2=(S_upd12+S_upd22)/2;
     P_collapse2=(P_upd12+P_upd22)/2;

 end

%
S_filtered(tt,:)=pp_collapse1*S_collapse1(:,tt)+pp_collapse2*S_collapse2(:,tt);
pp_filtered(tt)=pp_collapse1;


largest_eig1(tt)=0;
largest_eig2(tt)=0;
average_eig(tt)=0;
        pr_flag(tt)=0;
        beta_old=beta1_tt;
        alpha_old=alpha1;
        
%         if tt>2
%         for jj=[forward_indices]
%  
%  [alpha1(jj),beta1_tt(jj,jj),beta2_tt(jj,jj),rr(:,:,jj),pr_flag(tt)] =...
%            ar2_learning(S_filtered(tt,jj)',[1,S_filtered(tt-1,jj),S_filtered(tt-2,jj)]',...
%         alpha1(jj),beta1_tt(jj,jj),beta2_tt(jj,jj),rr(:,:,jj),gain);
%         end
%         end
    


gamma1_1=AA1_inv*(BB1+CC1*(beta1_tt^2+beta2_tt));
gamma2_1=AA1_inv*CC1*(eye(numVar)+beta1_tt)*alpha1;
gamma3_1=AA1_inv*DD1;
gamma4_1=AA1_inv*CC1*beta1_tt*beta2_tt;

gamma1_2=AA2_inv*(BB2+CC2*(beta1_tt^2+beta2_tt));
gamma2_2=AA2_inv*CC2*(eye(numVar)+beta1_tt)*alpha1;
gamma3_2=AA2_inv*DD2;
gamma4_2=AA2_inv*CC2*beta1_tt*beta2_tt;


for jj=forward_indices
learning_filtered(jj,tt-1,:)=[alpha1(jj),beta1_tt(jj,jj)]; 
end

end


likl=-sum(log(likl(burnIn+1:end)));


figure('Name','eigenvalues');
subplot(4,1,1);
plot(largest_eig1(2:end),'color','black','lineWidth',3);
hold on;plot(ones(length(largest_eig1(2:end)),1),'color','red');
title('Largest eigenvalue: normal regime');
subplot(4,1,2);
plot(largest_eig2(2:end),'color','black','lineWidth',3);
hold on;
plot(ones(length(largest_eig2(2:end)),1),'color','red');
title('Largest eigenvalue: zlb regime','lineWidth',3);
subplot(4,1,3);
plot(average_eig(2:end),'color','black','lineWidth',3);
hold on;
plot(ones(length(average_eig(2:end)),1),'color','red');
subplot(4,1,4);
title('Largest eigenvalue: weighted average');
plot(pr_flag(2:end),'color','blac');title('projection facility activity');

figure('Name','regime probability');
area(1-pp_filtered);
hold on;
plot(S_filtered(:,12),'lineWidth',3);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_ar1_regimeProb','-dpdf'); 



 figure('Name','AR(1) Learning Coef-Intercept','units','normalized','outerposition',[0 0 1 1]);
plot(squeeze(learning_filtered(:,:,1)'),'color','black');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_ar1_learning_alphas','-dpdf');  




 figure('Name','AR(1) Learning Coef-Betas','units','normalized','outerposition',[0 0 1 1]);
plot(squeeze(learning_filtered(:,:,2)'),'color','black');
hold on;
plot(ones(T,1),'--','color','red');
xlim([1 T]);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'sw_ar1_learning_betas','-dpdf'); 




% output_gap=S_filtered(:,8)-param(10)*S_filtered(:,18);
% gap_growth=diff(output_gap(2:end));
% figure('Name','Output gap and growth rate');
% subplot(2,1,1);
% plot(output_gap);
% subplot(2,1,2);
% plot(gap_growth);


% figure('Name','inflation intercept');
% subplot(2,1,1);
% plot(S_filtered(:,10));
% subplot(2,1,2);
% plot(learning_filtered(10,:,1));






