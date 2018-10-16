clear;clc;close all;
rng(11);
%m=round(1000*rand);rng(863)
%variable order: mc zcap rk k1   q c inve y lab pinf w r kp 
N=2000;numVar=24;numShocks=7;numEndo=17;numExo=7;numBackward=7;numForward=7;
eig_crit=.99;
eig_mult=(1-eig_crit)*2/pi;
%burn_in=round(0.9*N);
burn_in=2;
backward_indices=[6 7 8 10 11 12 13];
forward_indices=[3 5 6 7 9 10 11];
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [6.48,6.48] ;    %phi
parameters(7,:)  =[2.06,2.06] ;      %sigma_c
parameters(8,:)  =  [0.68,0.68];   %lambda 
parameters(9,:)  =  [0.57,0.57] ; %xi_w
parameters(10,:) = [1.78,1.78];        %sigma_l 
parameters(11,:) =  [0.7,0.7];      %xi_p 
parameters(12,:) =  [0.58,0.58] ;     %iota_w
parameters(13,:) =  [0.43,0.43];         %iota_p
parameters(14,:) = [0.69,0.69] ;   %psi
parameters(15,:) = [1.21,1.2]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.5,0];   %r_pi
parameters(17,:)    =  [0.87,0]; %rho
parameters(18,:)    =   [0.09,0];%0.0746;    %r_y
parameters(19,:)    =  [0.15,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0.43,0.43] ;    %pi_bar
parameters(21,:)    = [0.18,0.18];       %beta_const
parameters(22,:)    = [0,0];     %l_bar
parameters(23,:)    = [0.34,0.34];         %gamma_bar
parameters(24,:)    = [0.14,0.14] ;    %alpha

% %shock persistence
parameters(25,:) = [0.98,0.98];      %rho_a
parameters(26,:) =[0.1,0.1];  %rho_b
parameters(27,:) =[0.98,0.98];   %rho_g
parameters(28,:) =[0.89,0.89];   %rho_i
parameters(29,:) =  [0.65,0.65];   %rho_r
parameters(30,:) =  [0.33,0.33];  %rho_p
parameters(31,:) =[0.11,0.11];%rho_w 
% parameters(32,:) = [0.25,0.25] ;    %mu_p 
% parameters(33,:) =[0.25,0.25];    %mu_w

parameters(32,:) =  [0.5,0.5]      ;  %rho_ga

%shock standard deviations
parameters(33,:)=[0.42,0.42];  %sigma_a
parameters(34,:)=  [0.19,0.19];  %sigma_b
parameters(35,:)= [0.38,0.38]; %sigma_g
parameters(36,:)=  [0.29,0.29] ;   %sigma_i
parameters(37,:)=[0.09,0.01];   %sigma_r
parameters(38,:)= [0.12,0.12];  %sigma_p
parameters(39,:)=   [0.42,0.42];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
gain=0.01;
p_11=0.97;p_22=0.86; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];


[AA1, BB1, CC1, DD1, EE1 ,RHO1 ,FF1, GG1 ,E1 ,F1] =SW_sysmat_MSV_filter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2, RHO2, FF2 ,GG2, E2, F2]=SW_sysmat_MSV_filter(parameters(:,2));

AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

SIGMA1=diag(parameters(end-numShocks+1:end,1))^2;
SIGMA2=diag(parameters(end-numShocks+1:end,2))^2;

errors1=mvnrnd(zeros(numShocks,1),SIGMA1,N)';
errors2=mvnrnd(zeros(numShocks,1),SIGMA2,N)';



XX=nan(numVar,N);
XX(:,1) =zeros(numVar,1);
%----
beta_tt=zeros(numEndo,numEndo);
%  beta_tt=0.1*eye(numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=zeros(numEndo,numShocks);
rr_tt=1*ones(numBackward+numExo+1,numBackward+numExo+1);
rr_tt=nearestSPD(rr_tt);
% ----
load('initial_beliefs.mat');
beta_tt(forward_indices,backward_indices)=sqrtm(beta_init);
%cc_tt(forward_indices,:)=cc_init;
%rr_tt=rr_init;

theta=zeros(numExo+1,numForward);
pr_flag=zeros(N,1);
largest_eig1=zeros(N,1);
largest_eig2=zeros(N,1);
%----------------------------------
gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
gamma2_1_tilde=(AA1\CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_1_tilde=(AA1\CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+(AA1\DD1);

gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);
gamma2_2_tilde=(AA2\CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_2_tilde=(AA2\CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+(AA2\DD2);

aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_1=aux_matrix1\([gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1]);
gamma2_1=aux_matrix1\([gamma2_1_tilde;zeros(numExo,1)]);
gamma3_1=aux_matrix1\([(AA1\EE1);FF1]);

aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_2=aux_matrix2\([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
gamma2_2=aux_matrix2\([gamma2_2_tilde;zeros(numExo,1)]);
gamma3_2=aux_matrix2\([(AA2\EE2);FF2]);
%----------------
index=0;
eig_max=nan(N,1);
for tt=2:N


      regime(tt)=findRegime(regime(tt-1),p_11,p_22);
      
    XX(:,tt) = regime(tt)*( gamma1_1*XX(:,tt-1)+gamma2_1+gamma3_1*errors1(:,tt))+...
    (1-regime(tt))*( gamma1_2*XX(:,tt-1)+gamma2_2+gamma3_2*errors2(:,tt)); 





%--------------------------------store old expectation terms
alpha_old=alpha_tt;
beta_old=beta_tt;
cc_old=cc_tt;
rr_old=rr_tt;

thetaOld=[alpha_tt(forward_indices)...
    beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)];
% 
% thetaOld=[alpha_tt(forward_indices) cc_tt(forward_indices,:)];
%-------------------------------

% [theta,rr_tt ,largest_eig1(tt),largest_eig2(tt),pr_flag(tt)] =...
%     msv_learning_withoutBetas(XX(forward_indices,tt),...
%     [1;XX(numEndo+1:end,tt)],...
%     thetaOld,rr_tt,gain,...
%     forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_crit);

[theta,rr_tt ,largest_eig1(tt),largest_eig2(tt),pr_flag(tt)] =...
    msv_learning3(XX(forward_indices,tt),...
    [1;XX(backward_indices,tt-1);XX(numEndo+1:end,tt)],...
    thetaOld,rr_old,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_crit);

eig_max(tt)=max(largest_eig1(tt),largest_eig2(tt));
%-------------------------------update expectations
%if eig_max(tt)<eig_crit %update alpha beta and cc if eigenvalue ok
if largest_eig1(tt)<eig_crit
alpha_tt(forward_indices)=theta(1,:)';
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
else pr_flag(tt)=1;
    alpha_tt(forward_indices)=theta(1,:)';
    cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
    beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
    
    [VV,DD]=eig(beta_tt);
    DD=diag(DD);
    ii=find(DD>eig_crit);
    DD(ii)=0.99;
    beta_tt=real(VV*diag(DD)*pinv(VV));
    %flag the projection facility if eigenvalue not ok
   rr_tt=rr_old;
end

%------------------------------update system matrices
if largest_eig2(tt)<eig_crit
    
gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
gamma2_1_tilde=(AA1\CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_1_tilde=(AA1\CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+(AA1\DD1);

gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);
gamma2_2_tilde=(AA2\CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_2_tilde=(AA2\CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+(AA2\DD2);

aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_1=aux_matrix1\([gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1]);
gamma2_1=aux_matrix1\([gamma2_1_tilde;zeros(numExo,1)]);
gamma3_1=aux_matrix1\([(AA1\EE1);FF1]);

aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_2=aux_matrix2\([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
gamma2_2=aux_matrix2\([gamma2_2_tilde;zeros(numExo,1)]);
gamma3_2=aux_matrix2\([(AA2\EE2);FF2]);   
else
    pr_flag(tt)=1;
    
gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
    [VV,DD]=eig(gamma1_1_tilde);
    DD=diag(DD);
    ii=find(DD>eig_crit);
    DD(ii)=0.99;

gamma1_1_tilde=real(VV*diag(DD)*pinv(VV));
gamma2_1_tilde=(AA1\CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_1_tilde=(AA1\CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+(AA1\DD1);

gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);
    [VV,DD]=eig(gamma1_2_tilde);
    DD=diag(DD);
    ii=find(DD>eig_crit);
    DD(ii)=0.99;

gamma1_2_tilde=real(VV*diag(DD)*pinv(VV));
gamma2_2_tilde=(AA2\CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_2_tilde=(AA2\CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+(AA2\DD2);

aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_1=aux_matrix1\([gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1]);
gamma2_1=aux_matrix1\([gamma2_1_tilde;zeros(numExo,1)]);
gamma3_1=aux_matrix1\([(AA1\EE1);FF1]);

aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_2=aux_matrix2\([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
gamma2_2=aux_matrix2\([gamma2_2_tilde;zeros(numExo,1)]);
gamma3_2=aux_matrix2\([(AA2\EE2);FF2]);

end
%------------------------------------------------------------------
 index=index+1;
 if index==500
 disp(['projection facility activity:',num2str(mean(pr_flag(1:tt)))]);
%   disp(['Largest eigenvalue:',num2str(max(abs(largest_eig2(tt)),(abs(largest_eig(tt)))))]);
 disp(['Largest eigenvalue:',num2str(abs(eig_max(tt)))]);
  disp(tt);
   index=0;
 end

learning_matrix(tt,:,:)=[alpha_tt(forward_indices)...
    beta_tt(forward_indices,backward_indices) cc_tt(forward_indices,:)]';

end


figure('Name','intercepts');
plot(squeeze(learning_matrix(burn_in:end,1,:)),'color','black');

% figure('Name','first-order autocorrelations');
% plot(squeeze(learning_matrix(burn_in:end,:,2)),'color','black');

figure('Name','lagged inflation');
plot(squeeze(learning_matrix(burn_in:end,5,:)),'color','black');

figure('Name','lagged consumption');
plot(squeeze(learning_matrix(burn_in:end,2,:)),'color','black');

figure('Name','lagged interest rate');
plot(squeeze(learning_matrix(burn_in:end,7,:)),'color','black');

figure('Name','monetary policy shocks');
plot(squeeze(learning_matrix(burn_in:end,13,:)),'color','black');

figure('Name','largest eigenvalue of the system and projection facility');
subplot(5,1,1);
plot(largest_eig2,'lineWidth',3);title('largest eigenvalue, implied alm');
hold on;
plot(ones(N,1),'--');
subplot(5,1,2);
plot(largest_eig1,'lineWidth',3);title('largest eigenvalue, plm');
hold on;
plot(ones(N,1),'--');
subplot(5,1,3);
plot(eig_max,'lineWidth',3);title('largest eigenvalue, greater of the two');
hold on;
plot(ones(N,1),'--');
subplot(5,1,4);
plot(pr_flag,'lineWidth',3);title('projection flag');
subplot(5,1,5);
area(regime);

figure('Name','Inflation');
plot(XX(10,:),'lineWidth',3);