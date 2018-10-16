clear;clc;close all;
%rng(15);
m=round(1000*rand);rng(592)
%variable order: mc zcap rk k1   q c inve y lab pinf w r kp 
N=10000;numVar=24;numShocks=7;numEndo=17;numExo=7;numBackward=7;numForward=7;
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
parameters(6,:)  = [6.45,6.45] ;    %phi
parameters(7,:)  =[1.49,1.49] ;      %sigma_c
parameters(8,:)  =  [0.73,0.73];   %lambda 
parameters(9,:)  =  [0.79,0.79] ; %xi_w
parameters(10,:) = [1.54,1.54];        %sigma_l 
parameters(11,:) =  [0.76,0.76];      %xi_p 
parameters(12,:) =  [0.56,0.56] ;     %iota_w
parameters(13,:) =  [0.16,0.16];         %iota_p
parameters(14,:) = [0.78,0.78] ;   %psi
parameters(15,:) = [1.51,1.51]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.35,0];   %r_pi
parameters(17,:)    =  [0.86,0]; %rho
parameters(18,:)    =   [0.11,0];%0.0746;    %r_y
parameters(19,:)    =  [0.03,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0.6,0.6] ;    %pi_bar
parameters(21,:)    = [0.18,0.18];       %beta_const
parameters(22,:)    = [1.53,1.53];     %l_bar
parameters(23,:)    = [0.5,0.5];         %gamma_bar
parameters(24,:)    = [0.19,0.19] ;    %alpha

% %shock persistence
parameters(25,:) = [0.96,0.96];      %rho_a
parameters(26,:) =[0.58,0.58];  %rho_b
parameters(27,:) =[0.95,0.95];   %rho_g
parameters(28,:) =[0.80,0.80];   %rho_i
parameters(29,:) =  [0.59,0.59];   %rho_r
parameters(30,:) =  [0.92,0.92];  %rho_p
parameters(31,:) =[0.29,0.29];%rho_w 
% parameters(32,:) = [0.25,0.25] ;    %mu_p 
% parameters(33,:) =[0.25,0.25];    %mu_w

parameters(32,:) =  [0.45,0.45]      ;  %rho_ga

%shock standard deviations
parameters(33,:)=[0.37,0.37];  %sigma_a
parameters(34,:)=  [0.16,0.16];  %sigma_b
parameters(35,:)= [0.39,0.39]; %sigma_g
parameters(36,:)=  [0.28,0.28] ;   %sigma_i
parameters(37,:)=[0.09,0.01];   %sigma_r
parameters(38,:)= [0.07,0.07];  %sigma_p
parameters(39,:)=   [0.40,0.40];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
gain=0.0173;
p_11=0.9741;p_22=0.8682; 
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
alpha_tt=zeros(numEndo,1);
cc_tt=zeros(numEndo,numShocks);
rr_tt=0*ones(numBackward+numExo+1,numBackward+numExo+1);
rr_tt=triu(rr_tt);
rr_tt=(rr_tt+rr_tt')/2;
% % ----
 load('initial_beliefs.mat');
beta_tt=beta_init;
cc_tt=cc_init;
rr_tt=rr_init;
rr_tt=triu(rr_tt);
rr_tt=(rr_tt+rr_tt')/2;

pr_flag=zeros(N,1);
largest_eig1=zeros(N,1);
largest_eig2=zeros(N,1);
%----------------------------------
update_matrices
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


    

%-------------------------------
[theta rr_tt] =l_LS_version2(XX(forward_indices,tt),[1;XX(backward_indices,tt-1);XX(numEndo+1:end,tt)],...
    thetaOld,rr_tt,gain);
%------------------------------
alpha_tt(forward_indices)=theta(1,:)';
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
%------------------------------------------------------------------
update_matrices;
largest_eig(tt)=abs(eigs(gamma1_1,1));

if largest_eig(tt)>1
    pr_flag(tt)=1;
    alpha_tt=alpha_old;
    beta_tt=beta_old;
    cc_tt=cc_old;
    rr_tt=rr_old;
    update_matrices;
end


%------------------------------------------------------------------
 index=index+1;
 if index==500
 disp(['projection facility activity:',num2str(mean(pr_flag(1:tt)))]);
%   disp(['Largest eigenvalue:',num2str(max(abs(largest_eig2(tt)),(abs(largest_eig(tt)))))]);
 disp(['Largest eigenvalue:',num2str(abs(largest_eig(tt)))]);
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

figure('Name','Inflation');
plot(XX(10,:),'lineWidth',3);

figure('Name','largest eigenvalue of the system and projection facility');
plot(largest_eig(2:end),'lineWidth',3,'color','black');
hold on;
plot(ones(N-1,1),'--','color','red');
%subplot(5,1,1);
%plot(largest_eig2,'lineWidth',3);title('largest eigenvalue, implied alm');
%hold on;
%plot(ones(N,1),'--');
%subplot(5,1,2);
%plot(largest_eig1,'lineWidth',3);title('largest eigenvalue, plm');
%hold on;
%plot(ones(N,1),'--');
%subplot(5,1,3);
% plot(eig_max,'lineWidth',3);title('largest eigenvalue, greater of the two');
% hold on;
% plot(ones(N,1),'--');
% subplot(5,1,4);
% plot(pr_flag,'lineWidth',3);title('projection flag');
% subplot(5,1,5);
% area(regime);

