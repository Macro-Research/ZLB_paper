clear;clc;close all;
rng(111);
%m=round(1000*rand);rng(863)
%variable order: mc zcap rk k1   q c inve y lab pinf w r kp 
N=2000;numVar=24;numShocks=7;numEndo=17;numExo=7;numBackward=7;numForward=7;
%burn_in=round(0.9*N);
burn_in=1;
backward_indices=[6 7 8 10 11 12 13];
forward_indices=[3 5 6 7 9 10 11];
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [6.5,6.5] ;    %phi
parameters(7,:)  =[1.8,1.8] ;      %sigma_c
parameters(8,:)  =  [0.69,0.69];   %lambda 
parameters(9,:)  =  [0.75,0.75] ; %xi_w
parameters(10,:) = [1.78,1.78];        %sigma_l 
parameters(11,:) =  [0.75,0.75];      %xi_p 
parameters(12,:) =  [0.58,0.58] ;     %iota_w
parameters(13,:) =  [0.43,0.43];         %iota_p
parameters(14,:) = [0.7,0.7] ;   %psi
parameters(15,:) = [1.21,1.21]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.5,1.5];   %r_pi
parameters(17,:)    =  [0.85,0.85]; %rho
parameters(18,:)    =   [0.2,0.2];%0.0746;    %r_y
parameters(19,:)    =  [0.2,0.2];     %r_dy

%SS related parameters
parameters(20,:)    = [0.5,0.5] ;    %pi_bar
parameters(21,:)    = [0.13,0.13];       %beta_const
parameters(22,:)    = [0,0];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.17,0.17] ;    %alpha

% %shock persistence
parameters(25,:) = [0.8,0.8];      %rho_a
parameters(26,:) =[0.2,0.2];  %rho_b
parameters(27,:) =[0.8,0.8];   %rho_g
parameters(28,:) =[0.5,0.5];   %rho_i
parameters(29,:) =  [0.3,0.3];   %rho_r
parameters(30,:) =  [0.1,0.1];  %rho_p
parameters(31,:) =[0.1,0.1];%rho_w 
% parameters(32,:) = [0.25,0.25] ;    %mu_p 
% parameters(33,:) =[0.25,0.25];    %mu_w

parameters(32,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(33,:)=[0.5,0.5];  %sigma_a
parameters(34,:)=  [0.2,0.2];  %sigma_b
parameters(35,:)= [0.5,0.5]; %sigma_g
parameters(36,:)=  [0.3,0.3] ;   %sigma_i
parameters(37,:)=[0.1,0.1];   %sigma_r
parameters(38,:)= [0.2,0.2];  %sigma_p
parameters(39,:)=   [0.4,0.4];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
gain=0.01;
p_11=1;p_22=0; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];


[AA1, BB1, CC1, DD1, EE1 ,RHO1 ,FF1, GG1 ,E1 ,F1] =SW_sysmat_MSV_filter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2, RHO2, FF2 ,GG2, E2, F2]=SW_sysmat_MSV_filter(parameters(:,2));

AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

SIGMA1=diag(parameters(end-numShocks+1:end,1));
SIGMA2=diag(parameters(end-numShocks+1:end,2));

errors1=mvnrnd(zeros(numShocks,1),SIGMA1,N)';
errors2=mvnrnd(zeros(numShocks,1),SIGMA2,N)';



XX=nan(numVar,N);
XX(:,1) =zeros(numVar,1);
beta_tt=zeros(numEndo,numEndo);
alpha_tt=zeros(numEndo,1);
cc_tt=zeros(numEndo,numShocks);
rr_tt=10*eye(numBackward+numExo+1);
% rr_tt=99*repmat(eye(2),[1 1 numVar]);
theta=zeros(numBackward+numExo+1,numForward);
pr_flag=nan(N,1);
largest_eig=zeros(N,1);
largest_eig2=zeros(N,1);

update_matrices;
index=0;
for tt=2:N
% gain=1/tt;

      regime(tt)=findRegime(regime(tt-1),p_11,p_22);
      
    XX(:,tt) = regime(tt)*( gamma1_1*XX(:,tt-1)+gamma2_1+gamma3_1*errors1(:,tt))+...
    (1-regime(tt))*( gamma1_2*XX(:,tt-1)+gamma2_2+gamma3_2*errors2(:,tt)); 

if tt<1000
 update_beliefs;
end
% 
  update_matrices;
%  check_eigenvalue;
%  revert_if_explosive;
%  update_matrices;


 index=index+1;
 if index==500
 disp(['projection facility activity:',num2str(pr_flag(tt))]);
  disp(['Largest eigenvalue:',num2str(max(abs(largest_eig2(tt)),(abs(largest_eig(tt)))))]);
  disp(tt);
   index=0;
 end

learning_matrix(tt,:,:)=theta;
% largest_eig3(tt)=eigs(gamma1_1,1);
% learning_matrix(tt,:,:)=theta';
%largest_eig2(tt)=eigs(gamma1_1,1);
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
subplot(4,1,1);
plot(largest_eig2,'lineWidth',3);title('largest eigenvalue, implied alm');
hold on;
plot(ones(N,1),'--');
subplot(4,1,2);
plot(largest_eig,'lineWidth',3);title('largest eigenvalue, plm');
hold on;
plot(ones(N,1),'--');
subplot(4,1,3);
plot(pr_flag,'lineWidth',3);
subplot(4,1,4);
area(regime);

figure('Name','Inflation');
plot(XX(10,:),'lineWidth',3);