clear;clc;close all;
%m=round(1000*rand);rng(863)
N=10000;numVar=20;numShocks=7;numEndo=13;numExo=7;numBackward=6;numForward=7;
backward_indices=[6 7 8 10 12 13];
forward_indices=[3 4 5 6 7 9 10 11];
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [4.61,4.61] ;    %phi
parameters(7,:)  =[1.29,1.29] ;      %sigma_c
parameters(8,:)  =  [0.5,0.5];   %lambda 
parameters(9,:)  =  [0.5,0.5] ; %xi_w
parameters(10,:) = [1.85,1.85];        %sigma_l 
parameters(11,:) =  [0.5,0.5];      %xi_p 
parameters(12,:) =  [0,0] ;     %iota_w
parameters(13,:) =  [0,0];         %iota_p
parameters(14,:) = [0.74,0.74] ;   %psi
parameters(15,:) = [1.51,1.51]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.83,0];   %r_pi
parameters(17,:)    =  [0.5,0]; %rho
parameters(18,:)    =   [0.08,0];%0.0746;    %r_y
parameters(19,:)    =  [0.23,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0,0] ;    %pi_bar
parameters(21,:)    = [0.13,0.13];       %beta_const
parameters(22,:)    = [0,0];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.17,0.17] ;    %alpha

% %shock persistence
parameters(25,:) = [0.9,0.9];      %rho_a
parameters(26,:) =[0.2,0.2];  %rho_b
parameters(27,:) =[0.9,0.9];   %rho_g
parameters(28,:) =[0.2,0.2];   %rho_i
parameters(29,:) =  [0,0];   %rho_r
parameters(30,:) =  [0,0];  %rho_p
parameters(31,:) =[0,0];%rho_w 
parameters(32,:) = [0,0] ;    %mu_p 
parameters(33,:) =[0,0];    %mu_w
parameters(34,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[0.45,0.45];  %sigma_a
parameters(36,:)=  [1,1];  %sigma_b
parameters(37,:)= [0.5,0.5]; %sigma_g
parameters(38,:)=  [0.35,0.35] ;   %sigma_i
parameters(39,:)=[0.3,0.01];   %sigma_r
parameters(40,:)= [0.15,0.15];  %sigma_p
parameters(41,:)=   [0.4,0.4];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=0;
gain=0.01;
p_11=0.98;p_22=0.85; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

% [AA1, BB1, CC1, DD1, EE1 RHO1 FF1 GG1]=SW_matrixConverter_MSV_old(parameters(:,1));
% [AA2, BB2, CC2, DD2, EE2 RHO2 FF2 GG2]=SW_matrixConverter_MSV_old(parameters(:,2));

[AA1, BB1, CC1, DD1, EE1 RHO1 FF1 GG1]=SW_sysmat_MSV(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2 RHO2 FF2 GG2]=SW_sysmat_MSV(parameters(:,2));

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
rr_tt=100*eye(numBackward+numExo+1);

gamma1_1_tilde=AA1_inv*(BB1+CC1*beta_tt^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_1_tilde=(AA1_inv*CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta_tt^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_2_tilde=(AA2_inv*CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2_inv*DD2;



for tt=2:N

gamma1_1_tilde=AA1_inv*(BB1+CC1*beta_tt^2);gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_1_tilde=(AA1_inv*CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+AA1_inv*DD1;
gamma1_2_tilde=AA2_inv*(BB2+CC2*beta_tt^2);gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta_tt)*alpha_tt;gamma3_2_tilde=(AA2_inv*CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2_inv*DD2;


gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1];
gamma2_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_1_tilde;zeros(numExo,1)];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA1_inv*EE1;FF1];
gamma4_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(numEndo,numExo);GG1];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2];
gamma2_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma2_2_tilde;zeros(numExo,1)];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[AA2_inv*EE2;FF2];
gamma4_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(numEndo,numExo);GG2];
   
    
    
    

    disp(tt)
      regime(tt)=findRegime(regime(tt-1),p_11,p_22);
    XX(:,tt) = regime(tt)*( gamma1_1*XX(:,tt-1)+gamma2_1+gamma3_1*errors1(:,tt)+gamma4_1*errors1(:,tt-1))+...
    (1-regime(tt))*( gamma1_2*XX(:,tt-1)+gamma2_2+gamma3_2*errors2(:,tt)+gamma4_2*errors2(:,tt-1)); 

alpha_old=alpha_tt;beta_old=beta_tt;cc_old=cc_tt;    
thetaOld=[alpha_tt beta_tt(:,backward_indices) cc_tt];
[theta,rr_tt ,largestEig(tt),pr_flag(tt)] =msv_learning2(XX(1:numEndo,tt),[1;XX(backward_indices,tt-1);XX(numEndo+1:end,tt)],thetaOld,rr_tt,gain,numBackward,backward_indices);

alpha_tt=theta(1,:)';
beta_tt(:,backward_indices)=theta(2:numBackward+1,:)';
cc_tt=theta(numBackward+2:end,:)';

            try
largest_eig2(tt)=abs(eigs(AA1_inv*(BB1+CC1*beta_tt^2),1));
            catch
                largest_eig(tt)=1.01;
            end
            
    if largest_eig2(tt)>1
       alpha_tt=alpha_old;
    beta_tt=beta_old;
    cc_tt=cc_old;
    pr_flag(tt)=1;
    end

learning_matrix(tt,:,:)=theta;

largestEig2(tt)=eigs(gamma1_1,1);
end


