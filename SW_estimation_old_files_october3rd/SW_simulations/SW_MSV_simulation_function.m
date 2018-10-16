 function[pr_flag_mean,theta_final] = SW_MSV_simulation_function()

N=2000;numVar=24;numShocks=7;numEndo=17;numExo=7;numBackward=7;numForward=7;
eig_thr=.98;
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
parameters(6,:)  = [4,4] ;    %phi
parameters(7,:)  =[1.5,1.5] ;      %sigma_c
parameters(8,:)  =  [0.7,0.7];   %lambda 
parameters(9,:)  =  [0.75,0.75] ; %xi_w
parameters(10,:) = [2,2];        %sigma_l 
parameters(11,:) =  [0.75,0.75];      %xi_p 
parameters(12,:) =  [0,0] ;     %iota_w
parameters(13,:) =  [0,0];         %iota_p
parameters(14,:) = [0.5,0.5] ;   %psi
parameters(15,:) = [1.25,1.25]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.5,0];   %r_pi
parameters(17,:)    =  [0.75,0]; %rho
parameters(18,:)    =   [0.125,0];%0.0746;    %r_y
parameters(19,:)    =  [0.125,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0.5,0.5] ;    %pi_bar
parameters(21,:)    = [0.13,0.13];       %beta_const
parameters(22,:)    = [0,0];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.3,0.3] ;    %alpha

% %shock persistence
parameters(25,:) = [0.9,0.9];      %rho_a
parameters(26,:) =[0.5,0.5];  %rho_b
parameters(27,:) =[0.9,0.9];   %rho_g
parameters(28,:) =[0.5,0.5];   %rho_i
parameters(29,:) =  [0.5,0];   %rho_r
parameters(30,:) =  [0,0];  %rho_p
parameters(31,:) =[0,0];%rho_w 
% parameters(32,:) = [0.25,0.25] ;    %mu_p 
% parameters(33,:) =[0.25,0.25];    %mu_w

parameters(32,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(33,:)=[0.4,0.4];  %sigma_a
parameters(34,:)=  [0.18,0.18];  %sigma_b
parameters(35,:)= [0.36,0.36]; %sigma_g
parameters(36,:)=  [0.28,0.28] ;   %sigma_i
parameters(37,:)=[0.08,0.01];   %sigma_r
parameters(38,:)= [0.11,0.11];  %sigma_p
parameters(39,:)=   [0.4,0.4];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
gain=0.01;
p_11=0.98;p_22=0.85; 
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
rr_tt=10*eye(numBackward+numExo+1);

check=1.1;
while check>1
beta_tt(forward_indices,backward_indices)=0.5*rand(numForward,numBackward);
gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
check=abs(eigs(gamma1_1_tilde,1));
end
alpha_tt(forward_indices)=rand(numForward,1);
cc_tt(forward_indices,:)=rand(numForward,size(cc_tt,2));

% ----
load('initial_beliefs.mat');


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
gamma1_2=aux_matrix2^(-1)*([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
gamma2_2=aux_matrix2^(-1)*([gamma2_2_tilde;zeros(numExo,1)]);
gamma3_2=aux_matrix2^(-1)*([(AA2\EE2);FF2]);
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

[theta,rr_tt ,largest_eig1(tt),largest_eig2(tt),pr_flag(tt)] =...
    msv_learning3(XX(forward_indices,tt),...
    [1;XX(backward_indices,tt-1);XX(numEndo+1:end,tt)],...
    thetaOld,rr_tt,gain,...
    forward_indices,backward_indices,numEndo,AA1,BB1,CC1,AA2,BB2,CC2,eig_thr);

eig_max(tt)=max(largest_eig1(tt),largest_eig2(tt));
%-------------------------------update expectations
if eig_max(tt)<eig_thr %update alpha beta and cc if eigenvalue ok
alpha_tt(forward_indices)=theta(1,:)';
beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';
cc_tt(forward_indices,:)=theta(numBackward+2:end,:)';
else pr_flag(tt)=1;%flag the projection facility if eigenvalue not ok
  rr_tt=rr_old;
end


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



end

pr_flag_mean=mean(pr_flag);
theta_final=theta';


 end
