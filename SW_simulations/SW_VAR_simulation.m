clear;clc;close all;
% m=round(1000*rand);rng(863)
N=5000;numVar=24;numShocks=7;numEndo=13;numExo=7;
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [4.4,4.4] ;    %phi
parameters(7,:)  =[1.02,1.02] ;      %sigma_c
parameters(8,:)  =  [0.8,0.8];   %lambda 
parameters(9,:)  =  [0.67,0.67] ; %xi_w
parameters(10,:) = [2.36,2.36];        %sigma_l 
parameters(11,:) =  [0.74,0.74];      %xi_p 
parameters(12,:) =  [0,0] ;     %iota_w
parameters(13,:) =  [0,0];         %iota_p
parameters(14,:) = [0.46,0.46] ;   %psi
parameters(15,:) = [1.5,1.5]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.43,1.43];   %r_pi
parameters(17,:)    =  [0.9,0.9]; %rho
parameters(18,:)    =   [0.08,0.08];%0.0746;    %r_y
parameters(19,:)    =  [0.05,0.05];     %r_dy

%SS related parameters
parameters(20,:)    = [0.57,0.57] ;    %pi_bar
parameters(21,:)    = [0.18,0.18];       %beta_const
parameters(22,:)    = [0,0];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.17,0.17] ;    %alpha

% %shock persistence
parameters(25,:) = [0.97,0.97];      %rho_a
parameters(26,:) =[0.35,035];  %rho_b
parameters(27,:) =[0.98,0.98];   %rho_g
parameters(28,:) =[0.54,0.54];   %rho_i
parameters(29,:) =  [0.35,0.35];   %rho_r
parameters(30,:) =  [0.16,0.16];  %rho_p
parameters(31,:) =[0.04,0.04];%rho_w 
parameters(32,:) = [0,0] ;    %mu_p 
parameters(33,:) =[0,0];    %mu_w
parameters(34,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[0.41,0.41];  %sigma_a
parameters(36,:)=  [0.55,0.55];  %sigma_b
parameters(37,:)= [0.4,0.4]; %sigma_g
parameters(38,:)=  [1.24,   1.24] ;   %sigma_i
parameters(39,:)=[0.1,0.1];   %sigma_r
parameters(40,:)= [0.2,0.2];  %sigma_p
parameters(41,:)=   [0.8,0.8];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=0;
gain=0.01;
p_11=1;p_22=0; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
% 
% [AA1, BB1, CC1, DD1, EE1]=SW_matrixConverter_VAR(parameters(:,1));
% [AA2, BB2, CC2, DD2, EE2]=SW_matrixConverter_VAR(parameters(:,2));


[AA1, BB1, CC1, DD1, EE1]=SW_sysmat_VAR_filter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2]=SW_sysmat_VAR_filter(parameters(:,2));
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

SIGMA1=diag(parameters(end-numShocks+1:end,1));
SIGMA2=diag(parameters(end-numShocks+1:end,2));

errors1=mvnrnd(zeros(numShocks,1),SIGMA1,N)';
errors2=mvnrnd(zeros(numShocks,1),SIGMA2,N)';


XX=nan(numVar,N);
XX(:,1) =zeros(numVar,1);
beta_tt=0.1*eye(numVar);
alpha_tt=zeros(numVar,1);
rr_tt=99*repmat(eye(2),[1 1 numVar]);

for tt=2:N
% gain=0.001;
    disp(tt)
      regime(tt)=findRegime(regime(tt-1),p_11,p_22);
    XX(:,tt) = (regime(tt))* AA1_inv * (( BB1+CC1*beta_tt^2)*XX(:,tt-1)+CC1*(alpha_tt+beta_tt*alpha_tt)+DD1*errors1(:,tt)+EE1*errors1(:,tt-1))+...
    (1-regime(tt))* AA2_inv * (( BB2+CC2*beta_tt^2)*XX(:,tt-1)+CC2*(alpha_tt+beta_tt*alpha_tt)+DD2*errors2(:,tt)+EE2*errors2(:,tt-1));
beta_old=beta_tt;
            for jj=1:numVar
         [alpha_tt(jj) beta_tt(jj,jj) rr_tt(:,:,jj),largest_eig(tt),pr_flag(tt)] =...
              msv_learning(XX(jj,tt),[1,XX(jj,tt-1)]',...
           alpha_tt(jj),beta_tt(jj,jj),rr_tt(:,:,jj),gain);
%  [alpha_tt(jj),beta_tt(jj,jj)]=l_SAC_DGL(XX(jj,1:tt),tt);
% [alpha_tt(jj) beta_tt(jj,jj) rr_tt(1,1,jj)] =l_SAC_CGL(XX(jj,tt),XX(jj,tt-1),alpha_tt(jj),beta_tt(jj,jj),rr_tt(1,1,jj),gain);
            end 
          
            try
largest_eig2(tt)=abs(eigs(AA1_inv*(BB1+CC1*beta_tt^2),1));
            catch
                largest_eig(tt)=1.01;
            end
if largest_eig2(tt)>1
    beta_tt=beta_old;
    pr_flag(tt)=1;
end


matrix_learning(tt,:,:)=[alpha_tt,diag(beta_tt)];

end


