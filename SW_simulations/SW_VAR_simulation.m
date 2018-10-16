clear;clc;close all;
 %rng(1);
 %m=round(1000*rand);rng(863)
N=2000;numVar=24;numShocks=7;numEndo=13;numExo=7;
eig_crit=.99;
eig_mult=(1-eig_crit)*2/pi;

backward_indices=[6 7 8 10 11 12 13];
forward_indices=[3 5 6 7 9 10 11];
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [6,6] ;    %phi
parameters(7,:)  =[1.35,1.35] ;      %sigma_c
parameters(8,:)  =  [0.3,0.3];   %lambda 
parameters(9,:)  =  [0.9,0.9] ; %xi_w
parameters(10,:) = [2.14,2.14];        %sigma_l 
parameters(11,:) =  [0.9,0.9];      %xi_p 
parameters(12,:) =  [0.1,0.1] ;     %iota_w
parameters(13,:) =  [0.1,0.1];         %iota_p
parameters(14,:) = [0.66,0.66] ;   %psi
parameters(15,:) = [1.31,1.31]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.42,0];   %r_pi
parameters(17,:)    =  [0.91,0]; %rho
parameters(18,:)    =   [0.17,0];%0.0746;    %r_y
parameters(19,:)    =  [0.11,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0.68,0.68] ;    %pi_bar
parameters(21,:)    = [0.2,0.2];       %beta_const
parameters(22,:)    = [3.03,3.03];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.3,0.3] ;    %alpha

% %shock persistence
parameters(25,:) = [0.99,0.99];      %rho_a
parameters(26,:) =[0.4,0.4];  %rho_b
parameters(27,:) =[0.99,0.99];   %rho_g
parameters(28,:) =[0.6,0.6];   %rho_i
parameters(29,:) =  [0.3,0];   %rho_r
parameters(30,:) =  [0.03,0.03];  %rho_p
parameters(31,:) =[0.07,0.07];%rho_w 
parameters(32,:) = [0,0] ;    %mu_p 
parameters(33,:) =[0,0];    %mu_w
parameters(34,:) =  [0.39,0.39]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[0.4,0.4];  %sigma_a
parameters(36,:)=  [0.8,0.8];  %sigma_b
parameters(37,:)= [0.38,0.38]; %sigma_g
parameters(38,:)=  [1.15,   1.15] ;   %sigma_i
parameters(39,:)=[0.09,0.01];   %sigma_r
parameters(40,:)= [0.18,0.18];  %sigma_p
parameters(41,:)=   [0.86,0.86];   %sigma_w
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=1;
gain=0.01;
p_11=0.99;p_22=0.5; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];
% 
% [AA1, BB1, CC1, DD1, EE1]=SW_matrixConverter_VAR(parameters(:,1));
% [AA2, BB2, CC2, DD2, EE2]=SW_matrixConverter_VAR(parameters(:,2));
% % 
% 
[AA1, BB1, CC1, DD1, EE1]=SW_sysmat_VAR_filter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2]=SW_sysmat_VAR_filter(parameters(:,2));
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

SIGMA1=diag(parameters(end-numShocks+1:end,1))^2;
SIGMA2=diag(parameters(end-numShocks+1:end,2))^2;

errors1=mvnrnd(zeros(numShocks,1),SIGMA1,N)';
errors2=mvnrnd(zeros(numShocks,1),SIGMA2,N)';


XX=rand(numVar,N);
beta_tt=0.5*diag(ones(numVar,1));
alpha_tt=zeros(numVar,1);
rr_tt=1*repmat(eye(2),[1 1 numVar]);
load('ar1_initial_beliefs.mat');
%rr_tt=rr_init;
%beta_tt=beta_init;



largest_eig=zeros(N,1);largest_eig2=zeros(N,1);pr_flag=zeros(N,1);
matrix_learning=zeros(N,numVar,2);
pr_flag=zeros(N,1);

   index=0;
   
   var_update_matrices;
for tt=2:N
matrix_learning(tt-1,:,:)=[alpha_tt,diag(beta_tt)];
   % gain=1/tt;
   
   index=index+1;
   if index==500
   disp(['ACTIVITY OF PROJECTION FACILITY:',num2str(100*mean(pr_flag(1:tt))),'%']);
    disp(tt)
    index=0;
   end
%draw regime
      regime(tt)=findRegime(regime(tt-1),p_11,p_22);
%iterate forward
    XX(:,tt) = (regime(tt))*( gamma1_1*XX(:,tt-1)+gamma1_2+gamma1_3*errors1(:,tt))+...
    (1-regime(tt))* (gamma2_1*XX(:,tt-1)+gamma2_2+gamma2_3*errors2(:,tt));
%store old expectation terms
beta_old=beta_tt;
rr_old=rr_tt;
alpha_old=alpha_tt;
%update expectations
% 
            for jj=forward_indices
%                 thetaOld=[alpha_old(jj),beta_old(jj,jj)];
% [theta,rr_tt(:,:,jj)] =...
% l_LS_version2(XX(jj,tt),[1;XX(jj,tt-1)],thetaOld,rr_old(:,:,jj),gain);
% alpha_tt(jj)=theta(1);
% beta_tt(jj,jj)=theta(2);

%  [alpha_tt(jj) beta_tt(jj,jj) rr_tt(:,:,jj)] =...
%      l_SAC_CGL(XX(jj,tt),XX(jj,tt-1),alpha_tt(jj),beta_tt(jj,jj),rr_tt(:,:,jj),gain);
 [alpha_tt(jj) beta_tt(jj,jj) rr_tt(:,:,jj),pr_flag(tt)] =...
           msv_learning(XX(jj,tt)',[1,XX(jj,tt-1)]',...
        alpha_tt(jj),beta_tt(jj,jj),rr_tt(:,:,jj),gain);
 
%    
            end 

         

            
var_update_matrices;
          
                        try
%largest_eig2(tt)=max(abs(eigs(gamma1_1,1)),abs(eigs(gamma2_1,1)));
largest_eig(tt)=(abs(eigs(gamma1_1,1)));
% largest_eig2(tt)=ergodic_states(1)*(abs(eigs(gamma1_1,1)))+...
%     ergodic_states(2)* (abs(eigs(gamma2_1,1)));
                        catch
                largest_eig(tt)=1.01;
%                 largest_eig2(tt)=1.01;
                        end
            
                        
% if largest_eig2(tt)>1 || largest_eig(tt)>1
if  largest_eig(tt)>1
    alpha_tt=alpha_old;
    beta_tt=beta_old;
    rr_tt=rr_old;
%rr_tt=rr_init;

    pr_flag(tt)=1;

var_update_matrices;
gamma1_1=reduce_eigenvalue(gamma1_1,eig_crit,eig_mult);


end


eig_admitted(tt)=abs(eigs(gamma1_1,1));


rr_inf(:,:,tt)=rr_tt(:,:,10);
end
abs(eigs(AA1_inv*(BB1+CC1)))
disp(['ACTIVITY OF PROJECTION FACILITY:',num2str(100*mean(pr_flag)),'%']);

figure('Name','intercepts');
plot(matrix_learning(:,forward_indices,1),'color','black');
figure('Name','persistence');
plot(matrix_learning(:,forward_indices,2),'color','black');
figure('Name','projection facility');
plot(pr_flag);
figure('Name','Inflation');
plot(XX(10,:));

figure('Name','eigenvalues');
subplot(4,1,1);
plot(pr_flag(2:end),'lineWidth',3);title('activity of projection facility');
subplot(4,1,2);
plot(largest_eig(2:end),'lineWidth',3);title('largest eigenvalue, normal regime');
hold on;
plot(ones(N-1,1),'--');


%---------------------------------------
% subplot(4,1,3);
% plot(largest_eig2(2:end),'lineWidth',3);
% hold on;
% plot(ones(N-1,1),'--');title('largest eigenvalue, weighted average');
% subplot(4,1,4);
% plot(eig_admitted(2:end),'lineWidth',3);
% hold on;
% plot(ones(N,1,1),'--');title('admitted eigenvalue');


