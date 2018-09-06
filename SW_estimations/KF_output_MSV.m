clear;clc;close all;
numVar=24;numShocks=7;numEndo=13;numExo=7;
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [4.4,4.4] ;    %phi
parameters(7,:)  =[1.5,1.5] ;      %sigma_c
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
parameters(22,:)    = [1.81,1.81];     %l_bar
parameters(23,:)    = [0.4,0.4];         %gamma_bar
parameters(24,:)    = [0.13,0.13] ;    %alpha

% %shock persistence
parameters(25,:) = [0.97,0.97];      %rho_a
parameters(26,:) =[0.35,035];  %rho_b
parameters(27,:) =[0.98,0.98];   %rho_g
parameters(28,:) =[0.54,0.54];   %rho_i
parameters(29,:) =  [0,0];   %rho_r
parameters(30,:) =  [0.13,0.13];  %rho_p
parameters(31,:) =[0.04,0.04];%rho_w 
parameters(32,:) = [0,0] ;    %mu_p 
parameters(33,:) =[0,0];    %mu_w
parameters(34,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[0.41,0.41];  %sigma_a
parameters(36,:)=  [0.55,0.55];  %sigma_b
parameters(37,:)= [0.4,0.4]; %sigma_g
parameters(38,:)=  [1.3,1.3] ;   %sigma_i
parameters(39,:)=[100,0.15];   %sigma_r
parameters(40,:)= [0.2,0.2];  %sigma_p
parameters(41,:)=   [0.85,0.85];   %sigma_w
gain=0.01;
p_11=0.99;p_22=0.99; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

[ AA1 BB1 CC1 DD1 EE1 RHO1 FF1 GG1 E1 F1] = SW_sysmat_MSV_filter( parameters(:,1) );
[ AA2 BB2 CC1 DD1 EE1 RHO1 FF1 GG1 E1 F1] = SW_sysmat_MSV_filter( parameters(:,1) );

%E2(6)=0.01;
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

Sigma1=diag(parameters(end-numShocks+1:end,1));
Sigma2=diag(parameters(end-numShocks+1:end,2));
load('full_dataset.mat');
%load('raf_dataset.mat');
%dataset=[dy dc dw dinve pinfobs robs labobs];%-->when using sysmat_VAR_filter
dataset=[dy dc dinve dw pinfobs robs labobs];%-->when using sysmat_old
first_obs=44;
last_obs=size(dataset,1);
burn_in=1;
dataset=dataset(first_obs:last_obs,:);
T=size(dataset,1);numObs=7;l=7;
alpha1=0*ones(numVar,1);
beta1=0.1*eye(numVar);
beta1(3,3)=0.99;beta1(5,5)=0.8 ;beta1(6,6)=0.99;beta1(7,7)=0.99;beta1(9,9)=0.98;beta1(10,10)=0.51;beta1(11,11)=0.98;
% 
rr=eye(numVar);


gamma1_1=AA1_inv*(BB1+CC1*beta1^2);gamma2_1=AA1_inv*CC1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=AA1_inv*DD1;
gamma1_2=AA2_inv*(BB2+CC2*beta1^2);gamma2_2=AA2_inv*CC2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=AA1_inv*DD2;

H1=0*diag(var(dataset(1:150,:)));H2=0*diag(var(dataset(150:end,:)));



S_fore11=zeros(numVar,1);S_fore12=zeros(numVar,1);S_fore21=zeros(numVar,1);S_fore22=zeros(numVar,1);
P_fore11=zeros(numVar,numVar);P_fore12=zeros(numVar,numVar);P_fore21=zeros(numVar,numVar);P_fore22=zeros(numVar,numVar);
v11=zeros(numObs,1);v12=zeros(numObs,1);v21=zeros(numObs,1);v22=zeros(numObs,1);
Fe11=zeros(numObs,numObs);Fe12=zeros(numObs,numObs);Fe21=zeros(numObs,numObs);Fe22=zeros(numObs,numObs);
S_upd11=zeros(numVar,1);S_upd12=zeros(numVar,1);S_upd21=zeros(numVar,1);S_upd22=zeros(numVar,1);
P_upd11=zeros(numVar,numVar);P_upd12=zeros(numVar,numVar);P_upd21=zeros(numVar,numVar);P_upd22=zeros(numVar,numVar);
ml11=0;ml12=0;ml21=0;ml22=0;
likl=zeros(T,1);
%pp_fore11=0;pp_fore12=0;pp_fore21=0;pp_fore22=0;
pp_fore11=0;pp_fore12=0;pp_fore21=0;pp_fore22=1;
pp_upd11=0;pp_upd12=0;pp_upd21=0;pp_upd22=1;
%pp_collapse1=regime(1);pp_collapse2=1-regime(1);
pp_collapse1=0.1;pp_collapse2=0.9;
S_collapse1=zeros(numVar,1);S_collapse2=zeros(numVar,1);
P_collapse1=eye(numVar);P_collapse2=eye(numVar);
q_11=Q(1,1);q_12=Q(1,2);q_21=Q(2,1);q_22=Q(2,2);
S_filtered=zeros(T,numVar);
%pp_filtered=zeros(T,1);
pp_filtered=ones(T,1);


%----------------------------------------------------------------------------

for tt=2:T
    
    %when using sac_learning algorithms
%gamma1_1=A1_inv*(B1+C1*beta1^2);gamma2_1=A1_inv*C1*(eye(numVar)-beta1^2)*alpha1;gamma3_1=A1_inv*D1;
%gamma1_2=A2_inv*(B2+C2*beta1^2);gamma2_2=A2_inv*C2*(eye(numVar)-beta1^2)*alpha1;gamma3_2=A1_inv*D2;

   %when using msv_learning algorithms. Only intercepts are different
gamma1_1=AA1_inv*(BB1+CC1*beta1^2);gamma2_1=AA1_inv*CC1*(eye(numVar)+beta1)*alpha1;gamma3_1=AA1_inv*DD1;
gamma1_2=AA2_inv*(BB2+CC2*beta1^2);gamma2_2=AA2_inv*CC2*(eye(numVar)+beta1)*alpha1;gamma3_2=AA2_inv*DD2;

    x_tt=dataset(tt,:)';

%kalman block
S_fore11=gamma1_1*S_collapse1+gamma2_1;%
S_fore12=gamma1_2*S_collapse1+gamma2_2;%
S_fore21=gamma1_1*S_collapse2+gamma2_1;%
S_fore22=gamma1_2*S_collapse2+gamma2_2;%
%
P_fore11=gamma1_1*P_collapse1*gamma1_1'+gamma3_1*Sigma1*gamma3_1';
P_fore12=gamma1_2*P_collapse1*gamma1_2'+gamma3_2*Sigma2*gamma3_2';  
P_fore21=gamma1_1*P_collapse2*gamma1_1'+gamma3_1*Sigma1*gamma3_1';
P_fore22=gamma1_2*P_collapse2*gamma1_2'+gamma3_2*Sigma2*gamma3_2';
%    
v11=x_tt-E1-F1*S_fore11;
v12=x_tt-E2-F2*S_fore12;
v21=x_tt-E1-F1*S_fore21;
v22=x_tt-E2-F2*S_fore22;
%
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
if pp_collapse1>10e-5;
    S_collapse1=(pp_upd11*S_upd11+pp_upd21*S_upd21)/pp_collapse1;
P_collapse1=(pp_upd11*(P_upd11+(S_collapse1-S_upd11)*(S_collapse1-S_upd11)')+...
    pp_upd21*(P_upd21+(S_collapse1-S_upd21)*(S_collapse1-S_upd21)'))/pp_collapse1;
else S_collapse1=zeros(numVar,1);
    P_collapse1=1*eye(numVar);
end

if pp_collapse2>10e-5
P_collapse2=(pp_upd12*(P_upd12+(S_collapse2-S_upd12)*(S_collapse2-S_upd12)')+...
    pp_upd22*(P_upd22+(S_collapse2-S_upd22)*(S_collapse2-S_upd22)'))/pp_collapse2;
S_collapse2=(pp_upd12*S_upd12+pp_upd22*S_upd22)/pp_collapse2;
else S_collapse2=zeros(numVar,1);
    P_collapse2=1*eye(numVar);
end

%
S_filtered(tt,:)=pp_collapse1*S_collapse1+pp_collapse2*S_collapse2;
pp_filtered(tt)=pp_collapse1;


end


likl=-sum(log(likl(burn_in+1:end)));



