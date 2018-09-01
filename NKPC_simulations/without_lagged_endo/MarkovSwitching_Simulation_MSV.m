%Simulation of 3-equation NKPC with Markov-switching, MSV-learning with
%least squares.
clear;clc;close all;%tic
addpath('c:\users\tolga\desktop\zlb_paper\optimization_routines');
%------------------SIMULATION
seed=round(1000*rand);%rng(99)
param1=[0 0 0 0.01 3 1.5 0.5 0.5 0.5 0 0.7 0.3 0.001 ];
param2=[0 0 0 0.01 3 0    0  0.5 0.5 0 0.7 0.3 0.001 ];
numEndo=3;numExo=3;N=5000;
numVar=5;

sigma_y1 = param1(end-2);
sigma_pinf1=param1(end-1);
sigma_r1=param1(end);

sigma_y2 = param2(end-2);
sigma_pinf2=param2(end-1);
sigma_r2=param2(end);


[A1 B1 C1 D1]= NKPC_sysmat(param1);
[A2 B2 C2 D2]= NKPC_sysmat(param2);
RHO=zeros(numExo,numExo);RHO(1:2,1:2)=B1(numEndo+1:end,numEndo+1:end);
A1=A1(1:numEndo,1:numEndo);B1=B1(1:numEndo,1:numEndo);C1=C1(1:numEndo,1:numEndo);


A2=A2(1:numEndo,1:numEndo);B2=B2(1:numEndo,1:numEndo);C2=C2(1:numEndo,1:numEndo);



A1_inv=A1^(-1);A2_inv=A2^(-1);

X=zeros(numEndo,N);X(:,1)=zeros(numEndo,1);
EPS1=zeros(numExo,N);EPS1(:,1)=zeros(numExo,1);
EPS2=zeros(numExo,N);EPS2(:,1)=zeros(numExo,1);
EPS=zeros(numExo,N);EPS(:,1)=zeros(numExo,1);
eps_y(:,1) = normrnd(0,sigma_y1,[N,1]);
eps_y(:,2) = normrnd(0,sigma_y2,[N,1]);
eps_pinf(:,1) = normrnd(0,sigma_pinf1,[N,1]);    
eps_pinf(:,2) = normrnd(0,sigma_pinf2,[N,1]); 
eps_r(:,1) = normrnd(0,sigma_r1,[N,1]);
eps_r(:,2) = normrnd(0,sigma_r2,[N,1]);
errors1=[eps_y(:,1) eps_pinf(:,1) eps_r(:,1)]' ; 
errors2=[eps_y(:,2) eps_pinf(:,2) eps_r(:,2)]' ; 
regime=nan(N,1);%regime parameter: 0 if not at ZLB, 1 if at ZLB;
regime(1)=0;
p_11=0.99;p_22=0.95; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

aa_tt=zeros(numEndo,1);%constant. coef for learning
cc_tt=zeros(numEndo,numEndo);%coef. on lagged endo variables
dd_tt=zeros(numEndo,3);%coef on shocks
rr_tt=10*eye(3);%auxiliary learning matrix
expectations=zeros(numEndo,N);

d_REE= (eye(size(C1,1)^2)-kron(RHO',(ergodic_states(1)*...
A1_inv*C1+ergodic_states(2)*A2_inv*C2)))^(-1)*vec(A1_inv*ergodic_states(1)+A2_inv*ergodic_states(2));
d_REE=reshape(d_REE,[size(C1,1),size(C1,2)]);

for tt=2:N
    gain=0.01;
disp(tt);

EPS1(:,tt)=RHO*EPS1(:,tt-1)+errors1(:,tt);
EPS2(:,tt)=RHO*EPS2(:,tt-1)+errors2(:,tt);
regime(tt)=findRegime(regime(tt-1),p_11,p_22);
EPS(:,tt)=EPS1(:,tt)*regime(tt)+EPS2(:,tt)*(1-regime(tt));


expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*EPS(:,tt);

    
X(:,tt) =  regime(tt)*(A1_inv*(B1*X(:,tt-1)+C1*expectations(:,tt)+EPS1(:,tt)))+...;%state realized
       (1-regime(tt))*(A2_inv*(B2*X(:,tt-1)+C2*expectations(:,tt)+EPS2(:,tt)));

%thetaOld=[aa_tt cc_tt dd_tt];
thetaOld=[aa_tt dd_tt(:,1:2)];
[theta rr_tt] =l_LS_version2(X(:,tt),[1;EPS(1:2,tt)],thetaOld,rr_tt,gain,[1 3 3]);
aa_tt=theta(1,:)';%cc_tt=theta(2,:)';
dd_tt(:,1:2)=theta(2:3,:)';

aa(:,tt)=aa_tt;
%cc(tt,:,:)=cc_tt;
dd(tt,:,:)=dd_tt;
learningCovariance(tt,:,:)=rr_tt;

end

figure('Name','Endogenous variables');
subplot(2,1,1);
plot(X(1,:),'lineWidth',3);
subplot(2,1,2);
plot(X(2,:),'lineWidth',3);

figure('Name','learning coefficients-intercepts');
subplot(2,1,1);
plot(aa(1,:),'lineWidth',3);
hold on;
plot(zeros(N,1),'lineWidth',3);
xlim([0 N]);
subplot(2,1,2);
plot(aa(2,:),'lineWidth',3);
hold on;
plot(zeros(N,1),'lineWidth',3);
xlim([0 N]);
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MS_simulation_alphas','-dpdf');

% figure('Name','learning-lagged variable coefficients. should converge to zero');
% subplot(2,2,1);
% plot(cc(:,1,1),'lineWidth',3);
% hold on;
% plot(zeros(N,1),'lineWidth',3);
% subplot(2,2,2);
% plot(cc(:,1,2),'lineWidth',3);
% hold on;
% plot(zeros(N,1),'lineWidth',3);
% subplot(2,2,3);
% plot(cc(:,2,1),'lineWidth',3);
% hold on;
% plot(zeros(N,1),'lineWidth',3);
% subplot(2,2,4);
% plot(cc(:,2,2),'lineWidth',3);
% hold on;
% plot(zeros(N,1),'lineWidth',3);

figure('Name','learning:-shock coefficients');

subplot(2,2,1);
plot(dd(:,1,1),'lineWidth',3);
hold on;
plot(d_REE(1,1)*ones(N,1),'lineWidth',3);
xlim([0 N]);
subplot(2,2,2);
plot(dd(:,1,2),'lineWidth',3);
hold on;
plot(d_REE(1,2)*ones(N,1),'lineWidth',3);
xlim([0 N]);
xlim([0 N]);
subplot(2,2,3);
plot(dd(:,2,1),'lineWidth',3);
hold on;
plot(d_REE(2,1)*ones(N,1),'lineWidth',3);
xlim([0 N]);
subplot(2,2,4);
plot(dd(:,2,2),'lineWidth',3);
hold on;
plot(d_REE(2,2)*ones(N,1),'lineWidth',3);

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'MS_simulation_shockCoef','-dpdf');

figure;
M1=size(rr_tt,1);
M2=size(rr_tt,2);



for jj=1:M1
    for ii=1:M2
        subplot(M1,M2,(jj-1)*M1+ii);
        plot(learningCovariance(:,jj,ii),'lineWidth',3);
    end
end

figure;
area(regime);
