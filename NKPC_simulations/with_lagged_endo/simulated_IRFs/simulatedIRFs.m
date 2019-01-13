clear;clc;close all;%tic
%------------------SIMULATION
seed=round(1000*rand);
param1=[0 0 0 0.01 3 1.5 0.5 0.9 0.9 0 0.7 0.3 0.3 ];
param2=[0 0 0 0.01 3 0   0   0.9 0.9 0 0.7 0.3 0.03 ];

numEndo=3;numExo=3;
N=200;zlb_length=50;irf_length=50;
numVar=5;
burn_in=round(N/10);
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
regime=ones(N,1);%regime parameter: 1 if not at ZLB, 0 if at ZLB;
regime(end-zlb_length+1:end)=zeros(zlb_length,1);
p_11=.99;p_22=.9; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

aa_tt=zeros(numEndo,1);%constant. coef for learning
cc_tt=zeros(numEndo,numEndo);%coef. on lagged endo variables
dd_tt=zeros(numEndo,numExo);%coef on shocks
rr_tt=1*eye(numEndo+numExo+1);%auxiliary learning matrix
expectations=zeros(numEndo,N);

d_REE= (eye(size(C1,1)^2)-kron(RHO',(ergodic_states(1)*...
A1_inv*C1+ergodic_states(2)*A2_inv*C2)))^(-1)*vec(A1_inv*ergodic_states(1)+A2_inv*ergodic_states(2));
d_REE=reshape(d_REE,[size(C1,1),size(C1,2)]);
 gain=0.01;
 
dd_tt=d_REE;
 dd(1,:,:)=dd_tt;
for tt=2:N
  
disp(tt);

EPS1(:,tt)=RHO*EPS1(:,tt-1)+errors1(:,tt);
EPS2(:,tt)=RHO*EPS2(:,tt-1)+errors2(:,tt);
% regime(tt)=findRegime(regime(tt-1),p_11,p_22);
EPS(:,tt)=EPS1(:,tt)*regime(tt)+EPS2(:,tt)*(1-regime(tt));


expectations(1:numEndo,tt)=...
    (aa_tt+cc_tt*aa_tt)+cc_tt^2*X(1:numEndo,tt-1)+ (cc_tt*dd_tt+dd_tt*RHO)*EPS(:,tt);

    
X(:,tt) =  regime(tt)*(A1_inv*(B1*X(:,tt-1)+C1*expectations(:,tt)+EPS1(:,tt)))+...;%state realized
       (1-regime(tt))*(A2_inv*(B2*X(:,tt-1)+C2*expectations(:,tt)+EPS2(:,tt)));

thetaOld=[aa_tt cc_tt dd_tt];
[theta rr_tt] =l_LS_version2(X(:,tt),[1;X(:,tt-1);EPS(:,tt)],thetaOld,rr_tt,gain);
aa_tt=theta(1,:)';cc_tt=theta(2:4,:)';dd_tt=theta(5:7,:)';
cc_tt=zeros(3,3);




gamma1_1_tilde=A1_inv*(B1+C1*cc_tt^2);gamma3_1_tilde=(A1_inv*C1)*(cc_tt*dd_tt+dd_tt*RHO+eye(3));
gamma1_2_tilde=A2_inv*(B2+C2*cc_tt^2);gamma3_2_tilde=(A2_inv*C2)*(cc_tt*dd_tt+dd_tt*RHO+eye(3));

gamma1_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO];
gamma3_1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(3,3);eye(3)];

gamma1_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[gamma1_2_tilde,zeros(numEndo,numExo);zeros(numExo,numEndo),RHO];
gamma3_2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1)*[zeros(3,3);eye(3)];

gamma1_avg=regime(tt)*gamma1_1+(1-regime(tt))*gamma1_2;

gamma3_avg=regime(tt)*gamma3_1+(1-regime(tt))*gamma3_2;

imp(tt,:,:,:)=impulse_response(gamma1_avg,gamma3_avg,irf_length);

end


imp=imp(:,:,1:2,1:2);
start_irf=145;
end_irf=155;

for mm=1:2
    for nn=1:2
        subplot(2,2,(mm-1)*2+nn)
        for kk=start_irf:end_irf
   plot(imp(kk,:,mm,nn));
   hold on;
        end
    end
end




% refLine1=linspace(1,irf_length,irf_length);
% refLine1=repmat(refLine1,[N 1]);
% 
% refLine2=linspace(1,irf_length,irf_length);
% refLine2=repmat(refLine2,[N 1])';
% 
% index=0;
% for ii=1:2
%     for jj=1:2
% 
%         index=index+1;
%          subplot(2,2,index);
%  plot3(refLine1',refLine2',imp(:,:,ii,jj)','color','red','lineWidth',0.5);
% 
% 
%     end
% end