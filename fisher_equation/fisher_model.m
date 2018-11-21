clear;clc;%close all;
%Benchmark Fisher equation (no regime switching), msv learning. 
%Can be used for both CGL and DGL cases.
N=10000;
phi_pinf=.5;
rho=0.9;
eta_sigma=0.1;
eta=normrnd(0,eta_sigma,[N 1]);
r=zeros(N,1);
pinf=zeros(N,1);
rr=eye(2);
alpha=rand;
beta=500;
learning=nan(N,2);
beta_ree=(1/(phi_pinf-rho));
for jj=2:N
    gain=0.01;
    r(jj)=rho*r(jj-1)+eta(jj);
     learning(jj-1,:)=[alpha,beta];
   [alpha beta rr] =l_LS(pinf(jj-1),[1 r(jj-1)]',alpha,beta,rr,gain);
   alpha=0;
    pinf(jj)=(1/phi_pinf)*((alpha+beta*rho*r(jj))+r(jj));
    
   
end

figure;
subplot(2,1,1);
plot(learning(:,1),'lineWidth',3);
hold on;
plot(ones(N,1)*0,'--');
title('alpha');
subplot(2,1,2);
plot(learning(:,2),'lineWidth',3);
hold on;
plot(ones(N,1)*beta_ree,'--');
title('beta');
