clear;clc;close all;
%Benchmark Fisher equation (no regime switching), msv learning. 
%Can be used for both CGL and DGL cases.
N=10000;
alpha=1.5;
rho=0.9;
iota_p=0;
eta_sigma=0.1;uu_sigma=0.1;
eta=normrnd(0,eta_sigma,[N 1]);
uu=normrnd(0,uu_sigma,[N 1]);
r=zeros(N,1);
pinf=zeros(N,1);


b_1=(alpha-sqrt(alpha^2-4*iota_p))/2
b_2=(alpha+sqrt(alpha^2-4*iota_p))/2

a_1=2/(alpha-2*rho+sqrt(alpha^2-4*iota_p))
a_2=2/(alpha-2*rho-sqrt(alpha^2-4*iota_p))

Estab_b1=2*b_1/alpha
Estab_b2=2*b_2/alpha
Estab_a1=(rho+b_1)/alpha
Estab_a2=(rho+b_2)/alpha

aa=zeros(N,1);
bb=zeros(N,1);
rr=eye(2);
for jj=2:N
    gain=0.01;
    r(jj)=rho*r(jj-1)+eta(jj);

pinf(jj)=(1/alpha)*(aa(jj-1)*rho+bb(jj-1)*aa(jj-1)+1)*r(jj)+...
         (1/alpha)*(bb(jj-1)^2+iota_p)*pinf(jj-1)+...
         uu(jj);

[aa(jj) bb(jj) rr] =...
    l_LS(pinf(jj),[r(jj) pinf(jj-1)]',aa(jj-1),bb(jj-1),rr,gain);

end

figure('Name','learning coef.','units','normalized','outerposition',[0 0 1 1]);
subplot(2,1,1);
plot(aa,'lineWidth',3);
hold on;
plot(a_1*ones(N,1),'lineWidth',5);
title('coefficient on r_t');
subplot(2,1,2);
plot(bb,'lineWidth',3);
hold on;
plot(b_1*ones(N,1),'lineWidth',5);
title('coefficient on \pi_t');
