clear;clc;close all;
%rng(1);

%calibration
sigma=2; kappa=0.05; beta=0.99; rho=0; phi_y=.5;phi_pinf=1.5;eta_var=0.01;
i_star=(1/beta)-1;

%matrices
A= [ 1+phi_y/sigma, phi_pinf/sigma;-kappa 1];
B = [0;0];
C=[1,1/sigma;0,beta];
D=[1;0];

gamma1=A^(-1)*B;gamma2=A^(-1)*C;gamma3=A^(-1)*D;

a=gamma3;

%simulation
N=50000;
eta=normrnd(0,eta_var,[N 1]);
% 
eps_y=nan(N,1);
eps_y(1)=eta(1);
for jj=2:N
    eps_y(jj)=rho*eps_y(jj-1)+eta(jj);
end

z=eps_y*a';



z_learning=nan(N,2);
% z_learning(1,:)=zeros(2,1);
expectations=nan(N,2);
expectations(1,:)=rand(2,1);

for jj=1:N

    gain=.05;
    if jj>1
    expectations(jj,:)=gain*z_learning(jj-1,:)+(1-gain)*expectations(jj-1,:);
    end
   z_learning(jj,:)=(gamma2*expectations(jj,:)')+gamma3*eps_y(jj);
   
 
end

figure;
subplot(2,1,1);
plot(z(:,1),'lineWidth',3,'color','red');
hold on;
plot(z_learning(:,1),'lineWidth',3,'color','blue','lineStyle','--');
legend('msv','learning');
subplot(2,1,2);
plot(z(:,2),'lineWidth',3,'color','red');
hold on;
plot(z_learning(:,2),'lineWidth',3,'color','blue','lineStyle','--');
legend('msv','learning');


figure;
subplot(2,1,1);
plot(gamma1(1)*ones(N,1),'lineWidth',3);
hold on;
plot(expectations(:,1),'lineWidth',3);
legend('msv','learning');

subplot(2,1,2);
plot(gamma1(2)*ones(N,1),'lineWidth',3);
hold on;
plot(expectations(:,2),'lineWidth',3);
legend('msv','learning');


   
