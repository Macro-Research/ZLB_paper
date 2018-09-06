clear;clc;close all;
%rng(1);

%calibration
sigma=2; kappa=0.02; beta=0.99; rho=0.9; phi_y=0;phi_pinf=0;eta_var=0.1;
i_star=(1/beta)-1;

%matrices
A= [ 1+phi_y/sigma, phi_pinf/sigma;-kappa 1];
B = [0;0];
C=[1,1/sigma;0,beta];
D=[1;0];

gamma1=A^(-1)*B;gamma2=A^(-1)*C;gamma3=A^(-1)*D;

a=(eye(size(gamma2))-gamma2*rho)^(-1)*gamma3;

%simulation
N=500;
eta=normrnd(0,eta_var,[N 1]);

eps_y=nan(N,1);
eps_y(1)=eta(1);
for jj=2:N
    eps_y(jj)=rho*eps_y(jj-1)+eta(jj);
end

Z1=nan(N,2);Z2=nan(N,2);

% Z= [eps_y,eps_y].*repmat(a,[1 N])';
Z1(:,1)=a(1)*eps_y;aux=(gamma2*a*rho+gamma3);Z2(:,1)=aux(1)*eps_y;
Z1(:,2)=a(2)*eps_y;Z2(:,2)=aux(2)*eps_y;
% figure;
% plot(Z1(:,1));
% hold on;
% plot(Z2(:,2));
   
%learning
Z_learning=nan(N,2);
Z_learning(1,:)=100*rand;
a_hat=nan(2,N);
a_hat(:,1)=100*rand(2,1);
r_aux1=0;r_aux2=0;
for jj=1:N

   Z_learning(jj,:)=(gamma2*a_hat(:,jj)*rho*eps_y(jj))+(gamma3*eps_y(jj));
   
    X=eps_y(1:jj);
    Y1=Z_learning(1:jj,1);
    Y2=Z_learning(1:jj,2);
    a_hat(1,jj+1)=Y1(end)/X(end);
    a_hat(2,jj+1)=Y2(end)/X(end);
%     a_hat(1,jj+1)=msv_learningUpdate(Y1,X,a_hat(1,jj),r_aux1);
%     a_hat(2,jj+1)=msv_learningUpdate(Y2,X,a_hat(2,jj),r_aux2);
    
%     a_hat(1,jj+1)= (X'*X)^(-1)*(X'*Y1);
%     a_hat(2,jj+1)= (X'*X)^(-1)*(X'*Y2);
%     a_hat(1,jj+1)=sum(X.*Y1)/(sum(X.^2));
%     a_hat(2,jj+1)=sum(X.*Y2)/(sum(X.^2));



end
 
figure;
plot(Z_learning(:,1));
hold on;
plot(Z1(:,1),'*');
legend('learning','msv');

figure;
subplot(2,1,1);
plot(a(1)*ones(N,1));
hold on;
plot(a_hat(1,:));
subplot(2,1,2);
plot(a(2)*ones(N,1));
hold on;
plot(a_hat(2,:));
legend('true','learning');
