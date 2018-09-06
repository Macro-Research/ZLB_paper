clear;clc;close all;

% a_11=0.5; a_12 = -0.5; a_21=0.1 ;a_22 = 0.9;
% A=[a_11,a_12;a_21,a_22];
% disp(abs(eigs(A)))
N=10000;
numVar=3;


EPS=mvnrnd(zeros(numVar,1),eye(numVar),N)';

X=zeros(numVar,N);

A= [0.3709  , -0.0324 ,  -0.0614;
   -1.1020 ,   0.1661  , -0.1518;
   -0.0633   , 0.1341,    0.8438];
for jj=2:N
    
    X(:,jj)=A*X(:,jj-1)+EPS(:,jj);
    
end

figure;
subplot(2,1,1);
plot(X(1,:));
subplot(2,1,2);
plot(X(2,:));

acf1=autocorr(X(1,:));
acf2=autocorr(X(2,:));
acf3=autocorr(X(3,:));

gamma0=reshape( ((eye(length(A)^2)-kron(A,A))^(-1))*vec(eye(length(A))),size(A));
gamma1=A*gamma0;

rho=gamma1./gamma0
acf1(2)
acf2(2)
acf3(2)