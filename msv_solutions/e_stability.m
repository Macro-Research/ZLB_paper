clear;clc;close all;
%rng(1);

%calibration
sigma=2; kappa=0.1; beta=0.99; rho=0.99; phi_y=0.125;phi_pinf=1.5;eta_var=0.5;
i_star=(1/beta)-1;

A= [ 1+phi_y/sigma, phi_pinf/sigma;-kappa 1];
B = [0;0];
C=[1,1/sigma;0,beta];
D=[1;0];

gamma1=A^(-1)*B;gamma2=A^(-1)*C;gamma3=A^(-1)*D;

a=(eye(size(gamma2))-gamma2*rho)^(-1)*gamma3;

T_map=rho*gamma2*a+gamma3-a;