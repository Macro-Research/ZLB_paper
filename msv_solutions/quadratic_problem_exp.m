clear;clc;close all;
rng(1);
m_states=2;
PHI=[1,0;0,1];
LAMBDA=0.7*[1,0;0,1];
THETA=rand(m_states,m_states);

Xi_mat=[LAMBDA,THETA;eye(m_states),zeros(m_states,m_states)];
Delta_mat=[PHI,zeros(m_states,m_states);zeros(m_states,m_states),eye(m_states)];


[Xi_eigvec,Xi_eigval]=eig(Xi_mat,Delta_mat);
xx=Xi_eigvec(m_states+1:end,:);
OMEGA=licols(xx,10e-5);
Xi_eigval=Xi_eigval(Xi_eigval~=0);
Xi_eigval=diag(Xi_eigval);
PP=OMEGA*Xi_eigval*OMEGA^(-1);


       CHECK_SOL=PHI*PP^2-LAMBDA*PP-THETA