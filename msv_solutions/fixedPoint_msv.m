
clear;
clc;

parameters=[0.03; 0.82; 1.19; 0.2; 2.99;1.36 ; 0.48; 0.42; 0.31; 0.85; 0.73; 0.29; 0.29];

varCovar=[parameters(end-2)^2,0,0;0,parameters(end-1)^2,0;0,0,parameters(end)^2];
varCovar_vec=reshape(varCovar,[length(varCovar)^2,1]);
numVar=5;numEndo=3;numExo=2;
grid1=50;


[Atotal, Btotal, Ctotal, Dtotal]=NKPC_sysmat(parameters);
gamma1=Atotal^(-1)*Btotal;
gamma2=Atotal^(-1)*Ctotal;
gamma3=Atotal^(-1)*Dtotal;


% Convergence Analysis
N=500;

beta=zeros(numVar,numVar,N);
beta(:,:,1) = rand   *diag( ones(numVar,1));
for i=1:N
    dampingFactor=0.25;
    betaAux = beta(:,:,i);
    M=gamma1+gamma2*betaAux^2;
    
vec0=(eye(numVar^2)-kron(M,M))^(-1)*kron(gamma3,gamma3)*varCovar_vec;
vec1=(kron(eye(numVar),gamma1)+kron(eye(numVar),gamma2*betaAux^2))*vec0;

vv0=reshape(vec0,[5 5]);
vv1=reshape(vec1,[5 5]);
beta(1:3,1:3,i+1)=vv0(1:3,1:3)^(-1)*vv1(1:3,1:3);

%   for jj=1:numEndo
%         for kk=1:numEndo
%             if vec0((jj-1)*numVar+kk)==0 || vec0((jj-1)*numVar+kk)<0
%                 beta(jj,kk,i+1)=0;
%             else
%         beta(jj,kk,i+1)=beta(jj,kk,i)+dampingFactor*...
%             (vec1((jj-1)*numVar+kk)/vec0((jj-1)*numVar+kk)-beta(jj,kk,i));
%             end
%         end
%   end



  

    

end

disp(beta(:,:,end))