clear;clc;close all;
%dynare NKPC_REE_Simulation noclearall nolog;
%dd_dynare= [oo_.dr.ghu(4:5,1:2)' oo_.dr.ghu(1,1:2)']';
%bb_dynare=[oo_.dr.ghx(4:5) oo_.dr.ghx(1)]';
%clearvars -except bb_dynare dd_dynare;
%parameter indices-------
%y_bar=param(1); pi_bar=param(2);  r_bar=param(3);  
%gamma=param(4);  tau=param(5);  phi_pinf=param(6);  phi_y=param(7); 
%rho_y=param(8);  rho_pinf=param(9);  rho_r=param(10);  
%eps_y=param(11);  eps_pinf=param(12);  eps_r=param(13); 
parameters=[0 0 0 0.03 3 1.5 0.5 0.9 0.9 0 0.7 0.3 0.3 0 0.03  0 1];    
param(:,1)=parameters(1:13);
param(:,2)=parameters(1:13);
param(3,2)=parameters(14);

param(13,2)=parameters(15);

param(6,2)=0;param(7,2)=0;param(10,2)=0;
q_11=1-parameters(16);q_22=1-parameters(17); 

varCovar=[parameters(end-2)^2,0,0;0,parameters(end-1)^2,0;0,0,parameters(end)^2];
varCovar_vec=reshape(varCovar,[length(varCovar)^2,1]);
numVar=5;

[AA1 BB1 CC1 DD1 EE1 FF1 rho1 E1 F1 G1]=NKPC_sysmat_MSV(param(:,1));
[AA2 BB2 CC2 DD2 EE2 FF2 rho2 E2 F2 G2]=NKPC_sysmat_MSV(param(:,2));

Sigma1=diag([param(end-2,1)^2;param(end-1,1)^2;param(end,1)^2]);
Sigma2=diag([param(end-2,2)^2;param(end-1,2)^2;param(end,2)^2]);
ergodic_states=[(1-q_22)/(2-q_11-q_22);(1-q_11)/(2-q_11-q_22)];

A_tilde_inv = AA1^(-1)*ergodic_states(1)+AA2^(-1)*ergodic_states(2);
B_tilde= A_tilde_inv*(BB1*ergodic_states(1)+BB2*ergodic_states(2));
C_tilde=A_tilde_inv*(CC1*ergodic_states(1)+CC2*ergodic_states(2));
D_tilde=A_tilde_inv*(DD1*ergodic_states(1)+DD2*ergodic_states(2));
RHO_tilde=ergodic_states(1)*rho1+ergodic_states(2)*rho2;

m_states=3;
PHI=C_tilde;
LAMBDA=eye(m_states);%ok
THETA=-B_tilde;%ok
MANUAL_ROOTS=1;
DISPLAY_IMMEDIATELY=1;
TOL=1e-6;
warnings='';

Xi_mat=[LAMBDA,THETA;eye(m_states),zeros(m_states,m_states)];
Delta_mat=[PHI,zeros(m_states,m_states);zeros(m_states,m_states),eye(m_states)];


v=[1 2 3 4 5 6];k=3;
numComb=nchoosek(length(v),k);
allComb=nchoosek(v,k);

PP=zeros(5,5);
index=0;
for jj=1:numComb
    
    disp(jj);
    Xi_manual=allComb(jj,:);

 [Xi_eigvec,Xi_eigval] = eig(Xi_mat,Delta_mat);
     if rank(Xi_eigvec)<m_states,
        message = ['SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.         '
                   '         Try to run your program again with DO_QZ = 1.                 '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
     else
       [Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
       Xi_sortvec = Xi_eigvec(1:2*m_states,Xi_sortindex);
       Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
       Xi_select = 1 : m_states;
       if imag(Xi_sortval(m_states))~=0,
         if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) ) < TOL),
         % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE OCCURS MORE THAN ONCE!!
         % IF YOU HAVE THAT PROBLEM, PLEASE TRY MANUAL ROOT SELECTION.  
           drop_index = 1;
           while (abs(imag(Xi_sortval(drop_index)))>TOL) & (drop_index < m_states),
             drop_index = drop_index + 1;
           end;
           if drop_index >= m_states,
             message = ['SOLVE.M: You are in trouble. You have complex eigenvalues, and I cannot'
                        '   find a real eigenvalue to drop to only have conjugate-complex pairs.'
                        '   Put differently: your PP matrix will contain complex numbers. Sorry!'
                        '   Try increasing the dimension of your state space. You may then get  '
                        '   sunspots, too.                                                      '];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           else
             message = ['SOLVE.M: I will drop the lowest real eigenvalue to get real PP.        '
                        '         I hope that is ok. You may have sunspots.                     ']; 
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
             Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
           end; % if drop_index >= m_states,
         end; % if (abs( Xi_sortval(m_states) - ...
       end; % if imag(Xi_sortval(m_states))~=0,
       if MANUAL_ROOTS,
         message = ['SOLVE.M: You have chosen to select roots manually.  I am crossing my   '
                    '         fingers that you are doing it correctly.  In particular,      '
                    '         you should have defined Xi_manual.  Type help solve           '
                    '         and inspect SOLVE.M to get further information on how to do it'];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
         if exist('Xi_manual'),
            Xi_select = Xi_manual;
         else
            message = ['SOLVE.M: You have not defined Xi_manual.  Either define it or turn off '
                       '         the manual roots selection procedure with                     '
                       '         MANUAL_ROOTS = 0                                              '
                       '         Right now, I better let your calculations crash - sorry!      '
                       '         If you get results, they are based on previous calculations.  '];
            disp(message);
            warnings = [warnings;message];
         end; % if exist('Xi_manual'),
       else
         if max(Xi_select) < 2*m_states,
           if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
             message = ['SOLVE.M: You may be in trouble. There are stable roots NOT used for PP.'
                        '         I have used the smallest roots: I hope that is ok.            '  
                        '         If not, try manually selecting your favourite roots.          '
                        '         For manual root selection, take a look at the file solve.m    '
                        '         Watch out for sunspot solutions.                              '
                        '         Better yet: move the time index of some endogenous variables  '
                        '         back by one and turn them into (predetermined) state variables'];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           end; % if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
         end; % if max(Xi_select) < 2*m_states,
       end; % if MANUAL_ROOTS,
       if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
         message = ['SOLVE.M: You may be in trouble.  There are unstable roots used for PP. '
                    '         Keep your fingers crossed or change your model.               '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if max(abs(Xi_sortval(Xi_select))) ... 
       if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
         message = ['SOLVE.M: Your matrix PP contains a unit root. You probably do not have '
                    '         a unique steady state, do you?  Should not be a problem, but  '
                    '         you do not have convergence back to steady state after a shock'
                    '         and you should better not trust long simulations.             '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if abs( max(abs(Xi_sortval(Xi_select))) - 1 ... 
       Lambda_mat = diag(Xi_sortval(Xi_select));
    
       
       Omega_mat  = [Xi_sortvec((m_states+1):(2*m_states),Xi_select)];
       
    

       if rank(Omega_mat)<m_states,
         message = 'SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.          ';
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       else
         PP = Omega_mat*Lambda_mat/Omega_mat;
         PP_imag = imag(PP);
         PP = real(PP);
         if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
           message = ['SOLVE.M: PP is complex.  I proceed with the real part only.            '  
                      '         Hope that is ok, but you are probably really in trouble!!     '
                      '         You should better check everything carefully and be           '
                      '         distrustful of all results which follow now.                  '];
           if DISPLAY_IMMEDIATELY, disp(message); end;
          % warnings = [warnings;message];
         end; % if sum(sum(abs(PP_imag)))
      end; % if rank(Omega_mat)<m_states,
      % End of calculating the PP matrix.  Now comes the rest.
 
    end; % if rank(Xi_eigvec)<m_states,



       CHECK_SOL=PHI*PP^2-LAMBDA*PP-THETA;
       
       if abs(eigs(PP,1))<1
     index=index+1;
      PP_all(:,:,index)=PP;
       end
end

 bb_RPE=PP_all(:,:,1);

vec_dd=(eye(m_states*2)-kron(eye(m_states-1),C_tilde*bb_RPE)-kron(RHO_tilde',C_tilde))^(-1)*vec(D_tilde); 
dd_RPE=reshape(vec_dd,[3 2]); 
aa=zeros(m_states,1);


      disp('DYNARE OUTPUT:');
      disp('bb_dynare:');
   %disp(bb_dynare);
   disp('dd_dynare:');
   %disp(dd_dynare);
   disp('----------------');
   disp('TOOLBOX OUTPUT:');
   disp('bb:');
   disp(bb_RPE);
   disp('dd:');
   disp(dd_RPE);
