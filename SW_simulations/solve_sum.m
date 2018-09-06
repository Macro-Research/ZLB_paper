clear;clc;close all;
%fixed parameteters
parameters(1,:)   = [0.025,0.025];     %delta
parameters(2,:)      = [0.18,0.18];    % G 
parameters(3,:)    =[1.5,1.5];       % phi_w 
parameters(4,:)   = [10,10];        %curv_p
parameters(5,:)   = [10,10];         %curv_w
 
%nominal and real frictions
parameters(6,:)  = [2.0395,2.0395] ;    %phi
parameters(7,:)  =[1.2510,1.2510] ;      %sigma_c
parameters(8,:)  =  [0.66,0.66];   %lambda 
parameters(9,:)  =  [0.67,0.67] ; %xi_w
parameters(10,:) = [2.45,2.45];        %sigma_l 
parameters(11,:) =  [0.63,0.63];      %xi_p 
parameters(12,:) =  [0,0] ;     %iota_w
parameters(13,:) =  [0,0];         %iota_p
parameters(14,:) = [0.54,0.54] ;   %psi
parameters(15,:) = [1.53,1.53]; %phi_p

%policy related parameters

parameters(16,:)    =   [1.4434,0];   %r_pi
parameters(17,:)    =  [0.82,0]; %rho
parameters(18,:)    =   [0.12,0];%0.0746;    %r_y
parameters(19,:)    =  [0.15,0];     %r_dy

%SS related parameters
parameters(20,:)    = [0.71,0.71] ;    %pi_bar
parameters(21,:)    = [0.20,0.20];       %beta_const
parameters(22,:)    = [-4.37,-4.37];     %l_bar
parameters(23,:)    = [0.39,0.39];         %gamma_bar
parameters(24,:)    = [0.17,0.17] ;    %alpha

% %shock persistence
parameters(25,:) = [0.84,0.84];      %rho_a
parameters(26,:) =[0.21,0.21];  %rho_b
parameters(27,:) =[0.92,0.92];   %rho_g
parameters(28,:) =[0.51,0.51];   %rho_i
parameters(29,:) =  [0.06,0.06];   %rho_r
parameters(30,:) =  [0.18,0.18];  %rho_p
parameters(31,:) =[0.062,0.062];%rho_w 
parameters(32,:) = [0,0] ;    %mu_p 
parameters(33,:) =[0,0];    %mu_w
parameters(34,:) =  [0,0]      ;  %rho_ga

%shock standard deviations
parameters(35,:)=[0.55,0.1];  %sigma_a
parameters(36,:)=  [0.74,0.1];  %sigma_b
parameters(37,:)= [0.67,0.1]; %sigma_g
parameters(38,:)=  [1.74,0.5] ;   %sigma_i
parameters(39,:)=[0.44,0.1];   %sigma_r
parameters(40,:)= [0.35,0.1];  %sigma_p
parameters(41,:)=   [0.53,0.1];   %sigma_w
p_11=0.99;p_22=0.99; 
Q=[p_11,1-p_11;1-p_22,p_22];
ergodic_states=[(1-p_22)/(2-p_11-p_22);(1-p_11)/(2-p_11-p_22)];

[AA1, BB1, CC1, DD1, EE1]=SW_matrixConverter(parameters(:,1));
[AA2, BB2, CC2, DD2, EE2]=SW_matrixConverter(parameters(:,2));
AA1_inv=AA1^(-1);AA2_inv=AA2^(-1);

A_tilde_inv = AA1^(-1)*ergodic_states(1)+AA2^(-1)*ergodic_states(2);
B_tilde= A_tilde_inv*(BB1*ergodic_states(1)+BB2*ergodic_states(2));
C_tilde=A_tilde_inv*(CC1*ergodic_states(1)+CC2*ergodic_states(2));
D_tilde=A_tilde_inv*(DD1*ergodic_states(1)+DD2*ergodic_states(2));
% RHO_tilde=ergodic_states(1)*rho1+ergodic_states(2)*rho2;

m_states=20;
PHI=C_tilde;
LAMBDA=eye(m_states);%ok
THETA=-B_tilde;%ok
MANUAL_ROOTS=1;
DISPLAY_IMMEDIATELY=1;
TOL=1e-6;
warnings='';

Xi_mat=[LAMBDA,THETA;eye(m_states),zeros(m_states,m_states)];
Delta_mat=[PHI,zeros(m_states,m_states);zeros(m_states,m_states),eye(m_states)];


% v=linspace(1,2*m_states,2*m_states);k=m_states;
% numComb=nchoosek(length(v),k);
% allComb=nchoosek(v,k);
% 
% PP=zeros(20,20);
index=0;
% for jj=1:numComb
%     
%     disp(jj);
%     Xi_manual=allComb(jj,:);
Xi_manual=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];
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
% end

 bb_RPE=PP_all(:,:,1);

% vec_dd=(eye(m_states*2)-kron(eye(m_states-1),C_tilde*bb_RPE)-kron(RHO_tilde',C_tilde))^(-1)*vec(D_tilde); 
% dd_RPE=reshape(vec_dd,[3 2]); 
% aa=zeros(m_states,1);
% 
% 
%       disp('DYNARE OUTPUT:');
%       disp('bb_dynare:');
%    %disp(bb_dynare);
%    disp('dd_dynare:');
%    %disp(dd_dynare);
%    disp('----------------');
%    disp('TOOLBOX OUTPUT:');
%    disp('bb:');
%    disp(bb_RPE);
%    disp('dd:');
%    disp(dd_RPE);
