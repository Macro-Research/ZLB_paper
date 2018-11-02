%------------------------------------------------------------------------
aux_matrix1=[AA1, -(CC1*(beta_tt*cc_tt+cc_tt*RHO1)+DD1);zeros(numExo,numEndo),eye(numExo,numExo)]^(-1);
gamma1_1_tilde=[BB1+CC1*beta_tt^2,zeros(numEndo,numExo);zeros(numExo,numEndo),RHO1];
gamma2_1_tilde=[CC1*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
gamma3_1_tilde=[(AA1_inv*EE1);FF1];

aux_matrix2=[AA2, -(CC2*(beta_tt*cc_tt+cc_tt*RHO2)+DD2);zeros(numExo,numEndo),eye(numExo,numExo)]^(-1);
gamma1_2_tilde=[BB2+CC2*beta_tt^2,zeros(numEndo,numExo);zeros(numExo,numEndo),RHO2];
gamma2_2_tilde=[CC2*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
gamma3_2_tilde=[(AA2_inv*EE2);FF2];

gamma1_1=aux_matrix1*gamma1_1_tilde;
gamma2_1=aux_matrix1*gamma2_1_tilde;
gamma3_1=aux_matrix1*gamma3_1_tilde;

gamma1_2=aux_matrix2*gamma1_2_tilde;
gamma2_2=aux_matrix2*gamma2_2_tilde;
gamma3_2=aux_matrix2*gamma3_2_tilde;
%----------------------------------------------------------------------
% aux_matrix1=[AA1 zeros(numEndo,numExo);zeros(numExo,numEndo),eye(numExo)]^(-1);
% gamma1_1_tilde=[BB1+CC1*beta_tt^2,DD1*RHO1;zeros(numExo,numEndo),RHO1];
% gamma2_1_tilde=[CC1*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
% gamma3_1_tilde=[DD1*FF1;FF1];
% 
% aux_matrix2=[AA2 zeros(numEndo,numExo);zeros(numExo,numEndo),eye(numExo)]^(-1);
% gamma1_2_tilde=[BB2+CC2*beta_tt^2,DD2*RHO2;zeros(numExo,numEndo),RHO2];
% gamma2_2_tilde=[CC2*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
% gamma3_2_tilde=[DD2*FF2;FF2];
% 
% gamma1_1=aux_matrix1*gamma1_1_tilde;
% gamma2_1=aux_matrix1*gamma2_1_tilde;
% gamma3_1=aux_matrix1*gamma3_1_tilde;
% 
% gamma1_2=aux_matrix2*gamma1_2_tilde;
% gamma2_2=aux_matrix2*gamma2_2_tilde;
% gamma3_2=aux_matrix2*gamma3_2_tilde;
%%--------------------------------------------------------------------
% 
%  gamma1_1=[AA1_inv*(BB1+CC1*beta_tt^2),...
%      AA1_inv*CC1*(beta_tt*cc_tt*RHO1+cc_tt*RHO1^2)+AA1_inv*DD1*RHO1;
%      zeros(numExo,numEndo),RHO1];
%  gamma2_1=[AA1_inv*CC1*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
%  gamma3_1=[AA1_inv*CC1*(beta_tt*cc_tt*FF1+cc_tt*RHO1*FF1)+AA1_inv*DD1*FF1;FF1];
% 
%  gamma1_2=[AA2_inv*(BB2+CC2*beta_tt^2),...
%      AA2_inv*CC2*(beta_tt*cc_tt*RHO2+cc_tt*RHO2^2)+AA2_inv*DD2*RHO2;
%      zeros(numExo,numEndo),RHO2];
%  gamma2_2=[AA2_inv*CC2*(alpha_tt+beta_tt*alpha_tt);zeros(numExo,1)];
%  gamma3_2=[AA2_inv*CC2*(beta_tt*cc_tt*FF2+cc_tt*RHO2*FF2)+AA2_inv*DD2*FF2;FF2];



% gamma1_1_tilde=AA1_inv*(BB1+CC1*beta_tt^2);
% gamma2_1_tilde=(AA1_inv*CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
% gamma3_1_tilde=(AA1_inv*CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+(AA1_inv*DD1);
% % 
% gamma1_2_tilde=AA2_inv*(BB2+CC2*beta_tt^2);
% gamma2_2_tilde=(AA2_inv*CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
% gamma3_2_tilde=(AA2_inv*CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2_inv*DD2;
% 
% aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1);
% gamma1_1=aux_matrix1*([gamma1_1_tilde,zeros(numEndo,numExo);zeros(numExo,numEndo),RHO1]);
% gamma2_1=aux_matrix1*([gamma2_1_tilde;zeros(numExo,1)]);
% gamma3_1=aux_matrix1*([(AA1_inv*EE1);FF1]);
% 
% aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1);
% gamma1_2=aux_matrix2*([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
% gamma2_2=aux_matrix2*([gamma2_2_tilde;zeros(numExo,1)]);
% gamma3_2=aux_matrix2*([(AA2_inv*EE2);FF2]);


