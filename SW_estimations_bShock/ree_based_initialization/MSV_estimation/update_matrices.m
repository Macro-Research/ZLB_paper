msgstr=[];
msgid=[];
lastwarn('Success');

gamma1_1_tilde=AA1_inv*(BB1+CC1*beta_tt^2);
gamma2_1_tilde=AA1_inv*CC1*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_1_tilde=(AA1_inv*CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+AA1_inv*DD1;

gamma1_2_tilde=AA2_inv*(BB2+CC2*beta_tt^2);
gamma2_2_tilde=AA2_inv*CC2*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_2_tilde=(AA2_inv*CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2_inv*DD2;

aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1);
gamma1_1=aux_matrix1*[gamma1_1_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO1];
gamma2_1=aux_matrix1*[gamma2_1_tilde;zeros(numExo,1)];
gamma3_1=aux_matrix1*[AA1_inv*EE1;FF1];

aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)]^(-1);
gamma1_2=aux_matrix2*[gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2];
gamma2_2=aux_matrix2*[gamma2_2_tilde;zeros(numExo,1)];
gamma3_2=aux_matrix2*[AA2_inv*EE2;FF2];
 
[AA1_aux,BB1_aux,CC1_aux]=SW_sysmat_VAR_filter(parameters(1:end,1));
[AA2_aux,BB2_aux,CC2_aux]=SW_sysmat_VAR_filter(parameters(1:end,2));
aux1=AA1_aux^(-1)*(BB1_aux);
aux2=AA2_aux^(-1)*(BB2_aux);

gamma1_1(:,18)=aux1(:,18);
gamma1_2(:,18)=aux2(:,18);

% gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
% gamma2_1_tilde=(AA1\CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
% gamma3_1_tilde=(AA1\CC1)*(beta_tt*cc_tt+cc_tt*RHO1)+(AA1\DD1);
% 
% gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);
% gamma2_2_tilde=(AA2\CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
% gamma3_2_tilde=(AA2\CC2)*(beta_tt*cc_tt+cc_tt*RHO2)+AA2\DD2;
% 
% aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)];
% gamma1_1=aux_matrix1\([gamma1_1_tilde,zeros(numEndo,numExo);zeros(numExo,numEndo),RHO1]);
% gamma2_1=aux_matrix1\([gamma2_1_tilde;zeros(numExo,1)]);
% gamma3_1=aux_matrix1\([(AA1\EE1);FF1]);
% 
% aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)];
% gamma1_2=aux_matrix2\([gamma1_2_tilde,zeros(numEndo,numExo),;zeros(numExo,numEndo),RHO2]);
% gamma2_2=aux_matrix2\([gamma2_2_tilde;zeros(numExo,1)]);
% gamma3_2=aux_matrix2\([(AA2\EE2);FF2]);

[msgstr, msgid] = lastwarn;
