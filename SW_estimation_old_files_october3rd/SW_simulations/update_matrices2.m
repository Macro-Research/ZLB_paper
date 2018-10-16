msgstr=[];
msgid=[];
lastwarn('Success');

gamma1_1_tilde=AA1\(BB1+CC1*beta_tt^2);
gamma2_1_tilde=(AA1\CC1)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_1_tilde=AA1\(CC1*cc_tt+DD1);
gamma4_1_tilde=AA1\(CC1*beta_tt*cc_tt);

gamma1_2_tilde=AA2\(BB2+CC2*beta_tt^2);
gamma2_2_tilde=(AA2\CC2)*(eye(numEndo)+beta_tt)*alpha_tt;
gamma3_2_tilde=AA2\(CC2*cc_tt+DD2);
gamma4_2_tilde=AA2\(CC2*beta_tt*cc_tt);

aux_matrix1=[eye(numEndo),-gamma3_1_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_1=aux_matrix1\([gamma1_1_tilde,gamma4_1_tilde,;zeros(numExo,numEndo),RHO1]);
gamma2_1=aux_matrix1\([gamma2_1_tilde;zeros(numExo,1)]);
gamma3_1=aux_matrix1\([(AA1\EE1);FF1]);

aux_matrix2=[eye(numEndo),-gamma3_2_tilde;zeros(numExo,numEndo),eye(numExo)];
gamma1_2=aux_matrix2\([gamma1_2_tilde,gamma4_2_tilde,;zeros(numExo,numEndo),RHO2]);
gamma2_2=aux_matrix2\([gamma2_2_tilde;zeros(numExo,1)]);
gamma3_2=aux_matrix2\([(AA2\EE2);FF2]);

[msgstr, msgid] = lastwarn;