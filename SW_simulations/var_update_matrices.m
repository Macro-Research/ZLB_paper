
   gamma1_1=AA1_inv * (( BB1+CC1*beta_tt^2));
   gamma1_2=AA1_inv*CC1*(eye(numVar)+beta_tt)*alpha_tt;
   gamma1_3=AA1_inv * DD1;
   
   gamma2_1=AA2_inv * (( BB2+CC2*beta_tt^2));
   gamma2_2=AA2_inv*CC2*(eye(numVar)+beta_tt)*alpha_tt;
   gamma2_3=AA2_inv * DD2;