if largest_eig2(tt)>.9999
   alpha_tt=alpha_old;
    beta_tt=beta_old;
    cc_tt=cc_old;
    pr_flag(tt)=1;
    
elseif strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
    alpha_tt=alpha_old;
    beta_tt=beta_old;
    cc_tt=cc_old;
    pr_flag(tt)=1;
    
    
 end