if max(largest_eig2(tt),largest_eig(tt))>0.999
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