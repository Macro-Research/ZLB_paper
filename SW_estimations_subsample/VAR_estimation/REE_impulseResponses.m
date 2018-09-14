impREE(:,1,1)=c_eta_a;
impREE(:,1,2)=c_eta_b;
impREE(:,1,3)=c_eta_g;
impREE(:,1,4)=c_eta_i;
impREE(:,1,5)=c_eta_p;
impREE(:,1,6)=c_eta_w;

impREE(:,2,1)=inve_eta_a;
impREE(:,2,2)=inve_eta_b;
impREE(:,2,3)=inve_eta_g;
impREE(:,2,4)=inve_eta_i;
impREE(:,2,5)=inve_eta_p;
impREE(:,2,6)=inve_eta_w;

impREE(:,3,1)=y_eta_a;
impREE(:,3,2)=y_eta_b;
impREE(:,3,3)=y_eta_g;
impREE(:,3,4)=y_eta_i;
impREE(:,3,5)=y_eta_p;
impREE(:,3,6)=y_eta_w;

impREE(:,4,1)=pinf_eta_a;
impREE(:,4,2)=pinf_eta_b;
impREE(:,4,3)=pinf_eta_g;
impREE(:,4,4)=pinf_eta_i;
impREE(:,4,5)=pinf_eta_p;
impREE(:,4,6)=pinf_eta_w;


save impREE.mat impREE;
