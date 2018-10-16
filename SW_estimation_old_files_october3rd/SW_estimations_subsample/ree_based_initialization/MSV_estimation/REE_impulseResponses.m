st_dev=double([oo_.posterior_mode.shocks_std.eta_a,oo_.posterior_mode.shocks_std.eta_b,...
    oo_.posterior_mode.shocks_std.eta_g,oo_.posterior_mode.shocks_std.eta_i,...
    oo_.posterior_mode.shocks_std.eta_p,oo_.posterior_mode.shocks_std.eta_w]);

impREE(:,1,1)=c_eta_a/st_dev(1);
impREE(:,1,2)=c_eta_b/st_dev(2);
impREE(:,1,3)=c_eta_g/st_dev(3);
impREE(:,1,4)=c_eta_i/st_dev(4);
impREE(:,1,5)=c_eta_p/st_dev(5);
impREE(:,1,6)=c_eta_w/st_dev(6);

impREE(:,2,1)=inve_eta_a/st_dev(1);
impREE(:,2,2)=inve_eta_b/st_dev(2);
impREE(:,2,3)=inve_eta_g/st_dev(3);
impREE(:,2,4)=inve_eta_i/st_dev(4);
impREE(:,2,5)=inve_eta_p/st_dev(5);
impREE(:,2,6)=inve_eta_w/st_dev(6);

impREE(:,3,1)=y_eta_a/st_dev(1);
impREE(:,3,2)=y_eta_b/st_dev(2);
impREE(:,3,3)=y_eta_g/st_dev(3);
impREE(:,3,4)=y_eta_i/st_dev(4);
impREE(:,3,5)=y_eta_p/st_dev(5);
impREE(:,3,6)=y_eta_w/st_dev(6);

impREE(:,4,1)=pinf_eta_a/st_dev(1);
impREE(:,4,2)=pinf_eta_b/st_dev(2);
impREE(:,4,3)=pinf_eta_g/st_dev(3);
impREE(:,4,4)=pinf_eta_i/st_dev(4);
impREE(:,4,5)=pinf_eta_p/st_dev(5);
impREE(:,4,6)=pinf_eta_w/st_dev(6);


save impREE.mat impREE;
