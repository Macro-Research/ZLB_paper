st_dev=double([oo_.posterior_mode.shocks_std.eta_a,oo_.posterior_mode.shocks_std.eta_b,...
    oo_.posterior_mode.shocks_std.eta_g,oo_.posterior_mode.shocks_std.eta_i,...
    oo_.posterior_mode.shocks_std.eta_p,oo_.posterior_mode.shocks_std.eta_w]);

impREE_unitShock(:,1,1)=c_eta_a/st_dev(1);
impREE_unitShock(:,1,2)=c_eta_b/st_dev(2);
impREE_unitShock(:,1,3)=c_eta_g/st_dev(3);
impREE_unitShock(:,1,4)=c_eta_i/st_dev(4);
impREE_unitShock(:,1,5)=c_eta_p/st_dev(5);
impREE_unitShock(:,1,6)=c_eta_w/st_dev(6);

impREE_unitShock(:,2,1)=inve_eta_a/st_dev(1);
impREE_unitShock(:,2,2)=inve_eta_b/st_dev(2);
impREE_unitShock(:,2,3)=inve_eta_g/st_dev(3);
impREE_unitShock(:,2,4)=inve_eta_i/st_dev(4);
impREE_unitShock(:,2,5)=inve_eta_p/st_dev(5);
impREE_unitShock(:,2,6)=inve_eta_w/st_dev(6);

impREE_unitShock(:,3,1)=y_eta_a/st_dev(1);
impREE_unitShock(:,3,2)=y_eta_b/st_dev(2);
impREE_unitShock(:,3,3)=y_eta_g/st_dev(3);
impREE_unitShock(:,3,4)=y_eta_i/st_dev(4);
impREE_unitShock(:,3,5)=y_eta_p/st_dev(5);
impREE_unitShock(:,3,6)=y_eta_w/st_dev(6);

impREE_unitShock(:,4,1)=pinf_eta_a/st_dev(1);
impREE_unitShock(:,4,2)=pinf_eta_b/st_dev(2);
impREE_unitShock(:,4,3)=pinf_eta_g/st_dev(3);
impREE_unitShock(:,4,4)=pinf_eta_i/st_dev(4);
impREE_unitShock(:,4,5)=pinf_eta_p/st_dev(5);
impREE_unitShock(:,4,6)=pinf_eta_w/st_dev(6);
%------------------------------------------------------
impREE_stdev(:,1,1)=c_eta_a;
impREE_stdev(:,1,2)=c_eta_b;
impREE_stdev(:,1,3)=c_eta_g;
impREE_stdev(:,1,4)=c_eta_i;
impREE_stdev(:,1,5)=c_eta_p;
impREE_stdev(:,1,6)=c_eta_w;

impREE_stdev(:,2,1)=inve_eta_a;
impREE_stdev(:,2,2)=inve_eta_b;
impREE_stdev(:,2,3)=inve_eta_g;
impREE_stdev(:,2,4)=inve_eta_i;
impREE_stdev(:,2,5)=inve_eta_p;
impREE_stdev(:,2,6)=inve_eta_w;

impREE_stdev(:,3,1)=y_eta_a;
impREE_stdev(:,3,2)=y_eta_b;
impREE_stdev(:,3,3)=y_eta_g;
impREE_stdev(:,3,4)=y_eta_i;
impREE_stdev(:,3,5)=y_eta_p;
impREE_stdev(:,3,6)=y_eta_w;

impREE_stdev(:,4,1)=pinf_eta_a;
impREE_stdev(:,4,2)=pinf_eta_b;
impREE_stdev(:,4,3)=pinf_eta_g;
impREE_stdev(:,4,4)=pinf_eta_i;
impREE_stdev(:,4,5)=pinf_eta_p;
impREE_stdev(:,4,6)=pinf_eta_w;


save impREE.mat impREE_unitShock impREE_stdev;
