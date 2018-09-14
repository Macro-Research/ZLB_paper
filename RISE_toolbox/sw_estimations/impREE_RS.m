clear;clc;close all;
load('estimation_results_rs.mat');
myirfs0=irf(sw,'irf_periods',100);
clearvars -except myirfs0;

impREE_normal(:,1,1)=double(myirfs0.ea.c(:,1));
impREE_normal(:,1,2)=double(myirfs0.eb.c(:,1));
impREE_normal(:,1,3)=double(myirfs0.eg.c(:,1));
impREE_normal(:,1,4)=double(myirfs0.eqs.c(:,1));
impREE_normal(:,1,5)=double(myirfs0.epinf.c(:,1));
impREE_normal(:,1,6)=double(myirfs0.ew.c(:,1));

impREE_normal(:,2,1)=double(myirfs0.ea.inve(:,1));
impREE_normal(:,2,2)=double(myirfs0.eb.inve(:,1));
impREE_normal(:,2,3)=double(myirfs0.eg.inve(:,1));
impREE_normal(:,2,4)=double(myirfs0.eqs.inve(:,1));
impREE_normal(:,2,5)=double(myirfs0.epinf.inve(:,1));
impREE_normal(:,2,6)=double(myirfs0.ew.inve(:,1));

impREE_normal(:,3,1)=double(myirfs0.ea.y(:,1));
impREE_normal(:,3,2)=double(myirfs0.eb.y(:,1));
impREE_normal(:,3,3)=double(myirfs0.eg.y(:,1));
impREE_normal(:,3,4)=double(myirfs0.eqs.y(:,1));
impREE_normal(:,3,5)=double(myirfs0.epinf.y(:,1));
impREE_normal(:,3,6)=double(myirfs0.ew.y(:,1));

impREE_normal(:,4,1)=double(myirfs0.ea.pinf(:,1));
impREE_normal(:,4,2)=double(myirfs0.eb.pinf(:,1));
impREE_normal(:,4,3)=double(myirfs0.eg.pinf(:,1));
impREE_normal(:,4,4)=double(myirfs0.eqs.pinf(:,1));
impREE_normal(:,4,5)=double(myirfs0.epinf.pinf(:,1));
impREE_normal(:,4,6)=double(myirfs0.ew.pinf(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

impREE_zlb(:,1,1)=double(myirfs0.ea.c(:,2));
impREE_zlb(:,1,2)=double(myirfs0.eb.c(:,2));
impREE_zlb(:,1,3)=double(myirfs0.eg.c(:,2));
impREE_zlb(:,1,4)=double(myirfs0.eqs.c(:,2));
impREE_zlb(:,1,5)=double(myirfs0.epinf.c(:,2));
impREE_zlb(:,1,6)=double(myirfs0.ew.c(:,2));

impREE_zlb(:,2,1)=double(myirfs0.ea.inve(:,2));
impREE_zlb(:,2,2)=double(myirfs0.eb.inve(:,2));
impREE_zlb(:,2,3)=double(myirfs0.eg.inve(:,2));
impREE_zlb(:,2,4)=double(myirfs0.eqs.inve(:,2));
impREE_zlb(:,2,5)=double(myirfs0.epinf.inve(:,2));
impREE_zlb(:,2,6)=double(myirfs0.ew.inve(:,2));

impREE_zlb(:,3,1)=double(myirfs0.ea.y(:,2));
impREE_zlb(:,3,2)=double(myirfs0.eb.y(:,2));
impREE_zlb(:,3,3)=double(myirfs0.eg.y(:,2));
impREE_zlb(:,3,4)=double(myirfs0.eqs.y(:,2));
impREE_zlb(:,3,5)=double(myirfs0.epinf.y(:,2));
impREE_zlb(:,3,6)=double(myirfs0.ew.y(:,2));

impREE_zlb(:,4,1)=double(myirfs0.ea.pinf(:,2));
impREE_zlb(:,4,2)=double(myirfs0.eb.pinf(:,2));
impREE_zlb(:,4,3)=double(myirfs0.eg.pinf(:,2));
impREE_zlb(:,4,4)=double(myirfs0.eqs.pinf(:,2));
impREE_zlb(:,4,5)=double(myirfs0.epinf.pinf(:,2));
impREE_zlb(:,4,6)=double(myirfs0.ew.pinf(:,2));


save impREE.mat impREE_normal impREE_zlb;