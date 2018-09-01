function [likelihood]= likelihood(param)


pr=NKPC_prior(param);
likl=KF_MS_MSV(param);
likelihood=likl+pr;

end