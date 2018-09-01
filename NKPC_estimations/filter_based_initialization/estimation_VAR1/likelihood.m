function [likelihood]= likelihood(param)


pr=NKPC_prior(param);
likl=KF_MS_VAR1(param);
likelihood=likl+pr;

end