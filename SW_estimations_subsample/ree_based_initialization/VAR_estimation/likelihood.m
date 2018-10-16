function [likelihood]= likelihood(theta)

pr=prior_dist(theta);
likl=KF_VAR(theta);
likelihood=likl+pr;

end