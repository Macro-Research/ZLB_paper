function [likelihood,likl_increment]= likelihood_for_hessian(theta)

pr=prior_dist(theta);
[likl,likl_increment]=KF_VAR_for_hessian(theta);
likelihood=likl+pr;

end