function [likelihood]= likelihoodKalman(theta,dataset)

pr=prior_dist(theta);
likl=kalmanNKPC_measError(theta,dataset);
likelihood=likl+pr;

end