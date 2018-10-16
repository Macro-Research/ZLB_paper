function [likelihood]= likelihoodSPF(theta,dataset)

pr=priorDist_SPF_actualData(theta);
likl=sigmaPointFilter_NKPC(theta,dataset);
likelihood=likl+pr;

end