function[likl postLikl] = NKPC_objFcn(param, phi_smc,bounds,dataset)
    


likl=likelihoodKalman(param,dataset);
priorLikl=prior_dist(param);
postLikl=phi_smc*likl+priorLikl;

parabd_ind1 = param > bounds(:,1); % lower bounds
parabd_ind2 = param < bounds(:,2); % upper bounds
 parabd_ind1=sum(parabd_ind1);
 parabd_ind2=sum(parabd_ind2);

likl=-likl;postLikl=-postLikl;
% if parabd_ind1<length(param) || parabd_ind2<length(param)
% likl=-Inf;postLikl=-Inf;
% end



end




