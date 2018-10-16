function [likelihood]= likelihood(theta)

% penalty=50;
% load('bounds.mat');
% cond1=sum(double(theta>bounds(:,1)));
% cond2=sum(double(theta<bounds(:,2)));

pr=prior_dist(theta);
likl=KF_MSV(theta);
likelihood=likl+pr;



% 
% if cond1<length(theta)
%     likelihood=Inf;
% elseif cond2<length(theta)
%     likelihood=Inf;
% end

end