function [likelihood]= likelihood(theta)
% msgstr=[];
% msgid=[];
% penalty=50;
% load('bounds.mat');
% cond1=sum(double(theta>bounds(:,1)));
% cond2=sum(double(theta<bounds(:,2)));

pr=prior_dist(theta);
likl=KF_MSV(theta);
likelihood=likl+pr;


% [msgstr, msgid] = lastwarn;
%  if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1;
%      likl=likl+penalty;
%  end
% 
% if cond1<length(theta)
%     likelihood=Inf;
% elseif cond2<length(theta)
%     likelihood=Inf;
% end

end