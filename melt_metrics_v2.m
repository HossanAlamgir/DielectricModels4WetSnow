function [monset,mfreezup,mduration,maxsmelt] = melt_metrics_v2(days,LWA,min_thresh,yr)
% This function calculates melt onset, freeze dates, max summer melt and
% duration
%Alamgir Hossan, JPL, 9/10/2024

monset = nan;
mfreezup = nan;
mduration = nan;
maxsmelt = nan;

ind_melt = find(LWA>min_thresh);
if ~isempty(ind_melt)
    monset = days(ind_melt(1));
    mfreezup = days(ind_melt(end));
    mduration = mfreezup + 1 - monset;
    maxsmelt = max(LWA);

    % monset = datestr(monset,'mmm dd');
    % mfreezup = datestr(mfreezup,'mmm dd');
    % monset = monset - datenum(yr,1,1) + 1 ;
    % mfreezup = mfreezup - datenum(yr,1,1) + 1 ;

end

end