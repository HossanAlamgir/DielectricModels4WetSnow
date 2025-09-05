function [time1,data2_1] = matchup_time_series(data1,time1,data2,time2)
%%%
% This simple function matches two time series within a given temporal window
% Alamgir Hossan
% 07/24/2024
% alamgir.hossan@jpl.nasa.gov
%%%
if ~isempty(data1)
data2_1 = nan(size(data1,1), size(data2,2));

for i = 1:length(time1)
    deltaT = time2 - time1(i);
    ind = find(abs(deltaT) < 1);
    [m,si] = min(abs(deltaT(ind)));
    if ~isempty(ind)
        data2_1(i,:) = data2(ind(si),:);
    end

end
else
    time1 = time2;
    data2_1 = data2;


end