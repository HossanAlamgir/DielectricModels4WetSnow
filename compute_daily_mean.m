function [days,quantity_daily]=compute_daily_mean(dnum,quantity)

quantity_daily = [];
startdate = ceil(min(dnum));
enddate   = floor(max(dnum))-1;

days     = startdate:enddate;

for m = 1:length(days)
    b                        = days(m) == floor(dnum);
    quantity_daily(m,1) = nanmean(quantity(b));
end



end