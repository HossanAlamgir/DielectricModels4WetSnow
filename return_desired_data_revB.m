function [data_desiredm,dnum_desiredm]=return_desired_data_revB(data,dnum,desired_yr,desired_month)
if ~exist('desired_month','var')
    desired_month = [1 12];
end
dnum = double(dnum);
yr_initial = min(desired_yr);
yr_final = max(desired_yr);
mon_initial = min(desired_month);
mon_final = max(desired_month);
ind = find(dnum>=datenum(yr_initial,1,1) & dnum<=datenum(yr_final,12,31));
data_desired = data(ind,:);
dnum_desired = dnum(ind);

dvec = datevec(dnum(ind));

indm = find(dvec(:,2)>=mon_initial & dvec(:,2)<=mon_final);
data_desiredm = data_desired(indm,:);
dnum_desiredm = dnum_desired(indm);
end