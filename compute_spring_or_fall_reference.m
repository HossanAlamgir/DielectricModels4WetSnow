function [dTBV, dTBH, TBV_ref, TBH_ref, TBV_thres, TBH_thres,TBV_diff,TBH_diff] = compute_spring_or_fall_reference(TBV, TBH, dnum, yrs, Z,Tcap)

%%
for m = 1:length(yrs)

    dnum_spring1(m) = datenum(yrs(m),  1,  1); % Changed back to begining of the season
    dnum_spring2(m) = datenum(yrs(m),  3,  31);  % for taking SMAP case in 2015

    dnum_fall1(m)   = datenum(yrs(m), 11, 01);
    dnum_fall2(m)   = datenum(yrs(m), 12, 31); % Changed to the end of the season

end

%%
for n = 1:length(TBV)

    TBV_ref{n} = NaN(size(TBV{n}));
    TBH_ref{n} = NaN(size(TBH{n}));
    dTBV{n} = NaN(size(TBV{n}));
    dTBH{n} = NaN(size(TBH{n}));

    for m = 1:length(yrs)


        dnum1               = datenum(yrs(m),  1, 1);
        dnum2               = datenum(yrs(m), 12,31);

        b_yr                = dnum{n} >= dnum1 & dnum{n} < dnum2;

        b_spring{n}         = dnum{n} >= dnum_spring1(m) & dnum{n} <= dnum_spring2(m);
        b_fall{n}           = dnum{n} >= dnum_fall1(m)   & dnum{n} <= dnum_fall2(m);

        TBV_spring(m,n)     = nanmean(TBV{n}(b_spring{n}));
        TBH_spring(m,n)     = nanmean(TBH{n}(b_spring{n}));

        TBV_fall(m,n)       = nanmean(TBV{n}(b_fall{n}));
        TBH_fall(m,n)       = nanmean(TBH{n}(b_fall{n}));

        TBV_spring_std(m,n) = nanstd(TBV{n}(b_spring{n}));
        TBH_spring_std(m,n) = nanstd(TBH{n}(b_spring{n}));

        TBV_fall_std(m,n)   = nanstd(TBV{n}(b_fall{n}));
        TBH_fall_std(m,n)   = nanstd(TBH{n}(b_fall{n}));

        dnumA = dnum{n}(b_yr);
        if ~all(isnan(TBV_spring))
            TBV_ref{n}(:) = TBV_spring(m,n); % Initial Ref
            TBH_ref{n}(:) = TBH_spring(m,n);
        else
            TBV_ref{n}(:) = TBV_fall(m,n); % Initial Ref
            TBH_ref{n}(:) = TBH_fall(m,n);
        end

        dTBV{n}             = TBV{n} - TBV_ref{n}; % Initial anomaly
        dTBH{n}             = TBH{n} - TBH_ref{n};

        TBV_diff(m,n)   = TBV_spring(m,n) - TBV_fall(m,n);
        TBH_diff(m,n)   = TBH_spring(m,n) - TBH_fall(m,n);

        
    end
end

TBV_diff(TBV_diff< 0) = 0;
TBH_diff(TBH_diff< 0) = 0;

% TBV_thres = sqrt((TBV_spring_std*Z).^2 + (TBV_fall_std*Z).^2); % Misses the first melt
% TBH_thres = sqrt((TBH_spring_std*Z).^2 + (TBH_fall_std*Z).^2);

% TBV_thres = Z*(TBV_spring_std+ TBV_fall_std)/2; % Fall STD is too high
% TBH_thres = Z*(TBH_spring_std + TBH_fall_std)/2;

% TBV_thres = Z*TBV_spring_std; % Fall STD is too high
% TBH_thres = Z*TBH_spring_std;

% TBV_thres = min(Z*TBV_spring_std,Tcap); % Fall STD is too high
% TBH_thres = min(Z*TBH_spring_std,Tcap);
if ~all(isnan(TBV_spring_std))
    TBV_thres = TBV_spring_std; % Threshold to be computed later
    TBH_thres = TBH_spring_std;
else

    TBV_thres = TBV_fall_std; % Threshold to be computed later
    TBH_thres = TBH_fall_std;
end

end




