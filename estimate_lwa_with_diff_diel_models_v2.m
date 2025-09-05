function [mv_est,TBH_est,TBV_est,Ps,twet2_est,std_twet,eps3r,d_dry,ref_tbv] = estimate_lwa_with_diff_diel_models_v2(tbv,TBV_LUT,f,theta_0,par_range_LUT,ref_tbv,ref_std,fall_bias_tbv,fall_bias_tbh,meltf,mrho,model)

%This code estimates the mv2, TBH, TBV, and NPR for the case of NPR
%decrease. It used the "IceSheet_EM_Model_nlayer.m" function.
%
%
%%%--- INPUTS:
%%%--- dnum_pm_smap: time stamp from SMAP data
%%%--- npr_pm_smap: NPR from SMAP data
%%%--- tbh_pm_smap: TBH from SMAP data
%%%--- tbv_pm_smap: TBV from SMAP data
%%%--- T0_wet2: physical temperature the wet snow layer 2 (K)
%%%--- NPR_LUT_DECR: LUT NPR value
%%%--- TBH_LUT_DECR: LUT TBH value
%%%--- TBV_LUT_DECR: LUT TBV value
%%%--- mv2: LUT mv2 value
%%%--- Ps: LUT Ps value
%%%--- d_wet2: LUT d_wet2 value
%%%--- eps3r: LUT eps3r value
%%%--- eps3i: LUT eps3i value
%%%--- T0_dry3: LUT T0_dry3 value
%%%--- d_dry3: LUT d_dry3 value
%%%--- T0_dry4: LUT T0_dry4 value
%%%--- f: frequency of operation (GHz)
%%%--- theta_0: incident (observation) angle in degrees
%%%
%%%--- OUTPUTS:
%%%--- mv2_est_daily: estimated snow wetness percentage
%%%--- TBH_LUT_INCR: incoherent H-pol brightness temperature LUT-INCR (K)
%%%--- TBV_LUT_INCR: incoherent V-pol brightness temperature LUT-INCR (K)
%%%--- NPR_LUT_INCR: incoherent normalized polarization ratio LUT-INCR (NPR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,I_d_dry2] = min(abs(par_range_LUT.d_dry2 - 500)); % use a fixed dry layer thickness
[~,I_Ps] = min(abs(par_range_LUT.Ps - mrho/1000)); % use measured dry snow density
%%STEP 1: Retrieve Frozen Season Parameters
TBV_LUT_FS = squeeze(TBV_LUT(I_Ps,I_d_dry2,:,1)); % LUT with mv = 0;
cost_func_TBV=abs(TBV_LUT_FS-ref_tbv(1));
[delta,I_tbv] = min(cost_func_TBV(:));
[I_eps3r] = ind2sub(size(TBV_LUT_FS),I_tbv);
Ps=par_range_LUT.Ps(I_Ps); %test case
T0_dry2=par_range_LUT.T0_dry2;
d_dry=par_range_LUT.d_dry2(I_d_dry2);
eps3r=par_range_LUT.eps3r(I_eps3r);
eps3i=par_range_LUT.eps3i;
d3 = par_range_LUT.d3;
T0_dry4=par_range_LUT.T0_dry4;
% TB_dry_sim = TBV_LUT_FS(I_tbv);

ref_tbv_fall=ref_tbv(1)-fall_bias_tbv;
cost_func_TBV=abs(TBV_LUT_FS-ref_tbv_fall);
[delta,I_tbv_fall] = min(cost_func_TBV(:));
[I_eps3r_fall] = ind2sub(size(TBV_LUT_FS),I_tbv_fall);
eps3r_fall=par_range_LUT.eps3r(I_eps3r_fall);

clear  cost_func_TBV

%%STEP 2: Estimation of wet layer thickness from melt season TB, using frozen season
%%parameters as priori
tbvm=tbv(meltf);
% % tbvm=tbv;
for dd=1:length(tbvm)
    cost_func_TBV=abs(TBV_LUT(I_Ps,:,I_eps3r,:)-tbvm(dd));
    [delta_TB(dd),I_tbv] = min(cost_func_TBV(:));
    [~,I_dwet2,~,I_mv2] = ind2sub(size(TBV_LUT(I_Ps,:,I_eps3r,:)),I_tbv);
    mv2_est2(dd)=par_range_LUT.mv(I_mv2);
    twet2_est(dd)=par_range_LUT.d_dry2(I_dwet2);
end


% % Use const volume fraction
% I_mv2 =find(par_range_LUT.mv==2.0);% mv=2%
% for dd=1:length(tbvm)
%     test_func_TBV=abs(TBV_LUT(I_Ps,:,I_eps3r,I_mv2)-tbvm(dd));
%     [delta_TB(dd),I_dwet2] = min(test_func_TBV(:));
%     mv2_est2(dd)=par_range_LUT.mv(I_mv2);
%     dwet2_est2(dd)=par_range_LUT.d_dry2(I_dwet2);
% end

% % Use fixed twet
% dwet2_est2 = 200; % try fixed twet 2 m

% % Use avg twet
std_twet=nanstd(twet2_est);
twet2_est=nanmean(twet2_est);

% % Use cmedian twet
% dwet2_est2 = median(dwet2_est2,"omitmissing");

% % mv2_est_daily=mv2_est2;

%% STEP 3: Re-estimate mv2 for each day.
% tbv=tbv(melt_days);
% mv2_new=0:0.05:mv2_est2;
% Determine if there is negative difference between spring and fall TB, possibly caused by refrozen crust and
% account for that
delta_ref = tbv - ref_tbv(1);
[tbv_max,ind_tbv_max]=max(tbv);
% day_fall_bias = find(delta_ref<(-3*ref_std));
% if ~isempty(day_fall_bias) 
%     ind=find(day_fall_bias>ind_tbv_max,1);
%     day_fall_bias = day_fall_bias(ind);
% end


[epsr_ds, epsi_ds] = diel_dry_snow(T0_dry4-273, 0.917, f );
eps4=epsr_ds-1i*epsi_ds;

mv2_new=linspace(0,6,400);
for dd=1:length(tbv)
    if  dd<=ind_tbv_max
        eps3=eps3r-1i*eps3i;
    else
        eps3=eps3r_fall-1i*eps3i;
        ref_tbv(dd)=ref_tbv_fall;
    end

    for ii=1:length(mv2_new)
        if mv2_new(ii)==0
            T2=T0_dry2;
            T3=T0_dry2;
            [epsr_ws,epsi_ws] = diel_dry_snow(T0_dry2-273,  Ps, f );
            eps2=epsr_ws-1i*epsi_ws;
            d2 = d_dry;
        elseif mv2_new(ii)>0
            T2 = 273.15;
            T3=265;
            d2 = twet2_est;

            [epsr_ws,epsi_ws] = return_wet_snow_diel_const_v4(f, Ps, mv2_new(ii),model);
            eps2=epsr_ws-1i*epsi_ws;
        end       
        
        eps=[eps2,eps3,eps4];

        [TBH_incoh(dd,ii),TBV_incoh(dd,ii)] =return_memls_TB_v2(f,theta_0,[d2,d3,100000]',[T2,T3,T0_dry4]',[1000*Ps,1000*Ps,917]',[mv2_new(ii)/100,0,0]',eps');
    end

    cost_func_TBV=abs(TBV_incoh(dd,:)-(tbv(dd)));
    [~,I_tbv] = min(cost_func_TBV(:));
    mv_est(dd)=mv2_new(I_tbv);
    clear cost_func_TBV

    TBV_est(dd) = TBV_incoh(dd,I_tbv);
    TBH_est(dd) = TBH_incoh(dd,I_tbv);

end


end


