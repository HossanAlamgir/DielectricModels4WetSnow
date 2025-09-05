clear
cluster=1;

% Disable all warnings
warning('off', 'all');

if cluster==0
    input_drive = 'Z:';
    output_drive = 'Z:';

    addpath(fullfile([input_drive, '/MEMLS3&a/AH']))
    addpath(fullfile([input_drive, '/IceSheets/firn_models/Marshall_SJ_2023']))
    addpath(fullfile(input_drive, 'IceSheets', 'Documentation', 'Manuscripts_n_abstracts', 'wet_snow_diel_models', 'snow_mixing_formula'));

elseif cluster==1

    input_drive = '/projects/lotus';
    output_drive = '/projects/lotus';
    addpath('/projects/lotus/MEMLS3&a/AH/')
    addpath('/projects/lotus/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/codes/snow_mixing_formula/');
    addpath('/projects/lotus/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/codes/');
end

load(fullfile([input_drive, '/IceSheets/firn_models/Marshall_SJ_2023/z.mat'])) % EMB model vertical spacing

result_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/results/'];

retrieval=1;
pol='V';
% yrs = 2023:-1:2010;
% yrs = 2010:2023;
yrs = 2015:2023;

% Load SUMUP density data at PROMICE AWS
load( [input_drive,'/IceSheets/Weather Stations Data/Greenland/PROMICE/SUMup_2023_beta/SUMup_2023_density_GrIS_PROMICE_6AWS.mat']);
%% Load SMAP TB Data
% load([input_dir,sprintf('Ka_band_TBs_at_PROMICE_AWS.mat')])
load([result_dir,sprintf('SMAP_TBs_PROMICE_AWS.mat')])
stations = aws.name;
%% Load Look Up Tables
load([result_dir,'LUT_memlsv3.mat']); % RevC
theta_0=40;%SMAP Incident angle
f=1.41;
%%
for model = 1:10

    if pol=='V'
        TBV_LUT_memls = squeeze(tbv(:,:,:,:,model)); % L-band LUT, V-pol
    elseif pol=='H'
        TBH_LUT_memls = squeeze(tbh(:,:,:,model)); % L-band LUT, H-pol, not used
    end
    data_smap_ret = cell(1,length(stations));
    data_smap_retc = cell(1,length(stations));

    for stn = 1:length(stations)
        % try
        tic
        data_smap = [];
        data_smapc = [];

        smap_data = [smap.dnum(:,stn),smap.Tbv(:,stn),smap.Tbh(:,stn)];

        for yr=1:length(yrs)
            [smap_datay,dnum]=return_desired_data_revB(smap_data,smap_data(:,1),yrs(yr));
            Tbv = smap_datay(:,2);
            Tbh = smap_datay(:,3);

            if isempty(Tbv)
                continue
            else
                fprintf('%s(%d)\n',stations{stn},yrs(yr))
            end

            Z=10;
            Thresh_max =16;
            Tcap = 0;
            TBV_rSIR{1}=Tbv;
            TBH_rSIR{1}=Tbh;
            dnum_rSIR{1}=dnum;
            [tbv_anom_refw, ~, ref_tbv, ~, std_tbv, ~,TBV_fall_bias,TBH_fall_bias] ...
                = compute_spring_or_fall_reference(TBV_rSIR,TBH_rSIR,dnum_rSIR, yrs(yr), Z,Tcap);
            tbv_anom_refw=tbv_anom_refw{1};
            thresh_tbv = min(Thresh_max,Z*std_tbv);
            meltF = abs(tbv_anom_refw)> thresh_tbv;
            clear TBV_rSIR TBH_rSIR dnum_rSIR

            [rho0z,rhob0z] = return_density_profile_revA(density_aws{stn},yrs(yr),z);
            ind = find(z>3.0);
            mrho = nanmean(rho0z(1:ind(1)));
            %% Inversion
            y=1;
            if retrieval==1

                if pol=='V'
                    Tb = Tbv;
                elseif pol=='H'
                    Tb = Tbh;
                end

                [mv2_est{y},TBH_est{y},TBV_est{y},Ps_est{y},twet2_est{y},std_twet{y},...
                    eps2r_est{y},d_dry_est{y},ref_tbv_adj{y}]= estimate_lwa_with_diff_diel_models_v2(Tb,TBV_LUT_memls,f,theta_0,par_range_LUT,ref_tbv{y},std_tbv,TBV_fall_bias,TBH_fall_bias,meltF,mrho, model);

            end
            data = [dnum,Tbv,Tbh, ref_tbv_adj{y},ones(size(Tbv))*std_tbv,mv2_est{y}',ones(size(Tbv))*twet2_est{y},TBV_est{y}', TBH_est{y}'];
            data_smap = [data_smap;data];

            datac = [stn,eps2r_est{y},d_dry_est{y},Ps_est{y},twet2_est{y},std_twet{y}];
            data_smapc = [data_smapc;datac];


        end

        data_smap_ret{stn}=data_smap;
        data_smap_retc{stn}=data_smapc;
        stations_ret{stn}=char(stations{stn});
        fprintf('Elapsed Time for processing %s\n',char(stations{stn}))
        toc

        % catch
    end

    save([result_dir,'lwa/v4/smap_lwa_model_',num2str(model),'_avg_twet_avg_pstop3m.mat'],'data_smap_ret','data_smap_retc','stations_ret','par_range_LUT','-v7.3')

end




