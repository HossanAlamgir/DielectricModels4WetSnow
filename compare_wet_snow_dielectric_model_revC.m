clear
cluster=0;

% Disable all warnings
warning('off', 'all');

if cluster==0
    input_drive = 'Z:';
    output_drive = 'Z:';
elseif cluster==1
    input_drive = '/projects/lotus';
    output_drive = '/projects/lotus';
end

fig_dir = 'Z:\IceSheets\Documentation\Manuscripts_n_abstracts\wet_snow_diel_models\figs\revised\';
result_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/results/'];
models = [{'Mätzler'},{'Tinga'},{'Debye-like'},{'Hallikainen'},{'Ulaby'},...
    {'Colbeck'},{'Birchak'},{'Sihvola'},{'Looyenga'},{'Tiuri'}];
styleSet=plot_styles();

addpath([input_drive, '/IceSheets/matlab/'])
addpath([input_drive,'/IceSheets/EM_Model_SnoWR_Algorithm/TOOLS/'])
addpath 'Z:\IceSheets\EM_Model_SnoWR_Algorithm\snow_mixing_formula'

freqGHz = 1.4; % GHz
temperature = 273.15; %K
density = 400; % kg/cm^3
lwc = linspace(0.000001,0.06,100); % volume fraction 0-1

TdryC = -0.5;
[eps_eff_dr] = diel_dry_snow_matzler_06(TdryC, density/1000, freqGHz );

eps_eff=nan(length(lwc),10);


for j = 1:length(lwc)
    eps_eff(j,1) = wetsnow_permittivity_Matzler(freqGHz, density, lwc(j));
    eps_eff(j,2) = wetsnow_permittivity_tinga73(freqGHz, density, lwc(j));
    eps_eff(j,3) = wetsnow_permittivity_debyelike_hallikainen86(freqGHz,density/1000, lwc(j)*100);
    eps_eff(j,4) = wetsnow_permittivity_hallikainen86(freqGHz,density/1000, lwc(j)*100);
    eps_eff(j,5) = wetsnow_permittivity_hallikainen86_ulaby14(freqGHz,density/1000, lwc(j)*100);
    eps_eff(j,6) = wetsnow_permittivity_colbeck80(freqGHz, density, lwc(j));
    eps_eff(j,7) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/2); % Birchak 1974
    eps_eff(j,8) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/2.5); % Sihvola (1985)
    eps_eff(j,9) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/3); %Looyenga (1965)
    eps_eff(j,10) = wetsnow_permittivity_tiuri_1984(freqGHz,density/1000, lwc(j));
end

%% Figure 2: Real and imaginary parts of complex dielectric constant
selected_models = 1:10;
lw=3;
lw2=2;
ms=20;
fs = 18;
figure(Position=[3.9126e+03 -420.6000 1.6388e+03 1.2264e+03])
% figure
t=tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

nexttile(1);
hold on
count=1;
for c=selected_models
    s = styleSet{count};
    plot(lwc*100,real(eps_eff(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw);
    count = count+1;
end
plot([0,lwc*100],ones(length(lwc)+1,1)*real(eps_eff_dr),'--',Color=[0.5,0.5,0.5],LineWidth=lw2, DisplayName= '');
grid;
grid minor
text(0.015,2.7,'(a)','FontWeight','bold','FontSize', fs);
ylabel('Real Part of Dielectric Conatant','FontSize',14,FontWeight='bold')
set(gca,'FontSize',14,FontWeight='bold')

nexttile(2)
hold on
count=1;
for c=selected_models
    s = styleSet{count};
    plot(lwc*100,real(eps_eff(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw);
    count = count+1;
end
plot([0,lwc*100],ones(length(lwc)+1,1)*real(eps_eff_dr),'--',Color=[0.5,0.5,0.5],LineWidth=lw2, DisplayName= '');
grid;
grid minor
xlim([0 2])
text(0.015,2.15,'(b)','FontWeight','bold','FontSize', fs);
set(gca,'FontSize',14,FontWeight='bold')

count=1;
nexttile(3)
hold on
for c=selected_models
    s = styleSet{count};
    plot(lwc*100,imag(eps_eff(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw);
    count = count+1;
end
h=plot([0,lwc*100],ones(length(lwc)+1,1)*imag(eps_eff_dr),'--',Color=[0.5,0.5,0.5],LineWidth=lw2, DisplayName= '');
grid;
grid minor
ylabel('Imiginary Part of Dielectric Conatant','FontSize',14,FontWeight='bold')
set(gca,'FontSize',14,FontWeight='bold');
h.HandleVisibility = 'off';
text(0.015,0.11,'(c)','FontWeight','bold','FontSize', fs);


count=1;
nexttile(4)
hold on
for c=selected_models
    s = styleSet{count};
    plot(lwc*100,imag(eps_eff(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw,DisplayName=char(models(c)));
    count = count+1;
end
h = plot([0, lwc*100], ones(length(lwc)+1,1)*imag(eps_eff_dr), '--', ...
    'Color', [0.5, 0.5, 0.5], 'LineWidth', lw2);
set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

xlim([0 2])
grid;
grid minor
set(gca,'FontSize',14,FontWeight='bold')
text(0.015,0.0275,'(d)','FontWeight','bold','FontSize', fs);

lgd = legend('location','northwest','FontSize',14,FontWeight='bold');
lgd.NumColumns = 3;
lgd.Position = [0.0900    0.8354    0.2520    0.1023];
legend('boxoff')
title(lgd,sprintf('Wet Snow Dielectric Models'))

xlabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')
set(gca,'FontSize',14,FontWeight='bold')

filename = sprintf('%sFig2_wetsnow_diel_10_models_v2.png',fig_dir);
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Fig. 3: Penetration depth as fn of LWC
deltap = penetration(freqGHz,eps_eff);
deltap_max = max(abs(deltap'));
deltap_min = min(abs(deltap'));
spread = deltap_max - deltap_min;
mv = lwc*100;
styleSet=plot_styles();
selected_models = 1:10;
lw=3;
lw2=2;
ms=20;
fs = 18;
% figure(Position=[4.1882e+03 -256.2000 1.1408e+03 793.2000]);hold on;
figure(Position=[4.1674e+03 -315.8000 1.4544e+03 1008]);hold on;
t=tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
nexttile
hold on
count=1;
for p=selected_models
    s = styleSet{count};
    plot(lwc*100,abs(deltap(:,p)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw,DisplayName=sprintf('%s  ',char(models(p))));
    count = count+1;
end
grid;
grid minor
ylim([0 20])
text(5.5,19,'(a)','FontWeight','bold','FontSize', fs);
set(gca,"FontWeight","bold",FontSize=14)

nexttile
hold on
count=1;
for p=selected_models
    s = styleSet{count};
    plot(lwc*100,abs(deltap(:,p)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw,DisplayName=sprintf('%s  ',char(models(p))));
    count = count+1;
end
ylim([0 5])
xlim([3 6])
grid;
grid minor
text(3.1,4.73,'(b)','FontWeight','bold','FontSize', fs);
set(gca,"FontWeight","bold",FontSize=14)
lgd = legend();
lgd.NumColumns = 3;
title(lgd,sprintf('Wet Snow Dielectric Models'))
legend('boxoff')

title(t,sprintf('L-band Penetration Depth (m) at Snow Density = %d kg/m^3',density),'fontweight','bold')
xlabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')
ylabel(t,'Penetration Depth (m)','FontSize',14,FontWeight='bold')

filename = sprintf('%sFig3_spenetration_10_models_v2.png',fig_dir);
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Load MEMLS SIMULATION LUT
clearvars -except  input_drive fig_dir models styleSet
result_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/results/'];
memls = load([result_dir,'LUT_memlsv3.mat']);
par = memls.par_range_LUT;
%% Fig. 4: V- and H-pol TBs as fn of mv for different twet
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);
rho = 0.4; %g/cm^3
% d2 = 100;
eps3r = 23;

p=1;
for d2 = [100:100:300,100:100:300]

    if p<4
        TB = memls.tbv;
        yl= [140,280];
    else
        TB = memls.tbh;
        yl=[120,260];
    end


    nt=nexttile(t)
    hold on

    count=1;
    lh=[];
    tb =[];
    for model = 1:10
        [~,i] = min(abs(par.Ps - rho));
        [~,j] = min(abs(par.d_dry2 - d2));
        [~,k ]= min(abs(par.eps3r - eps3r));

        tb1 = squeeze(TB(i,j,k,:,model));
        s = styleSet{count};
        h=plot(par.mv, tb1, '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, DisplayName=char(models(count)));
        lh=[lh,h];

        count = count +1;

    end
    %

    grid;grid minor
    ax1 = gca; % current axes
    xlim(ax1,[0,6])
    ylim(ax1,yl)
    xticks(ax1,0:1:6)
    ax2 = axes('Position',ax1.Position,...
        'XAxisLocation','top',...
        'YAxisLocation','left', ...
        'Color','none','XLim',ax1.XLim,'YLim',ax1.YLim) ;
    xlabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')

    xticks(ax2,par.mv(1:10:61))
    xticklabels(ax2,par.mv(1:10:61)*d2/10)
    ax2.TickDir = 'out';
    ax1.TickDir = 'out';
    % ax1.TickDir = 'none';
    yticklabels(ax2,{})
    % xticks(ax1,lon_tr(4:8:100))
    % xticklabels(ax1,round(dist(4:8:100)))

    set(ax1,"FontWeight","bold")
    set(ax2,"FontWeight","bold")
    set(ax1,'FontSize',14)
    set(ax2,'FontSize',14)

    % ax1.YAxis(2).Color = 'r';
    % ax2.YAxis(1).Color = 'k';
    if p==2
        xlabel(ax2, 'Liquid Water Amount (mm)')
    end
    if p==1
        ylabel(ax1,'V-pol Brightness Temperature (TBV) [K]','FontSize',14,FontWeight='bold')
    end
    if p==4
        ylabel(ax1,'H-pol Brightness Temperature (TBH) [K]','FontSize',14,FontWeight='bold')
    end
    text(2,yl(1)+10,sprintf('twet = %d cm',d2),'FontSize',14,FontWeight='bold')
    p=p+1;
end
lgd =legend(lh,'location','southeast','FontSize',14,FontWeight='bold');

lgd.NumColumns = 2;
% title(lgd,sprintf('Wet Snow Dielectric Models'))
lgd.Box = 'off';
lgd.Position = [0.74 0.66 0.18 0.1];

filename = sprintf('%sFig4_TB_VnH_with_diel_models.png',fig_dir);
exportgraphics(gcf, filename, 'Resolution', 300);
%% Fig. 5: V- and H-pol TBs as fn of mv for different rho
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);

% rho = 0.4; %g/cm^3
d2 = 200;
% eps3r = 23;

p=1;
for rho = [200:200:600,200:200:600]

    if p<4
        TB = memls.tbv;
        yl= [140,280];
    else
        TB = memls.tbh;
        yl=[125,265];
    end


    nt=nexttile(t)
    hold on

    count=1;
    lh=[];
    tb =[];
    for model = 1:10
        [~,i] = min(abs(par.Ps - rho/1000));
        [~,j] = min(abs(par.d_dry2 - d2));
        [~,k ]= min(abs(par.eps3r - eps3r));


        tb1 = squeeze(TB(i,j,k,:,model));

        s = styleSet{count};
        h=plot(par.mv, tb1,  '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2,  DisplayName=char(models(count)));
        lh=[lh,h];

        count = count +1;

    end
    %

    grid;grid minor
    ax1 = gca; % current axes
    xlim(ax1,[0,6])
    ylim(ax1,yl)
    xticks(ax1,0:1:6)
    ax2 = axes('Position',ax1.Position,...
        'XAxisLocation','top',...
        'YAxisLocation','left', ...
        'Color','none','XLim',ax1.XLim,'YLim',ax1.YLim) ;
    xlabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')

    xticks(ax2,par.mv(1:10:61))
    xticklabels(ax2,par.mv(1:10:61)*d2/10)
    ax2.TickDir = 'out';
    ax1.TickDir = 'out';
    % ax1.TickDir = 'none';
    yticklabels(ax2,{})
    % xticks(ax1,lon_tr(4:8:100))
    % xticklabels(ax1,round(dist(4:8:100)))

    set(ax1,"FontWeight","bold")
    set(ax2,"FontWeight","bold")
    set(ax1,'FontSize',14)
    set(ax2,'FontSize',14)

    % ax1.YAxis(2).Color = 'r';
    % ax2.YAxis(1).Color = 'k';
    if p==2
        xlabel(ax2, 'Liquid Water Amount (mm)')
    end
    if p==1
        ylabel(ax1,'V-pol Brightness Temperature (TBV) [K]','FontSize',14,FontWeight='bold')
    end
    if p==4
        ylabel(ax1,'H-pol Brightness Temperature (TBH) [K]','FontSize',14,FontWeight='bold')
    end
    text(2,yl(1)+10,sprintf('{\\rho_s} = %d kg/m^3',rho),'FontSize',14,FontWeight='bold')
    p=p+1;
end
lgd =legend(lh,'location','southeast','FontSize',14,FontWeight='bold');

lgd.NumColumns = 2;
% title(lgd,sprintf('Wet Snow Dielectric Models'))
lgd.Box = 'off';
lgd.Position = [0.74 0.66 0.18 0.1];

filename = sprintf('%sFig5_TB_VnH_fnofmv_diff_rho_twet300cm.png',fig_dir);
exportgraphics(gcf, filename, 'Resolution', 300);

%% Fig. 6: V- and H-pol TBs as fn of twet for different mv
[col]=distinc_colors();
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);

rho = 0.4; %g/cm^3
% d2 = 100;
% eps3r = 23;

p=1;
for mv = [1:3,1:3]

    if p<4
        TB = memls.tbv;
        yl= [140,280];
    else
        TB = memls.tbh;
        yl=[120,260];
    end


    nt=nexttile(t)
    hold on

    count=1;
    lh=[];
    tb =[];
    for model = 1:10
        [~,i] = min(abs(par.Ps - rho));
        % [~,j] = min(abs(par.d_dry2 - d2));
        [~,k ]= min(abs(par.eps3r - eps3r));
        [~,m ]= min(abs(par.mv - mv));


        tb2 = squeeze(TB(i,:,k,m,model));
        s = styleSet{count};
        h=plot(par.d_dry2, tb2, '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, DisplayName=char(models(count)));
        lh=[lh,h];

        count = count +1;

    end
    %

    grid;grid minor
    ax1 = gca; % current axes
    xlim(ax1,[0 600])
    ylim(ax1,yl)
    xticks(ax1,0:200:600)
    xticklabels(ax1,[0:200:600]./100)
    ax2 = axes('Position',ax1.Position,...
        'XAxisLocation','top',...
        'YAxisLocation','left', ...
        'Color','none','XLim',ax1.XLim,'YLim',ax1.YLim) ;
    xlabel(t,'Thickness of Wet Layer (m)','FontSize',14,FontWeight='bold')

    xticks(ax2,[0:200:600])
    xticklabels(ax2,mv*[0:200:600]/10)
    ax2.TickDir = 'out';
    ax1.TickDir = 'out';
    % ax1.TickDir = 'none';
    yticklabels(ax2,{})
    % xticks(ax1,lon_tr(4:8:100))
    % xticklabels(ax1,round(dist(4:8:100)))

    set(ax1,"FontWeight","bold")
    set(ax2,"FontWeight","bold")
    set(ax1,'FontSize',14)
    set(ax2,'FontSize',14)

    % ax1.YAxis(2).Color = 'r';
    % ax2.YAxis(1).Color = 'k';
    if p==2
        xlabel(ax2, 'Liquid Water Amount (mm)')
    end
    if p==1
        ylabel(ax1,'V-pol Brightness Temperature (TBV) [K]','FontSize',14,FontWeight='bold')
    end
    if p==4
        ylabel(ax1,'H-pol Brightness Temperature (TBH) [K]','FontSize',14,FontWeight='bold')
    end
    text(50,yl(1)+10,sprintf('m_v = %d %%',mv),'FontSize',14,FontWeight='bold')
    p=p+1;
end
lgd =legend(lh,'location','southeast','FontSize',14,FontWeight='bold');

lgd.NumColumns = 2;
% title(lgd,sprintf('Wet Snow Dielectric Models'))
lgd.Box = 'off';
lgd.Position = [0.74 0.66 0.18 0.1];

filename = sprintf('%sFig6_TB_VnH_fnof_twet.png',fig_dir);
% exportgraphics(gcf, filename, 'Resolution', 300);
%% Fig. 7: TB sensitivity of LWC
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(1,2);

rho = 0.4; %g/cm^3
% d2 = 100;
% eps3r = 23;
p=1;
xl= [0,150];

mv=3;
p=1;
for d2 = [mv,mv]

    if p == 1 || p == 3
        TB = memls.tbv;

    else
        TB = memls.tbh;

    end


    nt=nexttile(t)
    hold on

    count=1;
    lh=[];
    tb =[];
    TB_sen = [];
    for model = 1:10
        [~,i] = min(abs(par.Ps - rho));
        % [~,j] = min(abs(par.d_dry2 - d2));
        [~,k ]= min(abs(par.eps3r - eps3r));
        [~,m ]= min(abs(par.mv - mv));


        tb2 = squeeze(TB(i,:,k,m,model));
        tb = [tb,tb2];

        LWA = par.d_dry2*mv/10;
        delta_lwa=diff(LWA'); %mm
        delta_TB=diff(tb2); %K
        TB_sen(:,model) = delta_TB'./delta_lwa; % rate of change of TB per mm change in lWA

        s = styleSet{count};
        h=plot(LWA(2:end), TB_sen(:,model), '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, DisplayName=char(models(count)));

        lh=[lh,h];
        yline(0,'--')

        count = count +1;

    end
    %

    grid;grid minor
    ax1 = gca; % current axes
    xlim(xl)
    ylim([-0.5 6])
    % xticks(ax1,0:1:6)

    set(ax1,"FontWeight","bold",'FontSize',14)

    xlabel(t, 'Liquid Water Amount (mm)','FontSize',14,FontWeight='bold')
    ylabel(t,'Brightness Temperature Sensitivity [K/mm]','FontSize',14,FontWeight='bold')

    if p==1
        title(ax1,'V-pol Brightness Temperature','FontSize',14,FontWeight='bold')
    end
    if p==2
        title(ax1,'H-pol Brightness Temperature','FontSize',14,FontWeight='bold')
    end

    p=p+1;
end


lgd =legend(lh,'location','northeast','FontSize',14,FontWeight='bold');

lgd.NumColumns = 2;
title(lgd,sprintf('Wet Snow Dielectric Models'))
lgd.Box = 'off';
lgd.Position = [ 0.7229    0.7795    0.1708    0.1278];

filename = sprintf('%sFig7_TB_VnH_sensitivity_v2_mv_3parcent.png',fig_dir);
exportgraphics(gcf, filename, 'Resolution', 300);


%% Load TB and LWA Estimates
clearvars -except input_drive fig_dir models styleSet result_dir
aws_filename1='combined_corrected_promice_data_n_model_output_promice_tsurf1.mat';
if ~exist('aws','var')
    aws1=load(['Z:\IceSheets\Weather Stations Data\Greenland\PROMICE\DATA\AWS\V10\promice_n_model_data\V3\combined_multyyr\' ...
        aws_filename1]);
end
if ~exist('GEMB_compare_to_PROMICE_AWS','var')
    fpath_GEMBgrid2 = ['Z:/IceSheets/data/GEMB-AH/daily_LWC_Nicole/promice_forced/'];
    fnameo2 = 'GEMB_TotalmmLWCMelt_Greenland_promice_daily_2010-2025_revA.mat';
    gemb2=load([fpath_GEMBgrid2, fnameo2]);
end

smap_dir = [result_dir, '/lwa/v3/'];

%% Fig. 8: Simulated and obs TB time series
stations=aws1.stations;
[col]=distinc_colors();
figure('Position', [5.7362e+03 -385 1.9588e+03 1.2392e+03]);
t=tiledlayout(3,2,'TileSpacing','compact');
yl = [140 273.15];
% yl = [0 6];
yr1 = 2023;
mons = [5 11];
fs = 14;
Z=6;
count = 1;
TBmax = nan(6,10);
for stn = 1:length(stations)

    nexttile();hold on;
    clear smap_datac

    for model =  [1:10]
        smap = load([smap_dir,'smap_lwa_model_',num2str(model),'_avg_twet_avg_pstop3m.mat']);
        smap_data0 = smap.data_smap_ret{stn};
        smap_datac(model,:) = smap.data_smap_retc{stn};
        [smap_data,smap_dnum]=return_desired_data_revB(smap_data0,smap_data0(:,1),yr1,mons);
        tbv = smap_data(:,2);
        % tbh = smap_data(:,3);
        tbv_ref = smap_data(:,4);
        tbv_sim = smap_data(:,8);
        delta_tb = tbv - tbv_sim;
        TBmax(stn,model) = max(tbv_sim);

        s = styleSet{model};
        plot(smap_dnum, tbv_sim, '', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, 'MarkerSize', 16, 'DisplayName', char(models(model)));
    end
    plot(smap_dnum, tbv, 'b-.', 'LineWidth', 2, 'MarkerSize', 16, 'DisplayName', 'TBV obs');
    plot(smap_dnum, tbv_ref,'--','LineWidth', 2, 'Color', [.5 .5 .5], 'DisplayName', 'Frozen ref');
    xlim([datenum(yr1,6,1),datenum(yr1,11,1)]);
    grid
    % datetick('x','dd/mm/yy','keepticks');
    tickDates = datenum(yr1, 5:11, 1);        % 1st of each month
    tickDates15 = datenum(yr1, 5:10, 15);     % 15th of each month
    allTicks = sort([tickDates, tickDates15]);  % Combine and sort
    xticks(allTicks);
    xticklabels(datestr(allTicks, 'mmm dd'));

    title(sprintf(char(stations{(stn)})),  'FontSize', 18,'Interpreter', 'none')

    % ylim(yl)
    ax = gca;
    ax.XAxis.FontWeight = 'bold'; % Make X-axis tick labels bold
    ax.YAxis.FontWeight = 'bold'; % Make Y-axis tick labels bold
    ax.XAxis.FontSize = fs; % Make X-axis tick labels bold
    ax.YAxis.FontSize = fs; % Make Y-axis tick labels bold

    box on


end

ylabel(t,'Brightness Temperature [K]','FontSize',18,'FontWeight','bold')
xlabel(t,'Dates','FontSize',18,'FontWeight','bold')
% ylabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')
nexttile(1)
lgd=legend('location','Northeast','FontSize',12,'FontWeight','bold');
legend('boxoff')
lgd.NumColumns = 2;

filename = sprintf('%sFig.8_TB_.png',fig_dir,yr1);
% exportgraphics(gcf, filename, 'Resolution', 300);


%% Fig. 9: LWA time series
figure('Position', [5.7362e+03 -385 1.9588e+03 1.2392e+03]);
t=tiledlayout(3,2,'TileSpacing','compact');
yl = [-1 120];
% yl = [0 6];
yr1 = 2023;
z=aws1.z;

LWAmax = nan(6,12);
d2 = nan(6,10);
eps3r_est= nan(6,1);
rho3m= nan(6,1);
tbv_max= nan(6,1);
mdwet_samimi = nan(6,1);

lwa_smap_daily = [];

monset = nan(6,12);
fonset = nan(6,12);
mduration = nan(6,12);
maxsmelt = nan(6,12);

min_thresh = 0;
Z = 10;

mons = [5 10];
fs = 14;

md = nan(6,11);
sigma = nan(6,11);
mad = nan(6,11);
r = nan(6,11);
rmsd = nan(6,11);

mT_air= nan(6,1);
mT_airjfm= nan(6,1);
mT_airjja= nan(6,1);

for stn = 1:length(stations)

    % figure('Position', [5.9902e+03 -385 1.7048e+03 1150]);
    % t=tiledlayout(1,1,'TileSpacing','compact');

    % all available data
    aws_data = aws1.data_awsnm{stn};
    aws_dnum = aws_data(:,1);

    % yearly (all available in a year) results
    [aws_datay,aws_dnumy]=return_desired_data_revB(aws_data,aws_dnum,yr1,1:12);
    mT_air(stn) = mean(aws_datay(:,5));
    % aws_T = aws_data(:,6:48);
    % aws_rho = aws_data(:,92:134);
    % ind = find(z>3.0);
    % mrho = nanmean(aws_rho(:,1:ind(1)),2);
    % mT= nanmean(aws_T(:,1:ind(1)),2);
    % mrho(1)

    % Summer and winter results
    [aws_datajja,aws_dnumjja]=return_desired_data_revB(aws_data,aws_dnum,yr1,1:2);
    mT_airjfm(stn) = mean(aws_datajja(:,5));

    [aws_datajja,aws_dnumjja]=return_desired_data_revB(aws_data,aws_dnum,yr1,6:8);
    mT_airjja(stn) = mean(aws_datajja(:,5));

    % Summer melt results
    [aws_data,aws_dnum]=return_desired_data_revB(aws_data,aws_dnum,yr1,mons);



    lwa_aws = aws_data(:,49:91);
    dz = repmat([0.1,diff(z)],[size(lwa_aws,1),1]);
    dw=nan(size(lwa_aws));
    dw(lwa_aws>0)=dz(lwa_aws>0);
    dwet_samimi = nansum(dw,2);
    dwet_samimi(dwet_samimi==0) = nan; % remove meas when not melting
    mdwet_samimi(stn) = nanmean(dwet_samimi)*100; % in cm
    aws_lwa_hourly = sum(aws_data(:,49:91),2)*1000;
    [aws_days,aws_mwa_daily]=compute_daily_mean(aws_dnum,aws_lwa_hourly);
    nexttile();hold on;
    clear smap_datac
    lwa_smap = [];
    for model =  [1:10]

        smap = load([smap_dir,'smap_lwa_model_',num2str(model),'_avg_twet_avg_pstop3m.mat']);
        smap_data0 = smap.data_smap_ret{stn};
        smap_datac(model,:) = smap.data_smap_retc{stn};
        [smap_data,smap_dnum]=return_desired_data_revB(smap_data0,smap_data0(:,1),yr1,mons);
        tbv = smap_data(:,2);
        % tbh = smap_data(:,3);
        tbv_ref = smap_data(:,4);
        tbv_std = smap_data(:,5);
        tbv_thresh = tbv_ref+min(Z*tbv_std,10);
        meltflag = tbv>=(tbv_thresh);
        mv_smap = smap_data(:,6); % Melt volume fraction, Mv [%]
        mv_smap(mv_smap>10) = NaN;
        dwet_smap = smap_data(:,7); % thickness of the wet layer in cm
        d2(stn,model) = dwet_smap(1);%Table 3: Average thickness of wet layer (in cm) 			
        mwa_smap = (mv_smap.*meltflag/100).*dwet_smap*10; % *dwet_smap in cm, threfore mwa_smap in mm
        % meltF_dy = abs(tbv-tbv_ref)> tbv_thresh;
        [smap_days,mwa_smap_daily]=compute_daily_mean(smap_dnum,mwa_smap);
        LWAmax(stn,model) = max(mwa_smap_daily(:)); %Table 4: Maximum summer melt (in mm) during 2023 	

        eps3r_est(stn)= smap_datac(1,2);
        rho3m(stn)= smap_datac(1,4);
        tbv_max(stn)= max(tbv);
        s = styleSet{model};
        plot(smap_days, mwa_smap_daily, '', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', char(models(model)));
        [smap_days,lwa_smap_daily] = matchup_time_series(aws_mwa_daily,aws_days,mwa_smap_daily,smap_days);
        [monset(stn,model),fonset(stn,model),mduration(stn,model),maxsmelt(stn,model)] = melt_metrics_v2(smap_days,lwa_smap_daily,min_thresh,yr1);
        lwa_smap(:,model) = lwa_smap_daily;
    end

    lwa_smap(:,model+1) = nanmean(lwa_smap,2);

    gemb_lwa_daily = gemb2.gemb.LWC(:,stn);

    plot(gemb2.gemb.dnum,gemb_lwa_daily,':','LineWidth',3,  'Color', col(12), 'DisplayName', 'GEMB');
    plot(aws_days,aws_mwa_daily,'-.','LineWidth',3, 'Color', col(15), 'DisplayName', "SAMIMI");
    xlim([datenum(yr1,6,1),datenum(yr1,11,1)]);
    grid
    % datetick('x','dd/mm/yy','keepticks');
    tickDates = datenum(yr1, 5:11, 1);        % 1st of each month
    tickDates15 = datenum(yr1, 5:10, 15);     % 15th of each month
    allTicks = sort([tickDates, tickDates15]);  % Combine and sort
    xticks(allTicks);
    xticklabels(datestr(allTicks, 'mmm dd'));

    title(sprintf(char(stations{(stn)})),'FontSize',18, 'Interpreter', 'none')
    LWAmax(stn,11) = max(aws_mwa_daily(:));
    LWAmax(stn,12) = max(gemb_lwa_daily(:));

    [gemb_days,gemb_lwa] = matchup_time_series(aws_mwa_daily,aws_days,gemb_lwa_daily ,gemb2.gemb.dnum);

    lwa_ref = (gemb_lwa+aws_mwa_daily)/2;

    % plot(aws_days,lwa_ref,'-.','LineWidth',4, 'Color', col(13), 'DisplayName', "SAMIMI");

    ylim(yl)
    ax = gca;
    ax.XAxis.FontWeight = 'bold'; % Make X-axis tick labels bold
    ax.YAxis.FontWeight = 'bold'; % Make Y-axis tick labels bold
    ax.XAxis.FontSize = fs; % Make X-axis tick labels bold
    ax.YAxis.FontSize = fs; % Make Y-axis tick labels bold

    md(stn,:) = nanmean(lwa_ref-lwa_smap);
    sigma(stn,:)= nanstd(lwa_ref-lwa_smap);
    mad(stn,:) = nanmean(abs(lwa_ref-lwa_smap));
    % table 5, make sure that data is May-Sept
    r(stn,:) = corr(lwa_ref,lwa_smap,'rows','complete');
    % table 6
    rmsd(stn,:) = rmse(lwa_ref,lwa_smap,"omitnan");


    [monset(stn,11),fonset(stn,11),mduration(stn,11),maxsmelt(stn,11)] = melt_metrics_v2(aws_days,aws_mwa_daily,min_thresh,yr1);
    [monset(stn,12),fonset(stn,12),mduration(stn,12),maxsmelt(stn,12)] = melt_metrics_v2(gemb_days,gemb_lwa,min_thresh,yr1);

    dates_monset=datestr(monset(stn,:)','DD-mmm-yyyy');
    dates_fonset=datestr(fonset(stn,:)','DD-mmm-yyyy');
    box on

    metrics{stn}=find_the_metrics(lwa_smap,lwa_ref);
end

ylabel(t,'Liquid Water Amount [mm]','FontSize',18,'FontWeight','bold')
xlabel(t,'Dates','FontSize',18,'FontWeight','bold')
% ylabel(t,'Volumetric Liquid Water Content (%)','FontSize',14,FontWeight='bold')
nexttile(1)
lgd=legend('location','Northeast','FontSize',12,'FontWeight','bold');
legend('boxoff')
lgd.NumColumns = 2;


filename = sprintf('%sFig9_LWA_avg_twet_given_ps_given.png',fig_dir,yr1);
% exportgraphics(gcf, filename, 'Resolution', 300);

pairwise_rmsd_all = [];
for k = 1:6
    mask = triu(true(10),1);
    vals = metrics{k}.rmsd_pairwise(mask);
    pairwise_rmsd_all = [pairwise_rmsd_all; vals(~isnan(vals))];
end

mean_pairwise_rmsd = mean(pairwise_rmsd_all);
std_pairwise_rmsd = std(pairwise_rmsd_all);
fprintf('Mean pairwise RMSD across AWS: %.3f ± %.3f\n', mean_pairwise_rmsd, std_pairwise_rmsd);

%%
function metrics = find_the_metrics(lwa, ref_lwa)
% FIND_THE_METRICS computes pairwise and reference comparison metrics
%
% INPUTS:
%   lwa      - [T x 11] matrix. Columns 1–10: model outputs; column 11: ensemble mean.
%   ref_lwa  - [T x 1] reference LWA (e.g., from SEB model)
%
% OUTPUT:
%   metrics - struct with fields:
%       rmsd_pairwise [10x10]
%       mad_pairwise  [10x10]
%       r_pairwise    [10x10]
%       rmsd_to_ref   [11x1]
%       mad_to_ref    [11x1]
%       r_to_ref      [11x1]

% ------------------------------
% Setup
% ------------------------------
n_models = 10;  % First 10 columns are models
T = size(lwa, 1);

% Preallocate pairwise comparison matrices
rmsd_pairwise = nan(n_models, n_models);
mad_pairwise  = nan(n_models, n_models);
r_pairwise    = nan(n_models, n_models);

% Preallocate model vs. reference vectors (includes ensemble)
rmsd_to_ref = nan(n_models+1, 1);  % 1–10 = models, 11 = ensemble
mad_to_ref  = nan(n_models+1, 1);
r_to_ref    = nan(n_models+1, 1);

% ------------------------------
% Pairwise comparisons (models 1–10)
% ------------------------------
min_overlap = 3;  % minimum required overlapping non-NaN points

for i = 1:n_models
    for j = i+1:n_models
        xi = lwa(:,i);
        xj = lwa(:,j);
        valid = ~isnan(xi) & ~isnan(xj);
        if nnz(valid) >= min_overlap
            diff = xi(valid) - xj(valid);
            rmsd_pairwise(i,j) = sqrt(mean(diff.^2));
            mad_pairwise(i,j)  = mean(abs(diff));
            r_pairwise(i,j)    = corr(xi(valid), xj(valid));
        end
    end
end

% Symmetrize matrices
rmsd_pairwise = rmsd_pairwise + rmsd_pairwise';
mad_pairwise  = mad_pairwise  + mad_pairwise';
r_pairwise    = r_pairwise    + r_pairwise';

% ------------------------------
% Model and ensemble (col 11) vs. reference
% ------------------------------
for i = 1:(n_models + 1)
    xi = lwa(:,i);
    valid = ~isnan(xi) & ~isnan(ref_lwa);
    if nnz(valid) >= min_overlap
        diff = xi(valid) - ref_lwa(valid);
        rmsd_to_ref(i) = sqrt(mean(diff.^2));
        mad_to_ref(i)  = mean(abs(diff));
        r_to_ref(i)    = corr(xi(valid), ref_lwa(valid));
    end
end

% ------------------------------
% Pack outputs into struct
% ------------------------------
metrics.rmsd_pairwise = rmsd_pairwise;
metrics.mad_pairwise  = mad_pairwise;
metrics.r_pairwise    = r_pairwise;
metrics.rmsd_to_ref   = rmsd_to_ref;
metrics.mad_to_ref    = mad_to_ref;
metrics.r_to_ref      = r_to_ref;

end
