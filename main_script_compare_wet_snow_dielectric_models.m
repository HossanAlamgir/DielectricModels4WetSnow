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

% Make sure you downloaded the data/results from
% https://zenodo.org/uploads/17195725 and add the path in final_result_dir
fig_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/figs/revised/'];
code_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/codes/'];
final_result_dir = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/final_results/'];
% Add the 'snow_mixing_formula' subfolder to the MATLAB path
addpath(fullfile(code_dir, 'snow_mixing_formula'));
models = [{'Mätzler'},{'Tinga'},{'Debye-like'},{'Hallikainen'},{'Ulaby'},...
    {'Colbeck'},{'Birchak'},{'Sihvola'},{'Looyenga'},{'Tiuri'}];
styleSet=plot_styles();
%% Compute effective permitivity and the depth of penetration
filename = fullfile(final_result_dir, 'Fig2n3_effective_diel_const_n_pen_depth.nc');
if isfile(filename)
    disp('File already exixts.')
else
    freqGHz = 1.4; % Frequency in GHz for the simulation
    temperature = 273.15; % Melting point temperature in K
    density = 400; % dry snow density kg/m^3 for the simulation
    lwc = linspace(0.000001,0.06,100); % volume fraction 0-6%

    TdryC = -0.5; % dry snow temperature in degree C for the simulation
    [eps_eff_dr] = diel_dry_snow_matzler_06(TdryC, density/1000, freqGHz );

    eps_eff_wet=nan(length(lwc),10);


    for j = 1:length(lwc)
        eps_eff_wet(j,1) = wetsnow_permittivity_Matzler(freqGHz, density, lwc(j));
        eps_eff_wet(j,2) = wetsnow_permittivity_tinga73(freqGHz, density, lwc(j));
        eps_eff_wet(j,3) = wetsnow_permittivity_debyelike_hallikainen86(freqGHz,density/1000, lwc(j)*100);
        eps_eff_wet(j,4) = wetsnow_permittivity_hallikainen86(freqGHz,density/1000, lwc(j)*100);
        eps_eff_wet(j,5) = wetsnow_permittivity_hallikainen86_ulaby14(freqGHz,density/1000, lwc(j)*100);
        eps_eff_wet(j,6) = wetsnow_permittivity_colbeck80(freqGHz, density, lwc(j));
        eps_eff_wet(j,7) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/2); % Birchak 1974
        eps_eff_wet(j,8) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/2.5); % Sihvola (1985)
        eps_eff_wet(j,9) = wetsnow_permittivity_power_law(freqGHz, density, lwc(j),1/3); %Looyenga (1965)
        eps_eff_wet(j,10) = wetsnow_permittivity_tiuri_1984(freqGHz,density/1000, lwc(j));
    end

    deltap = penetration(freqGHz,eps_eff_wet);

    % --- Wet snow effective permittivity ---
    nccreate(filename, 'eps_eff_wet_real', ...
        'Dimensions', {'lwc', length(lwc), 'model', length(models)}, ...
        'Datatype', 'double');
    ncwrite(filename, 'eps_eff_wet_real', real(eps_eff_wet));

    nccreate(filename, 'eps_eff_wet_imag', ...
        'Dimensions', {'lwc', length(lwc), 'model', length(models)}, ...
        'Datatype', 'double');
    ncwrite(filename, 'eps_eff_wet_imag', imag(eps_eff_wet));

    % --- Penetration depth ---
    nccreate(filename, 'penetration_depth', ...
        'Dimensions', {'lwc', length(lwc), 'model', length(models)}, ...
        'Datatype', 'double');
    ncwrite(filename, 'penetration_depth', deltap);
    ncwriteatt(filename, 'penetration_depth', 'units', 'm');

    % --- Liquid water content/volume fraction (mv) ---
    nccreate(filename, 'lwc', 'Dimensions', {'lwc', length(lwc)});
    ncwrite(filename, 'lwc', lwc);

    % --- Dry snow effective permittivity (single complex number) ---
    nccreate(filename, 'eps_eff_dry_real', 'Datatype', 'double');
    ncwrite(filename, 'eps_eff_dry_real', real(eps_eff_dr));

    nccreate(filename, 'eps_eff_dry_imag', 'Datatype', 'double');
    ncwrite(filename, 'eps_eff_dry_imag', imag(eps_eff_dr));

   % Store simulation constants
    
    nccreate(filename, 'freqGHz', 'Datatype', 'double');
    ncwrite(filename, 'freqGHz', freqGHz);

    nccreate(filename, 'melting_point_temperature_K', 'Datatype', 'double');
    ncwrite(filename, 'melting_point_temperature_K', temperature);

    nccreate(filename, 'dry_snow_density_kg_per_m3', 'Datatype', 'double');
    ncwrite(filename, 'dry_snow_density_kg_per_m3', density);

    nccreate(filename, 'dry_snow_temperature_C', 'Datatype', 'double');
    ncwrite(filename, 'dry_snow_temperature_C', TdryC);

    model_string = strjoin(models, ', ');
    ncwriteatt(filename, '/', 'models', model_string);
    ncwriteatt(filename, '/', 'freqGHz', ...
        'Simulation frequency in GHz');
    ncwriteatt(filename, '/', 'melting_point_temperature_K', ...
        'Melting point temperature in Kelvin');
    ncwriteatt(filename, '/', 'dry_snow_density_kg_per_m3', ...
        'Dry snow density in kg/m^3 ');
    ncwriteatt(filename, '/', 'dry_snow_temperature_C', ...
        'Dry snow temperature in degrees Celsius');
    ncwriteatt(filename, '/', 'description', ...
        'Effective dielectric constant (complex) of wet and dry snow; 10 models over 100 LWC values');

end
%% Figure 2: Real and imaginary parts of complex dielectric constant
if ~exist("eps_eff",'var')
    filename = fullfile(final_result_dir, 'Fig2n3_effective_diel_const_n_pen_depth.nc');
    epswr = ncread(filename,'eps_eff_wet_real');
    epswi = ncread(filename,'eps_eff_wet_imag');
    epsr_dr = ncread(filename,'eps_eff_dry_real');
    epsi_dr = ncread(filename,'eps_eff_dry_imag');
    lwc = ncread(filename,'lwc')';
    eps_eff_wet = epswr+1i*epswi;
    eps_eff_dr  = epsr_dr+1i*epsi_dr;
end
selected_models = 1:10;
lw1=3;
lw2=2;
ms=20;
fs = 18;
styleSet=plot_styles();
figure(Position=[3.9126e+03 -420.6000 1.6388e+03 1.2264e+03])
t=tiledlayout(2,2,"TileSpacing","compact","Padding","compact");

nexttile(1);
hold on
count=1;
for c=selected_models
    s = styleSet{count};
    plot(lwc*100,real(eps_eff_wet(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1);
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
    plot(lwc*100,real(eps_eff_wet(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1);
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
    plot(lwc*100,imag(eps_eff_wet(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1);
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
    plot(lwc*100,imag(eps_eff_wet(:,c)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1,DisplayName=char(models(c)));
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
if ~exist("deltap",'var')
    filename = fullfile(final_result_dir, 'Fig2n3_effective_diel_const_n_pen_depth.nc');
    deltap = ncread(filename,'penetration_depth');
    density = ncread(filename,'dry_snow_density_kg_per_m3');
end

deltap_max = max(abs(deltap'));
deltap_min = min(abs(deltap'));
spread = deltap_max - deltap_min;
mv = lwc*100;
selected_models = 1:10;
lw1=3;
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
    plot(lwc*100,abs(deltap(:,p)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1,DisplayName=sprintf('%s  ',char(models(p))));
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
    plot(lwc*100,abs(deltap(:,p)),'.', 'Color', s.Color, 'LineStyle', s.LineStyle,LineWidth=lw1,DisplayName=sprintf('%s  ',char(models(p))));
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
clearvars -except  input_drive fig_dir models styleSet final_result_dir
LUT = fullfile([final_result_dir,'LUT_memlsv4.nc']);
memls.tbv = ncread(LUT,'tbv');
memls.tbh = ncread(LUT,'tbh');
par.rhos = ncread(LUT,'rhos');
par.d_dry = ncread(LUT,'d_dry');
par.eps2r = ncread(LUT,'eps2r');
par.mv = ncread(LUT,'mv');
%% Fig. 4: V- and H-pol TBs as fn of mv for different twet
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);
rho = 0.4; %g/cm^3
eps2r = 23;
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
        [~,i] = min(abs(par.rhos - rho));
        [~,j] = min(abs(par.d_dry - d2));
        [~,k ]= min(abs(par.eps2r - eps2r));

        tb1 = squeeze(TB(i,j,k,:,model));
        s = styleSet{count};
        h=plot(par.mv, tb1, '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, DisplayName=char(models(count)));
        lh=[lh,h];

        count = count +1;

    end

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
% exportgraphics(gcf, filename, 'Resolution', 300);
%% Fig. 5: V- and H-pol TBs as fn of mv for different rho
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);
d2 = 200;
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
        [~,i] = min(abs(par.rhos - rho/1000));
        [~,j] = min(abs(par.d_dry - d2));
        [~,k ]= min(abs(par.eps2r - eps2r));

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
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Fig. 6: V- and H-pol TBs as fn of twet for different mv
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(2,3);

rho = 0.4; %g/cm^3
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
        [~,i] = min(abs(par.rhos - rho));
        [~,k ]= min(abs(par.eps2r - eps2r));
        [~,m ]= min(abs(par.mv - mv));

        tb2 = squeeze(TB(i,:,k,m,model));
        s = styleSet{count};
        h=plot(par.d_dry, tb2, '.', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, DisplayName=char(models(count)));
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
        [~,i] = min(abs(par.rhos - rho));
        [~,k ]= min(abs(par.eps2r - eps2r));
        [~,m ]= min(abs(par.mv - mv));

        tb2 = squeeze(TB(i,:,k,m,model));
        tb = [tb,tb2];

        LWA = par.d_dry*mv/10;
        delta_lwa=diff(LWA'); %mm
        delta_TB=diff(tb2); %K
        TB_sen(:,model) = delta_TB./delta_lwa; % rate of change of TB per mm change in lWA

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
% exportgraphics(gcf, filename, 'Resolution', 300);

%% Fig. 8: Simulated and obs TB time series
stations=readtable(sprintf('%s/AWS.xlsx',final_result_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
figure(Position=[4.0286e+03 -485.8000 1.6932e+03 1.1676e+03]);
t=tiledlayout(3,2,'TileSpacing','compact');
yl = [140 273.15];
yr1 = 2023;
mons = [5 11];
fs = 14;
Z=6;
count = 1;
TBmax = nan(6,10);
for stn = 1:length(stations)
    filename = fullfile(final_result_dir, sprintf('SMAP_TB_n_LWA_at_%s.nc', stations{(stn)}));
    smap_dnum = ncread(filename,'time');
    tbv = ncread(filename,'tbv_smap_obs');
    % tbh = ncread(filename,'tbh_smap_obs');
    tbv_ref = ncread(filename,'tbv_frozen_mean');
    tbv_sim = ncread(filename,'tbv_model');
    nexttile();hold on;
    clear smap_datac

    for model =  [1:10]

        delta_tb = tbv - tbv_sim(:,model);
        TBmax(stn,model) = max(tbv_sim(:,model));
        s = styleSet{model};
        plot(smap_dnum, tbv_sim(:,model), '', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, 'MarkerSize', 16, 'DisplayName', char(models(model)));
    end
    plot(smap_dnum, tbv, 'b-.', 'LineWidth', 2, 'MarkerSize', 16, 'DisplayName', 'TBV obs');
    plot(smap_dnum, tbv_ref*ones(length(smap_dnum),1),'--','LineWidth', 2, 'Color', [.5 .5 .5], 'DisplayName', 'Frozen ref');
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
clearvars -except input_drive fig_dir models styleSet result_dir final_result_dir
stations=readtable(sprintf('%s/AWS.xlsx',final_result_dir), 'ReadVariableNames', false);
stations = table2cell(stations);
[col]=distinc_colors();
figure('Position', [3.1802e+03 31 1.7292e+03 838.8000]);
t=tiledlayout(3,2,'TileSpacing','compact');
yl = [-1 120];
% yl = [0 6];
yr1 = 2023;
LWAmax = nan(6,12);
d2 = nan(6,10);
eps2r_est= nan(6,1);
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
    aws = stations{stn}; % station name
    aws_filename = fullfile(final_result_dir, sprintf('AWS_n_SAMIMI_model_data_at_%s.nc', aws));
    aws_dnum = ncread(aws_filename,'time');
    aws_lwc = ncread(aws_filename,'mv_emb');
    aws_lwa = ncread(aws_filename,'LWA_emb');
    aws_tair = ncread(aws_filename,'Tair');
    model_depth_z_samimi = ncread(aws_filename,'model_depth');
    aws_data = [aws_dnum,aws_lwa,aws_tair,aws_lwc];


    % 2023 results
    [aws_datay,aws_dnumy]=return_desired_data_revB(aws_data,aws_dnum,yr1,1:12);
    mT_air(stn) = mean(aws_datay(:,3));

    % Summer and winter results
    [aws_datajja,aws_dnumjja]=return_desired_data_revB(aws_data,aws_dnum,yr1,1:2);
    mT_airjfm(stn) = mean(aws_datajja(:,3));

    [aws_datajja,aws_dnumjja]=return_desired_data_revB(aws_data,aws_dnum,yr1,6:8);
    mT_airjja(stn) = mean(aws_datajja(:,3));

    % Summer melt results
    [aws_data,aws_dnum]=return_desired_data_revB(aws_data,aws_dnum,yr1,mons);

    lwc_aws = aws_data(:,4:46);
    dz = repmat([0.1,diff(model_depth_z_samimi')],[size(lwc_aws,1),1]);
    dw=nan(size(lwc_aws));
    dw(lwc_aws>0)=dz(lwc_aws>0);
    dwet_samimi = nansum(dw,2);
    dwet_samimi(dwet_samimi==0) = nan; % remove meas when not melting
    mdwet_samimi(stn) = nanmean(dwet_samimi)*100; % in cm
    aws_lwa_hourly = sum(lwc_aws,2)*1000;
    [aws_days,aws_mwa_daily]=compute_daily_mean(aws_dnum,aws_lwa_hourly);
    nexttile();hold on;
    clear smap_datac
    lwa_smap = [];
    for model =  [1:10]
        filename = fullfile(final_result_dir, sprintf('SMAP_TB_n_LWA_at_%s.nc', stations{(stn)}));
        smap_dnum = ncread(filename,'time');
        tbv = ncread(filename,'tbv_smap_obs');
        tbh = ncread(filename,'tbh_smap_obs');
        tbv_ref = ncread(filename,'tbv_frozen_mean');
        tbv_std = ncread(filename,'tbv_frozen_std');
        mv_smap = ncread(filename,'mv');
        mv_smap = mv_smap(:,model);
        dwet_smap = ncread(filename,'twet');
        dwet_smap = dwet_smap(:,model);
        eps2r_est(stn)= ncread(filename,'eps2r');
        rho3m(stn)= ncread(filename,'rho_top3m');
        smap_data0 = [smap_dnum, tbv, tbh, tbv_ref*ones(length(tbv),1), tbv_std*ones(length(tbv),1),mv_smap,dwet_smap];
        [smap_data,smap_dnum]=return_desired_data_revB(smap_data0,smap_data0(:,1),yr1,mons); % filter 2023 data
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

        tbv_max(stn)= max(tbv);
        s = styleSet{model};
        plot(smap_days, mwa_smap_daily, '', 'Color', s.Color, 'LineStyle', s.LineStyle, 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', char(models(model)));
        [smap_days,lwa_smap_daily] = matchup_time_series(aws_mwa_daily,aws_days,mwa_smap_daily,smap_days);
        [monset(stn,model),fonset(stn,model),mduration(stn,model),maxsmelt(stn,model)] = melt_metrics_v2(smap_days,lwa_smap_daily,min_thresh,yr1);
        lwa_smap(:,model) = lwa_smap_daily;
    end

    lwa_smap(:,model+1) = nanmean(lwa_smap,2);

    gemb_filename = fullfile(final_result_dir, sprintf('GEMB_model_data_at_%s.nc', aws));
    gemb_dnum = ncread(gemb_filename,'time');
    gemb_lwa_daily = ncread(gemb_filename,'LWA_daily');

    plot(gemb_dnum,gemb_lwa_daily,':','LineWidth',3,  'Color', col(12), 'DisplayName', 'GEMB');
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

    [gemb_days,gemb_lwa] = matchup_time_series(aws_mwa_daily,aws_days,gemb_lwa_daily ,gemb_dnum);

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

