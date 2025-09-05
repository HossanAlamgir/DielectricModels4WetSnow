
%% MEMLS simulation of Greenland Percolation Zone
%%% Simulation of brightness temperatures using MEMLSV3 in a simplified three-layer
% ice sheet configuration to represent equivalent snow and firn stratigraphy
% Contact: Alamgir Hossan, JPL, alamgir.hossan@jpl.nasa.gov, last updated 8/26/2025
% clear
cluster=1;
%Add the paths to MEMLS and snow mixing models
% The original version of MEMLS is available on GitHub at: https://github.com/akasurak/memls_TVC
% You’ll also find the document MEMLS3.pdf there, which includes a good user guide.
if cluster==0
    input_drive = 'Z:';
    output_drive = 'Z:';
    addpath ([input_drive,'/MEMLS3&a/AH'])    
    addpath([input_drive,'/IceSheets/EM_Model_SnoWR_Algorithm/snow_mixing_formula/'])
elseif cluster==1
    input_drive = '/projects/lotus';
    output_drive = '/projects/lotus';
    addpath('/projects/lotus/IceSheets/EM_Model_SnoWR_Algorithm/snow_mixing_formula/')    
    addpath('/projects/lotus/MEMLS3&a/AH')
end

% Sensor (SMPA) Obs parameters
fGHz=1.41; %GHz
fGHz        = [1.41];%GHz
thetad       = [40];

Tair = 260;% air temp (not used)

% Snow/firn layer parameters
Ps=single(linspace(0.2,0.9,36)); % firn density
d_dry2=single([10:10:100,120:20:300,340:40:500,600:100:2000]);
mv=[single([linspace(0,7,71)])];
eps3r = single(linspace(1,30,30));% Var real part of diel const of high-reflective layer
T0_dry2=250;

% Reflective layer parameters
% T0_dry3=T0_dry2;%use snow layer temp direcly or with a gradient
eps3i = 0.0002;% fixed Im part of diel const of high-reflective layer
d3 = 500; % thickness of the reflective layer, fixed (cm)

% Semi infinite bottom layer
% T0_dry4=255;
% T0_dry4=200;
T0_dry4 = 255;%single([linspace(200,255,12)]);
Pi = 0.917;% ice density

%Save LUT parameters
par_range_LUT.Ps= Ps;
par_range_LUT.d_dry2= d_dry2;
par_range_LUT.eps3r= eps3r;
par_range_LUT.mv= mv;
par_range_LUT.T0_dry2= T0_dry2;
par_range_LUT.d3= d3;
par_range_LUT.eps3i= eps3i;
par_range_LUT.T0_dry4= T0_dry4;

Nmodels = 10;

tbv = nan(length(Ps),length(d_dry2),length(eps3r),length(mv),Nmodels,length(fGHz));
tbh = nan(length(Ps),length(d_dry2),length(eps3r),length(mv),Nmodels,length(fGHz));

% Run MEMLS

for i=1:length(Ps)
    i
    for j = 1:length(d_dry2)

        for k = 1:length(eps3r)

            for m = 1:length(mv)

                for n = 1:Nmodels
                    if mv(m) == 0
                        T0_dry2=par_range_LUT.T0_dry2;
                        T0_dry3=T0_dry2;
                        [eps2r ,eps2i] = diel_dry_snow(T0_dry2-273, Ps(i), fGHz );
                    elseif mv(m)>0
                        T0_dry2 = 273.15;
                        T0_dry3=265;
                        [eps2r ,eps2i] = return_wet_snow_diel_const_v4(fGHz, Ps(i), mv(m),n);                       
                    end

                    eps2=eps2r-1i*eps2i;
                    eps3 = eps3r(k)-1i*eps3i;
                    [eps4r ,eps4i] = diel_dry_snow(T0_dry4-273, Pi, fGHz );
                    eps4=eps4r-1i*eps4i;                                     
                    eps = [1 ,eps2,eps3,eps4];                                    

                    [tbh(i,j,k,m,n),tbv(i,j,k,m,n)] =return_memls_TB_v2(fGHz,thetad,[d_dry2(j),d3,100000]',[T0_dry2,T0_dry3,T0_dry4]',[1000*Ps(i),1000*Ps(i),917]',[mv(m)/100,0,0]',eps(2:4)');
                end
            end
        end
    end
end


%% Save the LUT in matlab file
save([input_drive,'/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/results/LUT_memlsv3_v2.mat'],'tbv', 'tbh',...
    'par_range_LUT','-v7.3');

%% Save the LUT in netcdf
% Define output NetCDF file path
nc_filename = [input_drive, '/IceSheets/Documentation/Manuscripts_n_abstracts/wet_snow_diel_models/results/LUT_memlsv3_v2.nc'];

% Create NetCDF file
nccreate(nc_filename, 'tbv', ...
    'Dimensions', {'Ps', length(par_range_LUT.Ps), 'd_dry', length(par_range_LUT.d_dry2), 'eps2r', length(par_range_LUT.eps3r), 'mv', length(par_range_LUT.mv), 'model', Nmodels, 'freq', length(fGHz)}, ...
    'Datatype', 'single');

nccreate(nc_filename, 'tbh', ...
    'Dimensions', {'Ps', length(par_range_LUT.Ps), 'd_dry', length(par_range_LUT.d_dry2), 'eps2r', length(par_range_LUT.eps3r), 'mv', length(par_range_LUT.mv), 'model', Nmodels, 'freq', length(fGHz)}, ...
    'Datatype', 'single');

% Parameter dimensions
nccreate(nc_filename, 'Ps', 'Dimensions', {'Ps', length(Ps)}, 'Datatype', 'single');
nccreate(nc_filename, 'd_dry', 'Dimensions', {'d_dry', length(d_dry2)}, 'Datatype', 'single');
nccreate(nc_filename, 'eps2r', 'Dimensions', {'eps2r', length(eps3r)}, 'Datatype', 'single');
nccreate(nc_filename, 'mv', 'Dimensions', {'mv', length(mv)}, 'Datatype', 'single');
nccreate(nc_filename, 'fGHz', 'Dimensions', {'freq', length(fGHz)}, 'Datatype', 'single');

% Static parameters
nccreate(nc_filename, 'T0_dry', 'Dimensions', {}, 'Datatype', 'single');
nccreate(nc_filename, 'd2', 'Dimensions', {}, 'Datatype', 'single');
nccreate(nc_filename, 'eps2i', 'Dimensions', {}, 'Datatype', 'single');
nccreate(nc_filename, 'T0_dry3', 'Dimensions', {}, 'Datatype', 'single');

% Write data
ncwrite(nc_filename, 'tbv', single(tbv));
ncwrite(nc_filename, 'tbh', single(tbh));

ncwrite(nc_filename, 'Ps', Ps);
ncwrite(nc_filename, 'd_dry', d_dry2);
ncwrite(nc_filename, 'eps2r', eps3r);
ncwrite(nc_filename, 'mv', mv);
ncwrite(nc_filename, 'fGHz', fGHz);

ncwrite(nc_filename, 'T0_dry', T0_dry2);
ncwrite(nc_filename, 'd2', d3);
ncwrite(nc_filename, 'eps2i', eps3i);
ncwrite(nc_filename, 'T0_dry3', T0_dry4);

% Add units and metadata
ncwriteatt(nc_filename, 'tbv', 'units', 'K');
ncwriteatt(nc_filename, 'tbv', 'long_name', 'Brightness Temperature V-pol');

ncwriteatt(nc_filename, 'tbh', 'units', 'K');
ncwriteatt(nc_filename, 'tbh', 'long_name', 'Brightness Temperature H-pol');

ncwriteatt(nc_filename, 'Ps', 'units', 'g/cm^3');
ncwriteatt(nc_filename, 'Ps', 'long_name', 'Firn density');

ncwriteatt(nc_filename, 'd_dry', 'units', 'cm');
ncwriteatt(nc_filename, 'd_dry', 'long_name', 'Snow layer thickness');

ncwriteatt(nc_filename, 'eps2r', 'units', 'unitless');
ncwriteatt(nc_filename, 'eps2r', 'long_name', 'Real part of dielectric constant of the intermediate reflective layer');

ncwriteatt(nc_filename, 'mv', 'units', 'fractional vol. %');
ncwriteatt(nc_filename, 'mv', 'long_name', 'Liquid water content');

ncwriteatt(nc_filename, 'fGHz', 'units', 'GHz');
ncwriteatt(nc_filename, 'fGHz', 'long_name', 'Frequency');

ncwriteatt(nc_filename, 'T0_dry', 'units', 'K');
ncwriteatt(nc_filename, 'T0_dry', 'long_name', 'Snow layer temperature');

ncwriteatt(nc_filename, 'd2', 'units', 'cm');
ncwriteatt(nc_filename, 'd2', 'long_name', 'Thickness of reflective layer');

ncwriteatt(nc_filename, 'eps2i', 'units', 'unitless');
ncwriteatt(nc_filename, 'eps2i', 'long_name', 'Imaginary part of dielectric constant of the intermediate reflective layer');

ncwriteatt(nc_filename, 'T0_dry3', 'units', 'K');
ncwriteatt(nc_filename, 'T0_dry3', 'long_name', 'Temperature of semi-infinite bottom layer');

% Add global attributes
ncwriteatt(nc_filename, '/', 'title', 'MEMLS Simulation of Greenland Percolation Zone');
ncwriteatt(nc_filename, '/', 'description', 'Brightness temperatures from MEMLS V3 using three-layer firn stratigraphy');
ncwriteatt(nc_filename, '/', 'institution', 'JPL');
ncwriteatt(nc_filename, '/', 'author', 'Alamgir Hossan');
ncwriteatt(nc_filename, '/', 'contact', 'alamgir.hossan@jpl.nasa.gov');
ncwriteatt(nc_filename, '/', 'date_created', datestr(now, 'yyyy-mm-dd'));
