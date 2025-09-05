function [Tbh,Tbv] =return_memls_TB_v2(fGHz,thetad,thickness,temp,rho,mv,eps_ref)
% NOTE: This code is modified from the main MEMLS main code named as memlsmain.
if size(thickness,2)>1 | size(temp,2)>1 | size(rho,2)>1 | size(mv,2)>1
    disp('Inputs must be colunm vectors')
    return
end

%%
n_days = 1;
n_layer = length(thickness);
temp_layer = flip(temp);
lwc_layer_volfrac = flip(mv);
rhos_layer = flip(rho);
d_layer = flip(thickness);
salin_layer = 0 + zeros(n_layer,1);
corr_len_layer = 0 + zeros(n_layer,1);

bottom = 0;
showmessage = 0;


%% Substrate

if bottom == 0 % ground

    Tgnd        = 270;  % Tgnd:  ground temperature [K]

    Ts           = Tgnd-273.15;
    Ts(Ts < -30) = -30;

    m0 = 35;     % content of organic matter in percentages
    mv = 0.30;   % soil bulk dry density g/cm^3
    rd = 0.8;    % volumetric soil moisture cm^3/cm^3
    mg = mv/rd;  % gravimetric soil moisture g/g


    es = GRMDM_O(m0, mg, rd, Ts);

    [rv_ground_ice, rh_ground_ice] = fresnel_reflectivity([real(es), imag(es)], 40, [3.2, 0]);

    s0h         = rh_ground_ice;                % ice-ground reflectivity (H-pol)
    s0v         = rv_ground_ice;                % ice-ground reflectivity (V-pol)

elseif bottom == 1 % ocean
    % To simulated OCEAN-ICE interface

    [rv_water_ice, rh_water_ice] = fresnel_reflectivity([80, 0], 40, [3.2, 0]);

    s0h         = rh_water_ice;            % ice-ocean reflectivity (H-pol)
    s0v         = rv_water_ice;            % ice-ocean reflectivity (V-pol)

    Tgnd        = temp_layer(1);      % Tgnd:  water temperature [K]

    %         wetness_layer(2) = 0.10;
    %         wetness_layer(1) = 0.15;

end
%%
Tsky        = 2.7;                %     Tsky:  sky brightness temperature [K]

sccho       = 11;               %     sccho: type of scattering coefficient (11 recommended)

%  Aux. input parameter for dielectric model of snow:
graintype   = 1;	% 1, 2 or 3 for different assumptions:
%"1" empirical snow measurements, "2" small spheres, "3" thin shells

c0          = 0.299793;               % vac. speed of light in m/ns

n_angles    = length(thetad);

for ii=1:length(fGHz)

    freq = fGHz(ii);

    theta = thetad(ii)/180*pi;

    [layer_char] = combine_daily_ground_truth(n_days, n_layer, temp_layer, lwc_layer_volfrac, rhos_layer,...
        d_layer, salin_layer, corr_len_layer);

    y           = layer_char;%(:,:,jj);

    di          = y(:,5);
    y1          = find(di>0);
    num         = y(y1,1);
    Ti          = y(:,2);
    Wi          = y(:,3);
    roi         = y(:,4);
    Sppt        = y(:,6);
    pci         = y(:,7);
    N           = length(num);
    if N == 0
        return                                      % test if there is there a layer at all
    end
    roi         = roi./1000;                        % transforms density to units of g/cm3
    di          = di./100;                          % transforms thicknesses to units of m
    cc          = vK2eps(roi, graintype);
    v           = cc(:,1);                          % ice volume fraction
    kq          = cc(:,2);                          % K^2 ratio
    epsid       = cc(:,3);                          % epsilon of dry snow
    nid         = sqrt(epsid);                      % real refract.index of dry snow
    eii         = epsaliceimag(freq, Ti, Sppt);       % imag epsilon of saline ice
    epsiid      = v.*eii.*kq.*nid;                  % imag epsilon of dry snow component
    epsd        = epsid+1i*epsiid;

    %% Insert the reflective layer
    eps = flip(eps_ref);
    % epsd(2) = eps_ref;
    %%

    % eps         = epswet(freq, Ti, Wi, epsd);
    epsi        = real(eps);                        % real epsilon of snow, dry or wet
    epsii       = imag(eps);                        % imag epsilon of snow, dry or wet

    % epsi_new(:,ii,jj)  = epsi;
    % epsii_new(:,ii,jj) = epsii;

    gai         = (4*pi*freq).*imag(sqrt(eps))./c0; % absorption coeff (1/m)

    ns          = sqrt(epsi);                       % approx. real refractive index of snow
    tei         = [asin(sin(theta)./ns);theta];
    dei         = pfadi(tei, di);

    [sih, siv]   = fresnelyc(tei, [epsi;1]);
    %         reflec_h(:,=sih

    [rnum, rroi, repsi, repsii, rtei, rsih, rsiv, rdi, rdei, rTi, rpci, rWi, rgai, rkq] ...
        = slred(num, roi, epsi, epsii, tei, sih, siv, di, dei, Ti, pci, freq, Wi, gai, kq, showmessage);

    [gbih,gbiv,gs6,ga2i] ...
        = sccoeff(rroi, rTi, rpci, freq, repsi, rgai, sccho, rkq);
    %         ga2i_plot(:,ii,jj)=ga2i;
    %         gbih_new(:,ii,jj)=gbih;
    %         gbiv_new(:,ii,jj)=gbiv;
    %         gs6_new(:,ii,jj)=gbiv;
    %
    [rdei, rtei, tscat] ...
        = pfadc(theta, rdi, repsi, gs6);
    rsih        = [s0h; rsih];
    rsiv        = [s0v; rsiv];
    [rsih,rsiv] = polmix(tscat, rsih, rsiv);



    N                    = length(rroi);

    rsih_top(ii) = rsih(N+1);
    rsiv_top(ii) = rsiv(N+1);

    [rih, tih]       = rt(ga2i, gbih, rdei);
    Dh               = layer(rih, rsih, tih, rTi, Tgnd, Tsky);
    Tbh(ii)  = (1-rsih(N+1))*Dh(N) + rsih(N+1)*Tsky;

    [riv, tiv]       = rt(ga2i, gbiv, rdei);
    Dv               = layer(riv, rsiv, tiv, rTi, Tgnd, Tsky);
    Tbv(ii)  = (1-rsiv(N+1))*Dv(N) + rsiv(N+1)*Tsky;

end

end
