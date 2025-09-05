function [epsr,epsi] = return_wet_snow_diel_const_v2(freqGHz,density_gpcm3, mv_percent, modeln )

%Input Variables:
%freqGHz: frequency in GHz
%density_gpcm3: Dry Snow Density (g/cm^3)
%mv: volumetric water content (0<mv<6%)
%Output Products:
%eps:complex dielectric constant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
density_gpcm3 = double(density_gpcm3);
mv_percent = double(mv_percent);
density_kgpm3 = density_gpcm3*1000;
lwc = mv_percent/100;

switch modeln

    case 1
        eps_eff = wetsnow_permittivity_hallikainen86_ulaby14(freqGHz,density_gpcm3, mv_percent);
    case 2
        eps_eff = wetsnow_permittivity_hallikainen86(freqGHz,density_gpcm3, mv_percent);
    case 3
        eps_eff = wetsnow_permittivity_debyelike_hallikainen86(freqGHz,density_gpcm3, mv_percent);
    case 4
        eps_eff = wetsnow_permittivity_memls(freqGHz, density_kgpm3, lwc);
    case 5
        eps_eff = wetsnow_permittivity_wiesmann99(freqGHz, density_kgpm3, lwc);
    case 6
        eps_eff = wetsnow_permittivity_tinga73(freqGHz, density_kgpm3, lwc);
    case 7
        eps_eff = wetsnow_permittivity_colbeck80_caseI(freqGHz, density_kgpm3, lwc);
    case 8
        eps_eff = wetsnow_permittivity_colbeck80_caseIII(freqGHz, density_kgpm3, lwc);
    case 9
        eps_eff = wetsnow_permittivity_tiuri_1984(freqGHz,density_gpcm3, lwc);
    case 10
        eps_eff = wetsnow_permittivity_colbeck80_caseII(freqGHz, density_kgpm3, lwc);
    case 11
        eps_eff = wetsnow_permittivity_MG(freqGHz, density_kgpm3, lwc);
    case 12
        eps_eff = wetsnow_permittivity_three_component_polder_van_santen(freqGHz, density_kgpm3, lwc);
    case 13
        eps_eff = wetsnow_permittivity_power_law(freqGHz, density_kgpm3, lwc,1/2);
    case 14
        eps_eff = wetsnow_permittivity_power_law(freqGHz, density_kgpm3, lwc,1/3);
    case 15
        eps_eff = wetsnow_permittivity_power_law(freqGHz, density_kgpm3, lwc,0);
end

epsr = real(eps_eff);
epsi = imag(eps_eff);


end