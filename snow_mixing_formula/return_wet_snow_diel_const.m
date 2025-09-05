function [eps_eff] = return_wet_snow_diel_const(freqGHz,density_gpcm3, mv_percent )

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


eps_eff(1) = wetsnow_permittivity_hallikainen86_ulaby14(freqGHz,density_gpcm3, mv_percent);
eps_eff(2) = wetsnow_permittivity_hallikainen86(freqGHz,density_gpcm3, mv_percent);
eps_eff(3) = wetsnow_permittivity_debyelike_hallikainen86(freqGHz,density_gpcm3, mv_percent);
eps_eff(4) = wetsnow_permittivity_memls(freqGHz, density_kgpm3, lwc);
eps_eff(5) = wetsnow_permittivity_wiesmann99(freqGHz, density_kgpm3, lwc);
eps_eff(6) = wetsnow_permittivity_tinga73(freqGHz, density_kgpm3, lwc);
eps_eff(7) = wetsnow_permittivity_colbeck80_caseI(freqGHz, density_kgpm3, lwc);
eps_eff(8) = wetsnow_permittivity_colbeck80_caseII(freqGHz, density_kgpm3, lwc);
eps_eff(9) = wetsnow_permittivity_colbeck80_caseIII(freqGHz, density_kgpm3, lwc);
eps_eff(10) = wetsnow_permittivity_three_component_polder_van_santen(freqGHz,density_kgpm3, lwc);

end