
function [eps_eff]=wetsnow_permittivity_colbeck80(freqGHz, density, lwc)

if density<=550 && lwc <=0.07
    eps_eff = wetsnow_permittivity_colbeck80_caseI(freqGHz, density, lwc);
elseif density>550 && lwc <=0.07
    eps_eff = wetsnow_permittivity_colbeck80_caseIII(freqGHz, density, lwc);
elseif density<550 && lwc >0.07
    eps_eff = wetsnow_permittivity_colbeck80_caseII(freqGHz, density, lwc);
end

end