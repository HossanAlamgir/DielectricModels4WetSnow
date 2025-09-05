
function Es = wetsnow_permittivity_tinga73(frequency, density, liquid_water, ice_permittivity_model, water_permittivity_model)

    %     The function, converted from smrt permitivity module with the same
    %     func name,
    %     """Computes the effective permittivity proposed by Tinga et al. 1973 for three-component mixing. The component 1 is the background ("a" here),
    %     the compoment 2 ("w" here) is a spherical shell surrounding the component 3 ("i" here).
    %
    %      It was used by Tiuri as well as T. Mote to compute wet snolw permittivity.
    %
    %     Tinga, W.R., Voss, W.A.G. and Blossey, D. F.: General approach to multiphase dielectric mixture theory.
    %     Journal of Applied Physics, Vol.44(1973) No.9,pp.3897-3902.
    %     doi: /10.1063/1.1662868
        
    %     Tiuri, M. and Schultz, H., Theoretical and experimental studies of microwave radiation from a natural snow field. In Rango, A. , ed.
    %     Microwave remote sensing of snowpack properties. Proceedings of a workshop ... Fort Collins, Colorado, May 20-22, 1980.
    %     Washington, DC, National Aeronautics and Space Center, 225-234. (Conference Publication 2153.)
    %
    % frequency GHz
    % temperature K,
    % density kg/m^3
    % liquid_water ,
    % ice_permittivity_model, water_permittivity_model
    % """

%Constatnts
DENSITY_OF_ICE = 916.7;
DENSITY_OF_WATER = 1000;
FREEZING_POINT = 273.15;

if ~exist('temperature','var')
    temperature = FREEZING_POINT;
end

if temperature < FREEZING_POINT && any(liquid_water > 0)
    error('Liquid water is positive and temperature is set below freezing. This seems incompatible.');
end

% wetness W is the weight percentage of liquid water contained in the snow
W = liquid_water * DENSITY_OF_WATER / (liquid_water * DENSITY_OF_WATER + (1 - liquid_water) * DENSITY_OF_ICE);

% equation for spheres. Here we rather defined V to avoid the exponentiation
% ri = 0.5e-3  # the result is independent on this value, because only ratio rw/ri or ra/ri or rw/ra are used

% rw = ri * (1 + DENSITY_OF_ICE / DENSITY_OF_WATER * W / (1 - W))**(1 / 3)

% ra = ri * ((DENSITY_OF_ICE / density) * (1 + W / (1 - W)))**(1 / 3)

Vw_i = 1 + DENSITY_OF_ICE / DENSITY_OF_WATER * W / (1 - W);
Va_i = (DENSITY_OF_ICE / density) * (1 + W / (1 - W));

% Set default models for water and ice permittivity if not provided
if nargin < 5 || isempty(water_permittivity_model)
    water_permittivity_model = @water_permittivity_maetzler87; % Use a default model
end
if nargin < 6 || isempty(ice_permittivity_model)
    ice_permittivity_model = @ice_permittivity_maetzler06; % Use a default model
end

eps_a = 1;
eps_w = water_permittivity_model(frequency, temperature);
eps_i = ice_permittivity_model(frequency, temperature); % this must be dry ice !

alpha = 2 * eps_w + eps_i;
diff_wi = eps_w - eps_i;
diff_wa = eps_w - eps_a;

denominator = (2 * eps_a + eps_w) * alpha - 2 * (1 / Vw_i) * diff_wa * diff_wi ...
    - (Vw_i / Va_i) * diff_wa * alpha ...
    + (1 / Va_i) * diff_wi * (2 * eps_w + eps_a);

Es = eps_a * (1 + 3 * ((Vw_i / Va_i) * diff_wa * alpha - (1 / Va_i) * diff_wi * (2 * eps_w + eps_a)) / denominator);

% # it is possible to compute the square_field_ratio_tinga73
% # using np.abs(eps_w * eps_a / denominator)**2
% # to be implemented

end