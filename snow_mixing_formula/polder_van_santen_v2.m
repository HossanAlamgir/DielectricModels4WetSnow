function eps_eff  = polder_van_santen_v2(frac_volume, e0, eps, inclusion_shape)
% """ Calculates effective permittivity of snow by solution of quadratic Polder Van Santen equation for spherical inclusion.
%
% :param frac_volume: Fractional volume of inclusions
% :param e0: Permittivity of background (default is 1)
% :param eps: Permittivity of scattering material (default is 3.185 to compare with MEMLS)
% :param depol_xyz: [Optional] Depolarization factors, spherical isotropy is default. It is not taken into account here.
% :param length_ratio: Length_ratio. Used to estimate depolarization factors when they are not given.
% :param inclusion_shape: Assumption for shape(s) of brine inclusions. Can be a string for single shape, or a list/tuple/dict of strings for mixture of shapes. So far, we have the following shapes: "spheres" and "random_needles" (i.e. randomly-oriented elongated ellipsoidal inclusions).
%     If the argument is a dict, the keys are the shapes and the values are the mixing ratio. If it is a list, the mixing_ratio argument is required.
%     :param mixing_ratio: The mixing ratio of the shapes. This is only relevant when inclusion_shape is a list/tuple. Mixing ratio must be a sequence with length len(inclusion_shape)-1. The mixing ratio of the last shapes is deduced as the sum of the ratios must equal to 1.
%     :returns: Effective permittivity
%
%     **Usage example:**
%
%     ::
%
%     from smrt.permittivity.generic_mixing_formula import polder_van_santen
%     effective_permittivity = polder_van_santen(frac_volume, e0, eps)
%
%     # for a mixture of 30% spheres and 70% needles
%         effective_permittivity = polder_van_santen(frac_volume, e0, eps, inclusion_shape={"spheres": 0.3, "random_needles": 0.7})
%         # or
%         effective_permittivity = polder_van_santen(frac_volume, e0, eps, inclusion_shape=("spheres", "random_needles"), mixing_ratio=0.3)
%
%         .. todo::
%
%         Extend Polder Van Santen model to account for ellipsoidal inclusions
%
%         """

% Set default values for optional arguments
if nargin < 2 || isempty(e0)
    e0 = 1;
end
if nargin < 3 || isempty(eps)
    eps = 3.185;
end
if nargin < 4 || isempty(inclusion_shape)
    inclusion_shape = [];
end


% Check if frac_volume is valid (should be <= 1)
assert(all(frac_volume <= 1), 'The fractional volume is larger than 1: %g', max(frac_volume));


% Optional: You could implement depolarization factors if necessary
% depol_xyz = depolarization_factors(length_ratio);



% # Polder Van Santen / de Loor / Böttcher / Bruggeman formula
% # Solution of quadratic equation arising from eqn 9.2. in Sihvola: Electromagnetic Mixing Formulas and Applications
if (isempty(inclusion_shape)) || (inclusion_shape == "spheres")
    a_quad = 2;
    b_quad = eps - 2 * e0 - 3. * frac_volume * (eps - e0);
    c_quad = - eps * e0;

    % # Polder and Van Santen model, modified by de Loor (according to Shokr (1998) simplified by Hoekstra and Capillino (1971))
    % # Solution of quadratic equation arising from eqn (18) in Shokr (1998): 'Field Observations and Model Calculations of Dielectric Properties of Arctic Sea Ice in the Microwave C-Band', IEEE.
elseif inclusion_shape == "random_needles"
    a_quad = 1;
    b_quad = eps - e0 - 5. / 3. * frac_volume * (eps - e0);
    c_quad = - eps * (e0 + 1. / 3. * frac_volume * (eps - e0));

else
    error("inclusion_shape must be one of (or a list of) the following: 'spheres' (default) or 'random_needles'.")
end


eps_eff = (-1*b_quad + sqrt(b_quad.^2 - 4. * a_quad * c_quad)) / (2. * a_quad);




end
