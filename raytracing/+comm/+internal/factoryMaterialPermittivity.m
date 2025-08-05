function [epsilon, sigma, complexEpsilon] = factoryMaterialPermittivity(mtl, fc)
%FACTORYMATERIALPERMITTIVITY Calculate permittivity, conductivity for factory materials
%   [EPSILON, SIGMA, COMPLEXEPSILON] =
%   COMM.INTERNAL.FACTORYMATERIALPERMITTIVITY(MLT, FC) calculates the real
%   relative permittivity, EPSILON, conductivity, SIGMA, and complex
%   relative permittivity, COMPLEXEPSILON, for factory material, MTL, at
%   the carrier frequency, FC, in Hz. MLT must be one of "plasterboard" |
%   "ceiling-board" | "chipboard" | "floorboard" | "concrete" | "brick" |
%   "wood" | "glass" | "metal" | water" | "vegetation" | "loam" |
%   "perfect-reflector".
%   
%   For "water", "loam" and "vegetation", 20 degree temperature in Celsius
%   is assumed. For "loam", sand percentage of 41.96%, clay percentage of
%   8.53%, silt percentage of 49.51%, gravity of 2.70, volumetric water
%   content of 0.5 and bulk density of 1.5781 are assumed. For
%   "vegetation", gravimetric water content of 0.5 is assumed. For
%   "perfect-reflector", the relative permittivity and conductivity are
%   NaNs.
%
%   The output EPSILON is real-valued relative permittivity. The output
%   SIGMA is real-valued conductivity in S/m. The output COMPLEXEPSILON is
%   complex-valued relative permittivity, calculated by
%
%       COMPLEXEPSILON = EPSILON - 1i*SIGMA/(2*pi*FC*EPSILON0)
%
%   where EPSILON0 is the electric constant (permittivity of free space) in
%   F/m which is equal to 8.854187817e-12. All the outputs are scalars. 
%
%   References:
%   [1] ITU-R Recommendation P.2040-3, "Effects of building materials and
%       structures on radiowave propagation above about 100 MHz".
%       International Telecommunications Union - Radiocommunications Sector
%       (ITU-R), August 2023.
%   [2] ITU-R P.527-5, Electrical Characteristics of the Surface of the
%       Earth. International Telecommunications Union - Radiocommunications
%       Sector (ITU-R), August 2019.
%   [3] ITU-R P.527-6, Electrical Characteristics of the Surface of the
%       Earth. International Telecommunications Union - Radiocommunications
%       Sector (ITU-R), September 2021.
%
%   See also buildingMaterialPermittivity, earthSurfacePermittivity.

% Copyright 2019-2023 The MathWorks, Inc.
%#codegen

switch lower(mtl)
    case {"concrete","brick","wood","glass","metal", ... % ITU-R P.2040-1.
          "plasterboard","ceiling-board","chipboard","floorboard"} 
        [epsilon, sigma, complexEpsilon] = ...
            buildingMaterialPermittivity(mtl, fc);
    case "water"            % ITU-R P.527-4
        T = 20;
        [epsilon, sigma, complexEpsilon] = ...
            earthSurfacePermittivity('pure-water', fc, T);
    case "loam"             % ITU-R P.527-4, Table 1
        T = 20; 
        P_SAND = 41.96;
        P_CLAY = 8.53; 
        G = 2.70; 
        WC = 0.5; 
        BD = 1.5781;
        [epsilon, sigma, complexEpsilon] = ...
            earthSurfacePermittivity('soil', fc, T, ...
            P_SAND, P_CLAY, G, WC, BD);
    case "vegetation"   % ITU-R P.527-4
        T = 20;
        WC = 0.5;
        [epsilon, sigma, complexEpsilon] = ...
            earthSurfacePermittivity('vegetation', fc, T, WC); 
    otherwise % 'perfect-reflector' or any unknown material
        epsilon = nan;
        sigma = nan;
        complexEpsilon = nan - 1i*nan;
end

end

% [EOF]