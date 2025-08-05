function [G, az, el] = gainPattern(site,fq)
%gainPattern   Antenna gain pattern
%   [G, AZ, EL] = gainPattern(SITE, FQ) returns the gain pattern G of
%   SITE's antenna for the frequency FQ. The gain is calculated over
%   the angle vectors (in degrees) AZ and EL.

%   Copyright 2021 The MathWorks, Inc.

ant = site.Antenna;

% Return early for 'isotropic' antenna
if ischar(ant)
    G = 0;
    az = 0;
    el = 0;
    return
end

% Compute angles to compute gain over
res = rfprop.Constants.GainAntennaPatternResolution;
azSweep = -180:res:180;
if azSweep(end) ~= 180
    azSweep(end + 1) = 180;
end
elSweep = -90:res:90;
if elSweep(end) ~= 90
    elSweep(end + 1) = 90;
end

[G, az, el] = pattern(ant, fq, azSweep, elSweep);
