function [JV, isPol] = calcJonesVector(antenna, fc, ang, isPol)
%CALCJONESVECTOR Calculate Jones vector for an antenna element or array
%   [JV, ISPOL] = COMM.INTERNAL.CALCJONESVECTOR(ANTENNA, FC, ANGLE, ISPOL)
%   calculates the Jones vector for an antenna element or array, ANTENNA,
%   given a frequency FC in Hz and direction ANGLE in degrees. ANTENNA is
%   from Phased System Toolbox or Antenna Toolbox. FC is real, positive
%   scalar. ANGLE is a 2-by-M matrix with each column in [azimuth;
%   elevation] form. The azimuth angle must be between [-180 180] degrees
%   and the elevation angle must be between [-90 90] degrees. ISPOL is a
%   logical scalar that is an optional input if already known for Phased
%   System antennas.
%   
%   JV is 2-by-M matrix with ech column being the Jones vector for the
%   corresponding column in ANGLE. The first and second rows represent
%   horizontal and vertical polarizations respectively. 
%   
%   ISPOL is a logical scalar indicating whether or not the antenna
%   supports polarization. Each column of JV is [1/sqrt(2); 1/sqrt(2)] for
%   an unpolarized antenna.  Unpolarized antennas are actually
%   randomly-polarized; practical RF antennas are polarized.
%
%   Example 1:
%       comm.internal.calcJonesVector(phased.CrossedDipoleAntennaElement, ...
%           28e9, [100 60 -110; 30 15 -35])
%       comm.internal.calcJonesVector(phased.IsotropicAntennaElement, ...
%           28e9, [100 60 -110; 30 15 -35])
%   Example 2:
%       comm.internal.calcJonesVector( ...
%           phased.URA('Element', phased.CrossedDipoleAntennaElement), ...
%           9e9, [-20 30 -10.2; 3 10 -32])
%   Example 3: 
%       fc = 100e8;
%       element = design(dipoleCrossed, fc);
%       comm.internal.calcJonesVector(element, ...
%           fc, [-20 30 -10.2; 3 10 -32])

% Copyright 2019-2024 The MathWorks, Inc.

% Jones vector will be normalized by magnitude and by dominant
% polarization's phase. By normalizing for phase, the Jones vector will
% match the nominal form (e.g. [0;1] for a vertical dipole); otherwise, the
% math requires all transmitter/receiver pairs use the same phase basis for
% their polarization vector.

isAbstractAntennaElement = isa(antenna, 'phased.internal.AbstractAntennaElement');
if isAbstractAntennaElement || isa(antenna, 'phased.internal.AbstractArray') || ...
        isa(antenna, 'phased.internal.AbstractSubarray')
    % Elements or arrays in PST

    if nargin < 4
        isPol = isPolarizationCapable(antenna);
    end
    if isPol
        if isAbstractAntennaElement
            resp = step(antenna, fc, ang);
        else % array
            ar = phased.ArrayResponse( ...
                'SensorArray', antenna, ...
                'EnablePolarization', true);
            resp = step(ar, fc, ang);
        end
        JV = normalizeJonesVector([resp.H.'; resp.V.']);
    else
        JV = ones(2, size(ang, 2))/sqrt(2);
    end
elseif isa(antenna, 'em.Antenna') || isa(antenna, 'em.Array') || ...
       isa(antenna,'installedAntenna') || isa(antenna,'customAntennaStl')
    % Elements or arrays in Antenna Toolbox

    isPol = true;

    % Calculate electric field in far-field @ 100 wavelengths
    az = ang(1,:);
    el = ang(2,:);
    r = 100*(299792458/fc); % = 100*c/f = 100 wavelengths
    [x, y, z] = sph2cart(az*(pi/180), el*(pi/180), r); 
    e = EHfields(antenna, fc, [x; y; z], "CoordinateSystem","Spherical");
    JV = normalizeJonesVector(e(1:2,:));
else % Default 'isotropic' or arrayConfig object
    isPol = false;
    JV = ones(2, size(ang, 2))/sqrt(2);
end

end

function JV = normalizeJonesVector(JV)
for col = 1:size(JV,2)
    thisJV = JV(:,col);
    if isequal(thisJV, [0;0])
        % Assume polarized equally in both directions
        JV(:,col) = sqrt(2)/2;
    else
        % Normalize by dominant polarization's phase, then magnitude
        [~,I] = max(abs(thisJV));
        thisJV = thisJV./thisJV(I);
        JV(:,col) = thisJV./norm(thisJV);
    end
end
end