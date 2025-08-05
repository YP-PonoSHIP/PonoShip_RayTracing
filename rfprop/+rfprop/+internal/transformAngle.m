function [taz, tel] = transformAngle(locAng, az, el)
%transformAngle   Transform azimuth/elevation angles to local coordinate system angles
%   [TAZ, TEL] = transformAngle(LOCALANG, AZ, EL) transforms the
%   azimuth/elevation angles AZ and EL (degrees) to azimuth/elevation
%   angles TAZ and TEL (degrees) in the local coordinate system defined by
%   LOCALANG, which defines the angle of the x-axis of the local coordinate
%   system. LOCALANG is specified as a numeric scalar or 2-by-1 vector. If
%   specified as a scalar, LOCALANG is azimuth angle to the local x-axis.
%   If specified as a 2-by-1 vector, the first element is azimuth angle and
%   the second element is elevation angle.
%
%   % Example: Transform east to local reference
%
%   % Define local reference angle oriented northeast with downtilt
%   localAngle = [45;-10];
%
%   % Transform east direction
%   [taz,tel] = rfprop.internal.transformAngle(localAngle,0,0)

%   Copyright 2018-2020 The MathWorks, Inc.

% Calculate local axes based local angles
R = rfprop.internal.coordinateTransformationMatrix(locAng);

% Convert spherical angles in degrees from GCS to LCS
angLCS = global2localcoord( ...
    [az(:), el(:), ones(length(az), 1)]', 'ss', zeros(3,1), R);

% Get azimuth and elevation angles in LCS and convert them to columns
taz = angLCS(1, :)'; 
tel = angLCS(2, :)';
