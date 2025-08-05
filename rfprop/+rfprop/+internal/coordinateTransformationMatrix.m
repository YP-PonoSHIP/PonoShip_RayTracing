function R = coordinateTransformationMatrix(localAng)
%coordinateTransformationMatrix Calculate rotation matrix from LCS to GCS
%   R = coordinateTransformationMatrix(LOCALANG) calculates the 3-by-3
%   unitary rotation matrix (local axes), R, from LCS to GCS, given the LCS
%   angles, LOCALANG, in the form of [alpha] or [alpha beta]. Alpha is the
%   azimuth angle, measured counterclockwise from x-axis to y-axis, around
%   z-axis. Beta is the elevation angle, measured counterclockwise from
%   x-axis to z-axis, around y-axis.
% 
%   See also sph2cart, cart2sph

%   Copyright 2020 The MathWorks, Inc.

% Azimuth (bearing) angle in degrees
alpha = localAng(1); 
% Elevation (downtilt) angle in degrees
if numel(localAng) > 1 
    beta = localAng(2);
else
    beta = 0;
end

% Azimuth rotation around a-axis: x --> y, y --> x, z stays
azRot = [cosd(alpha) sind(alpha) 0; -sind(alpha) cosd(alpha) 0; 0 0 1];

% Elevation rotation around y-axis: x --> z, z --> x, y stays
elRot = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];

% Local axes from the two rotations 
R = (elRot * azRot)';
 
end