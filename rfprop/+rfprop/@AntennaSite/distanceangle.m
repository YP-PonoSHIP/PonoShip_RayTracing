function [ds, az, el] = distanceangle(sourceSites, targetSites, varargin)
%distanceangle   Distance and angle between sites
%   [D,AZ,EL] = distanceangle(SITE1,SITE2) returns the distance in meters
%   and azimuth and elevation angles in degrees between SITE1 and SITE2,
%   using Euclidean path.

%   Copyright 2017-2020 The MathWorks, Inc.

% Function allows for getting distance and angles between sites in a single
% call, which is useful for better performance than separate calls to
% distance and angle.

% Process optional name/value pairs
p = inputParser;
p.addParameter('Map', []);
p.addParameter('SourceAntennaSiteCoordinates', []);
p.addParameter('TargetAntennaSiteCoordinates', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.SourceAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(sourceSites, targetSites);
else
    usingCartesian = strcmp(p.Results.SourceAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

% Get antenna site coordinates
if usingCartesian
    map = rfprop.internal.Validators.validateCartesianMap(p);
    coordSys = 'cartesian';
else
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'distanceangle');
    coordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, coordSys);
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSites, map, 'distanceangle');
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'distanceangle');

% Get Cartesian position of target sites with respect to source sites. Each
% returned variable is coordinate matrix with one column for each source
% site.
[X,Y,Z] = targetCoords.position(sourceCoords);

% Get distance and angles from Cartesian coordinates
[az,el,ds] = cart2sph(X,Y,Z);

% Convert angles to degrees in proper range
az = wrapTo180(rad2deg(az));
el = rad2deg(el);