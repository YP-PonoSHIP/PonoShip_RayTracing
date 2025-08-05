function [d, heading] = greatCircleDistance(sourceSites, targetSites, varargin)
%greatCircleDistance   Great circle distance and heading between sites
%   [D,HEADING] = greatCircleDistance(SITE1,SITE2) returns the great circle
%   distance and heading from SITE1 to SITE2. The distance D is expressed
%   as meters, and the heading HEADING is expressed in degrees clockwise
%   from north.
%
%   The inputs SITE1 and SITE2 can be scalars or arrays. The outputs are
%   arrays of size M-by-N, where M is the number of sites in SITE2 and N is
%   the number of sites in SITE1.

%   Copyright 2017-2019 The MathWorks, Inc.

% Process optional name/value pairs
p = inputParser;
p.addParameter('SourceAntennaSiteCoordinates', []);
p.addParameter('TargetAntennaSiteCoordinates', []);
p.parse(varargin{:});

% Get source antenna site coordinates
defaultTerrain = siteviewer.defaultTerrainName;
if (strcmp(defaultTerrain, 'none'))
    map = defaultTerrain;
else
    map = terrain.internal.TerrainSource.createFromSettings(siteviewer.defaultTerrainName);
end
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSites, map, 'greatCircleDistance');
[srclat, srclon] = sourceCoords.geodetic;

% Get target antenna site coordinates
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'greatCircleDistance');
[tgtlat, tgtlon] = targetCoords.geodetic;

d = zeros(numel(targetSites),numel(sourceSites));
heading = d;
for srcInd = 1:numel(srclat)
    [d(:,srcInd), heading(:,srcInd)] = rfprop.internal.MapUtils.greatCircleDistance(...
        srclat(srcInd), srclon(srcInd), tgtlat, tgtlon);
end