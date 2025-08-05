function [lats, lons, actualResolutions, d] = sampleGreatCircle(sourceSite, targetSites, resolution, varargin)
%sampleGreatCircle   Sample great circle path between sites
%   [LATS,LONS,ACTRES,D] = sampleGreatCircle(SITE1,SITE2,RES) returns
%   latitudes and longitudes of sample points along the great circle path
%   with fixed step size. The actual step size resolution ACTRES is
%   computed as the nearest step size below RES that can produce a fixed
%   step size. The values ACTRES, RES, and total distance of the great
%   circle path D are expressed in meters.
%
%   The inputs SITE1 and SITE2 must be antenna sites, where SITE1 must be
%   scalar and SITE2 can be an array. The outputs ACTRES and D are double
%   arrays of size M-by-1, where M is the number of sites in SITE2. If
%   SITE2 is scalar, then the outputs LATS and LONS are double column
%   vectors. Otherwise, LATS and LONS are cell arrays of size M-by-1, where
%   each cell contains the double column vector corresponding to the sample
%   locations between SITE1 and the corresponding element of SITE2.

%   Copyright 2018-2019 The MathWorks, Inc.

validateattributes(sourceSite,{'txsite','rxsite'},{'scalar'},'sampleGreatCircle','',1);
validateattributes(targetSites,{'txsite','rxsite'},{'nonempty'},'sampleGreatCircle','',2);
validateattributes(resolution,{'numeric'},{'positive','finite'});

% Process optional name/value pairs
p = inputParser;
p.addParameter('Map', []);
p.addParameter('SourceAntennaSiteCoordinates', []);
p.addParameter('TargetAntennaSiteCoordinates', []);
p.parse(varargin{:});

% Get source site antenna coordinates
map = rfprop.internal.Validators.validateMapTerrainSource(p, 'sampleGreatCircle');
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSite, map, 'sampleGreatCircle');
[srclat, srclon] = sourceCoords.geodetic;

% Get target site antenna coordinates
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'sampleGreatCircle');
[tgtlat, tgtlon] = targetCoords.geodetic;

% Compute path locations from source to all targets
[lats,lons,actualResolutions,d] = rfprop.internal.MapUtils.sampleGreatCircle(...
    srclat, srclon, tgtlat, tgtlon, resolution);
end
