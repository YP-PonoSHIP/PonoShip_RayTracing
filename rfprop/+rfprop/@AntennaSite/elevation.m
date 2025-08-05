function Z = elevation(sites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate inputs and terrain setting
validateattributes(sites,{'txsite','rxsite'},{'nonempty'}, ...
    'elevation','',1);

% Validate sites are all geographic
rfprop.internal.Validators.validateGeographicSites(sites, 'elevation')

% Process optional name/value pairs
p = inputParser;
p.addParameter('Map', []);
p.parse(varargin{:});

map = rfprop.internal.Validators.validateMapTerrainSource(p, 'elevation');
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'geographic');
% Get terrain
coords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(sites, map);
Z = coords.SurfaceHeightAboveGeoid;

% Snap very small values to 0
tolerance = 1e-3;
Z(abs(Z)<tolerance) = 0;