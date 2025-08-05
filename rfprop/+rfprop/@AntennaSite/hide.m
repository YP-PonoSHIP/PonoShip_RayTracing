function hide(sites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Process optional name/value pairs
p = inputParser;
p.addParameter('Map', []);
p.parse(varargin{:});

% Validate that all sites use the same coordinate system
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(sites);

% Get viewer and hide
viewer = rfprop.internal.Validators.validateMap(p, 'hide', usingCartesian);
if usingCartesian
    siteCoordSys = 'cartesian';
else
    siteCoordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(viewer, siteCoordSys);
% Pack hide() parameters into same-sized cell arrays
numSites = numel(sites);
% Can't preallocate because there's no guarantee each site in "sites" has
% been plotted
IDs = {};
for k = 1:numSites
    if isfield(viewer.SiteGraphics, sites(k).UID)
        IDs = [IDs; viewer.SiteGraphics.(sites(k).UID).marker.ID]; %#ok<AGROW>
    end
end
viewer.remove(IDs);
end