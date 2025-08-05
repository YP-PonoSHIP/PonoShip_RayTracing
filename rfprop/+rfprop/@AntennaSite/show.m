function show(sites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Check if too many sites to show
if numel(sites) > rfprop.Constants.MaxNumSitesShow
    error(message('shared_channel:rfprop:ShowTooManySites', rfprop.Constants.MaxNumSitesShow));
end

% Process optional name/value pairs
p = inputParser;
p.addParameter('Icon', []);
p.addParameter('Animation', '');
p.addParameter('Map', []);
p.addParameter('EnableWindowLaunch', true);
p.addParameter('ClusterMarkers', false);
p.addParameter('ShowAntennaHeight', true);
p.addParameter('IconSize', [36,36]);
p.addParameter('IconAlignment', 'bottom');
p.addParameter('AntennaSiteCoordinates', []);
p.parse(varargin{:});

% Validate that all sites use the same coordinate system
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(sites);

% Get Site Viewer visibility state and validate web graphics
viewer = rfprop.internal.Validators.validateMap(p, 'show', usingCartesian);
if usingCartesian
    siteCoordSys = 'cartesian';
else
    siteCoordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(viewer, siteCoordSys);
isViewerInitiallyVisible = viewer.Visible;

% Validate and get parameters
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'show');
[isCustomIcon, icon] = validateIcon(p);
[isCustomCluster, clusterMarkers] = validateClusterMarkers(p);
[isDefaultShowAntennaHeight, showStem] = validateShowAntennaHeight(p);
[isCustomSize, iconSize] = validateIconSize(p);
[isCustomAlignment, iconAlignment] = validateIconAlignment(p);
coords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.AntennaSiteCoordinates, sites, viewer, 'show');

% Pack show() parameters into same-sized cell arrays
numSites = numel(sites);
IDs = cell(1, numSites);
iconSizes = cell(1, numSites);
iconAlignments = cell(1, numSites);
icons = cell(1, numSites);
clusterMarkerSites = cell(1, numSites);
descriptions = cell(1, numSites);
showStems = true(1, numSites);
stemBases = zeros(1, numSites);
names = IDs;
graphicsToRemove = {};
dec = rfprop.Constants.MaxLatLonInfoDecimalPlaces;
% Use the same numericFormat as mat2str
numericFormat = "%.15g";
for k = 1:numSites
    siteID = sites(k).UID;
    IDs{k} = siteID;

    % If the site already exists in the siteviewer's graphics, then we may
    % need to update stale graphics.
    viewer.initializeSiteGraphics(siteID);
    siteGraphicLocation = viewer.SiteGraphics.(siteID).marker.Location;
    if usingCartesian
        currentLocation = coords.AntennaPosition(k, :);
    else
        currentLocation = [coords.LatitudeLongitude(k,:) coords.AntennaHeight(k)];
    end
    if (~all(isnan(siteGraphicLocation)) && ...
            ~all(siteGraphicLocation == currentLocation))
        % If the site was moved, then we must clear graphics associated
        % with the site.
        graphicsToRemove = [graphicsToRemove; viewer.getGraphicsIDs(siteID)];
        viewer.clearSiteGraphics(siteID);
        viewer.SiteGraphics.(siteID) = rfprop.AntennaSite.SiteGraphicsTemplate;
    end
    viewer.SiteGraphics.(siteID).marker.ID = IDs{k};
    viewer.SiteGraphics.(siteID).marker.Location = currentLocation;
    names{k} = sites(k).Name;
    % Set the hidden visualization properties of the site
    if (isCustomAlignment)
        sites(k).IconAlignment = iconAlignment;
    end
    if (isCustomSize)
        sites(k).IconSize = iconSize;
    end
    if (isCustomIcon)
        sites(k).Icon = icon;
    end
    if (isCustomCluster)
        sites(k).ClusterMarkers = clusterMarkers;
    end
    
    % Modify the antenna's ShowAntennaHeight value if the user explicitly
    % specified it (it is not default)
    stemVisibility = sites(k).ShowAntennaHeight;
    if ~isDefaultShowAntennaHeight
        sites(k).ShowAntennaHeight = showStem;
        stemVisibility = showStem;
    elseif isDefaultShowAntennaHeight && strcmp(sites(k).ShowAntennaHeight, 'auto')
        % Otherwise, if the default antenna height is being used AND it was
        % never changed (which means it still uses 'auto'), then update the
        % value based on various conditions.
        stemVisibility = viewer.UseTerrain || viewer.HasBuildings || viewer.IsCartesian;
        % If showAntennaHeight was set by the user previously, ShowAntennaHeight will not be 'auto'
        % Meaning we must honor the user-specified value even if
        % isDefaultShowAntennaHeight is true.
    end
    iconSizes{k} = sites(k).IconSize;
    iconAlignments{k} = sites(k).IconAlignment;
    icons{k} = sites(k).Icon;
    clusterMarkerSites{k} = sites(k).ClusterMarkers;
    showStems(k) = stemVisibility;
    % Description
    % Use 'compose' to avoid performance cost of 'mat2str'
    if usingCartesian
        descriptions{k} = ['Position: [' char(compose(numericFormat, round(sites(k).AntennaPosition(1), dec))) ', ' ...
            char(compose(numericFormat, round(sites(k).AntennaPosition(2), dec))) ', ' ...
            char(compose(numericFormat, round(sites(k).AntennaPosition(3), dec))) ']'];
        if showStems(k)
            stemBases(k) = computeStemBase(viewer, sites(k));
        end
    else
        descriptions{k} = [char(compose(numericFormat,round(coords.LatitudeLongitude(k,1), dec))) ',' char(compose(numericFormat,round(coords.LatitudeLongitude(k,2), dec))) ...
            '<br />' 'Antenna height: ' char(compose(numericFormat,round(sites(k).AntennaHeight, dec))) ' m', ...
            '<br />' 'Surface elevation: ' char(compose(numericFormat,round(coords.SurfaceHeightAboveGeoid(k)))) ' m'];
    end
end

viewer.remove(graphicsToRemove);

iconUrls = cell(1, numSites);
defaultIcon = sites(1).DefaultMarkerIcon;
for h = 1:numSites
    if  (isempty(icons{h}))
        % Empty means the default icon is being used.
        iconUrls{h} = defaultIcon;
    elseif (isempty(which(icons{h})))
        % Icon not found on path
        if (exist(icons{h}, 'file') == 2)
            % Full path given
            iconUrls{h} = globe.internal.ConnectorServiceProvider.getResourceURL(icons{h}, "rfprop" + sites(k).UID);
        else
            % No path given, cannot find icon. Switch to default icon.
            % This case only occurs if a site is saved with an icon used,
            % and a different enviroment is entered where the icon is no
            % longer on the path.
            warning(message('shared_channel:rfprop:IconFileNotFound', icons{h}));
            iconUrls{h} = defaultIcon;
        end
    else
        % Custom icon in current folder
        iconUrls{h} = globe.internal.ConnectorServiceProvider.getResourceURL(which(icons{h}), "rfprop" + sites(k).UID);
    end
end
icons = iconUrls;

if isViewerInitiallyVisible && ~viewer.Visible
    return % Abort if Site Viewer has been closed (test before show)
end

% Package inputs
if usingCartesian
    position = coords.AntennaPosition;
    if isempty(animation)
        animation = 'fly';
    end
    if ~viewer.Visible
        animation = 'zoom';
    end
    args = {...
        'IconSize', iconSizes, ...
        'IconAlignment', iconAlignments, ...
        'ID', IDs, ...
        'ValidateIcon', false, ...
        'Name', names, ...
        'Description', descriptions, ...
        'ShowStem', showStems, ...
        'Animation', animation, ...
        'StemBase', stemBases};
else
    position = [coords.LatitudeLongitude, coords.AntennaHeight];
    args = {...
        'ID', IDs, ...
        'Name', names, ...
        'GroundElevation', coords.SurfaceHeightAboveTerrainReference, ...
        'IconSize', iconSizes, ...
        'IconAlignment', iconAlignments ,...
        'ShowStem', showStems, ...
        'Animation', animation, ...
        'Description', descriptions, ...
        'EnableWindowLaunch', enableWindowLaunch, ...
        'ClusterMarkers', clusterMarkerSites, ...
        'ValidateIcon', false};
end

viewer.Visualizer.marker(position, icons, args{:});

if (enableWindowLaunch)
    if ~viewer.Visible
        viewer.Visible = true;
    end
    viewer.bringToFront;
end
end

function [isCustomIcon, icon] = validateIcon(p)

try
    isCustomIcon = ~ismember('Icon', p.UsingDefaults);
    icon = [];
    % Validate user-defined icon
    if isCustomIcon
        icon = p.Results.Icon;

        % Validate type/size
        validateattributes(icon, {'char','string'}, {'scalartext'}, ...
            'show', 'Icon');
        if isstring(icon)
            icon = char(icon);
        end
        if strcmpi(icon, 'default')
            icon = [];
        else
            % Validate file extension
            [~, ~, ext] = fileparts(icon);
            if ~ismember(ext, {'.png', '.svg'})
                error(message('shared_channel:rfprop:IconFormatNotSupported'));
            end

            % Validate file existence
            iconPath = which(icon);
            if isempty(iconPath) % Not found on MATLAB path
                if exist(icon, 'file') ~= 2 % Full or relative pathname
                    error(message('shared_channel:rfprop:IconDoesNotExist'));
                end
            end
        end
    end

catch e
    throwAsCaller(e)
end
end

function [isCustomCluster, clusterMarkers] = validateClusterMarkers(p)
try
    clusterMarkers = p.Results.ClusterMarkers;
    validateattributes(clusterMarkers, {'logical'}, {'nonsparse','scalar'}, ...
        'show', 'ClusterMarkers');

    isCustomCluster = ~ismember('ClusterMarkers', p.UsingDefaults);
catch e
    throwAsCaller(e)
end
end

function [isCustomSize, iconSize] = validateIconSize(p)
try
    isCustomSize = ~ismember('IconSize', p.UsingDefaults);
    iconSize = p.Results.IconSize;
    validateattributes(iconSize, {'numeric'}, {'size', [1,2], 'positive', 'finite', 'real', 'nonsparse', '<=', 1000}, 'show', 'IconSize');
catch e
    throwAsCaller(e);
end
end


function [isCustomAlignment, iconAlignment] = validateIconAlignment(p)
try
    isCustomAlignment = ~ismember('IconAlignment', p.UsingDefaults);
    iconAlignment = p.Results.IconAlignment;
    validAlignments = ["bottom", "center", "top"];
    iconAlignment = validatestring(iconAlignment, validAlignments, 'show', 'IconAlignment');
catch e
    throwAsCaller(e);
end
end

function [isDefaultShowAntennaHeight, showStem] = validateShowAntennaHeight(p)
try
    isDefaultShowAntennaHeight = ismember('ShowAntennaHeight', p.UsingDefaults);
    showStem = p.Results.ShowAntennaHeight;
    if ~isDefaultShowAntennaHeight
        validateattributes(showStem, {'logical'}, {'nonsparse','scalar'}, ...
            'show', 'ShowAntennaHeight');
    end
catch e
    throwAsCaller(e)
end
end

function stemBase = computeStemBase(viewer, site)
% Compute the Z height of the base of the stem.
map = viewer.ModelMesh;
if ~strcmp(map, 'none') && ~isempty(map)
    warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    warnCleanup = onCleanup(@()warning(warnState));
    tri = Geometry.meshToTriangulation(map);
    % Shoot a ray downwards to see if we hit any part of the model. If we
    % do, the stem will be drawn down to that surface.
    env = struct( ...
        "Triangulation", tri, ...
        "RayTracer", matlabshared.internal.StaticSceneRayTracer(tri), ...
        "SharpEdgeFlags", comm.internal.geometry.getSharpEdges(tri));
    dir = [0 0 -1];
    [firstpt, ~,isobstructed] = comm.internal.geometry.firstIntersect( ...
        env, site.AntennaPosition', dir, 'ray');
else % If there is no model, then the stem will be drawn down to (0,0)
    isobstructed = false;
end

if (isobstructed)
    stemBase = firstpt(3);
else
    stemBase = 0;
end
end