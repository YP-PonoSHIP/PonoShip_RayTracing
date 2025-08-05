function varargout = coverage(txs, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate number of output arguments
nargoutchk(0,1)

% Validate site
validateattributes(txs,{'txsite'},{'nonempty'},'coverage','',1);

% Validate sites are all geographic
rfprop.internal.Validators.validateGeographicSites(txs, 'coverage')

% Process optional rxsite input
rx = [];
numTx = numel(txs);
args = varargin;
defaultRxGain = 2.1;
defaultRxHeight = 1;
if numel(args) > 0
    secondInput = args{1};
    if isa(secondInput,'rxsite')
        rx = secondInput;
        validateattributes(rx,{'rxsite'},{'scalar'},'coverage','',2);
        args = args(2:end); % Remove from input list
        
        % Set default receiver parameters from input. Generate default rx
        % gain for each tx, since they may have different frequencies.
        defaultRxGain = zeros(1,numTx);
        for k = 1:numTx
            defaultRxGain(k) = gain(rx,txs(k).TransmitterFrequency) - rx.SystemLoss;
        end
        defaultRxHeight = rx.AntennaHeight;
    end
end

% Add parameters
p = inputParser;
if mod(numel(args),2) % Odd number of inputs
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Type', 'power');
p.addParameter('SignalStrengths', []);
p.addParameter('Resolution', 'auto');
p.addParameter('ReceiverGain', defaultRxGain);
p.addParameter('ReceiverAntennaHeight', defaultRxHeight);
p.addParameter('Animation', '');
p.addParameter('MaxRange', 'auto');
p.addParameter('Colormap', 'jet');
p.addParameter('Colors', [])
p.addParameter('ColorLimits', []);
p.addParameter('Transparency', 0.4);
p.addParameter('ShowLegend', true);
p.addParameter('ReceiverLocationsLayout', []);
p.addParameter('MaxImageResolutionFactor', 5);
p.addParameter('RadialResolutionFactor', 2);
p.addParameter('Map', []);
p.parse(args{:});

% Get Site Viewer and validate web graphics
outputRequested = nargout > 0;
if ~outputRequested
    viewer = rfprop.internal.Validators.validateMap(p, 'coverage');
    isViewerInitiallyVisible = viewer.Visible;
else
    isViewerInitiallyVisible = false;
end

% Get coverage map type and corresponding values 
allowedTypes = {'efield','power'};
[type, defaultStrengths, defaultColorLimits, legendTitle] = ...
    rfprop.internal.Validators.validateType(p, allowedTypes, 'coverage');

% Validate parameters
[map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'coverage');
terrainSource = rfprop.internal.Validators.validateTerrainSource(map, 'coverage');
animation = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'coverage');
strengths = validateSignalStrengths(p, defaultStrengths);
colors = rfprop.internal.Validators.validateColors(p, 'coverage');
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'coverage');
cmap = rfprop.internal.Validators.validateColorMap(p, 'coverage');
transparency = rfprop.internal.Validators.validateTransparency(p, 'coverage');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'coverage');
pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'coverage');
rfprop.internal.Validators.validateMaxNumReflections(pm, 'coverage');
rfprop.internal.Validators.validateMaxNumDiffractions(pm, 'coverage'); % Disable 2 order diffractions
rxGain = rfprop.internal.Validators.validateReceiverGain(p, 'coverage');
rxAntennaHeight = rfprop.internal.Validators.validateReceiverAntennaHeight(p, 'coverage');
maxImageResFactor = p.Results.MaxImageResolutionFactor;
radialResFactor = p.Results.RadialResolutionFactor;

% Get site coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, map);
txslatlon = txsCoords.LatitudeLongitude;

% Validate dependent parameters
if isscalar(rxGain) % Expand scalar gain to match length of tx
    rxGain = repmat(rxGain,1,numTx);
end
maxrange = validateMaxRange(p, type, pm, txs, txslatlon, strengths, rxGain, map);
[res, isAutoRes] = rfprop.internal.Validators.validateResolution(p, maxrange, 'coverage');
datarange = rfprop.internal.Validators.validateDataRange(txslatlon, maxrange, res, ~strcmp(terrainSource,'none'));
rxLocationsLayout = rfprop.internal.Validators.validateReceiverLocationsLayout(p, pm, txslatlon, 'coverage');

if outputRequested
    % Do not show progress dialog
    generatingMapMsg = '';
    computingDataMsg = '';
else
    % Get color info
    colorsSpecified = ~ismember('Colors', p.UsingDefaults);
    colormapSpecified = ~all(ismember({'Colormap','ColorLimits'}, p.UsingDefaults));
    useColors = colorsSpecified || (numel(strengths) < 2 && ~colormapSpecified);
    if useColors
        if isempty(colors) && ismember(type,{'power','efield'})
            colors = [0 1 0]; % Default coverage map color (green)
        end
        colorData = struct('Levels',strengths,'Colors',colors);
    else
        colorData = struct('Colormap',cmap,'ColorLimits',clim);
    end
    
    % Clear old contour color info for each site. This is because a site can
    % never conflict with itself (since a site will replot it's own contour if
    % coverage is called on a site multiple times)
    % Also mark the contour graphics to be removed
    graphicsToRemove = {};
    for k = 1:numTx
        oldContourID = viewer.getGraphic(txs(k).UID, 'contour');
        if (~isempty(oldContourID))
            viewer.removeColorData(oldContourID);
            graphicsToRemove = [graphicsToRemove; oldContourID]; %#ok<AGROW>
        end
    end
    if ~isempty(viewer.LegendID)
        graphicsToRemove = [graphicsToRemove; viewer.LegendID];
    end
    
    % Generate a new ID so old IDs don't conflict with new images
    idNum = viewer.getId(1);
    contourID = ['contour' num2str(idNum{1})];
    
    % Validate image for color conflicts. This must be done before calculations
    % occur to avoid wasting the user's time when color conflict occurs.
    oldColorGraphics = viewer.ColorGraphics;
    viewer.checkForGraphicsConflict(type, contourID, colorData);
    resetColorGraphics = rfprop.internal.onExit(@()setColorGraphics(viewer,oldColorGraphics));
    
    % Show txsites. If Site Viewer is already open, keep current camera view.
    if isViewerInitiallyVisible && ~viewer.Visible
        return % Abort if Site Viewer has been closed
    end
    if isViewerInitiallyVisible
        showAnimation = 'none';
    else
        showAnimation = 'zoom';
    end
    show(txs,'Map',viewer,'Animation',showAnimation, ...
        'AntennaSiteCoordinates', txsCoords);
    
    % Show progress dialog. Use cleanup object to close dialog if a forced
    % exit occurs (error or Ctrl-C). Use three-stage dialog if using ray
    % tracing and single-stage dialog if not a terrain model. If using a
    % terrain model, dialog is launched in radialReceiverLocationsLayoutData.
    resetProgressDialog = rfprop.internal.onExit(@()hideProgressDialog(viewer));
    showWaitbar = pm.isMultipathModel;
    generatingMapMsg = message('shared_channel:rfprop:ProgressDialogGeneratingCoverageMap').getString;
    computingDataMsg = message('shared_channel:rfprop:ProgressDialogComputingCoverage').getString;
    if ~pm.requiresTerrain || strcmp(rxLocationsLayout,'grid')
        if showWaitbar
            msg = message('shared_channel:rfprop:ProgressDialogPreparingMapData').getString;
        else
            msg = generatingMapMsg;
        end
        viewer.showProgressDialog('Message', msg, ...
            'Indeterminate', true, ...
            'Cancelable', false);
    end
end

% Generate location grid containing data range from each transmitter site
if isa(map,'siteviewer')
    maxImageSize = map.MaxImageSize;
else
    maxImageSize = rfprop.Constants.DefaultMaxImageSize;
end
[latNorth, latSouth, lonEast, lonWest, animation] = ...
    rfprop.internal.MapUtils.geobounds(txslatlon, datarange, animation);
[gridlats, gridlons, res] = rfprop.internal.MapUtils.geogrid(...
    latNorth, latSouth, lonEast, lonWest, res, isAutoRes, maxrange, maxImageSize, 'coverage');
gridSize = size(gridlats);

% Compute and validate coverage map image size
imageSize = rfprop.internal.Validators.validateImageSize(...
    gridSize, maxImageResFactor, maxImageSize, res, 'coverage');

% Trim grid locations to those which are within data range
[latitude, longitude] = rfprop.internal.MapUtils.georange(...
    txs, gridlats(:), gridlons(:), datarange, terrainSource);

% Trim grid locations to those which are not inside buildings
if pm.isMultipathModel && isa(map,'siteviewer')
    isWithinBldg = isWithinAnyBuilding(map.BuildingsArray, latitude, longitude);
    latitude = latitude(~isWithinBldg);
    longitude = longitude(~isWithinBldg);
end

if ~outputRequested && siteViewerWasClosed(viewer)
    return % Abort if Site Viewer has been closed
end
if isempty(latitude)
    warning(message('shared_channel:rfprop:NoCoverageMapArea'))
    if outputRequested
        varargout = {[]};
    else
        viewer.removeColorData(viewer.SiteGraphics.(txs(k).UID).contour);
        viewer.remove(graphicsToRemove);
    end
    return
end

% Use try/catch in case cancel error thrown from waitbar dialog
try
    if strcmp(rxLocationsLayout,'grid')
        % Define rxsites at grid locations within data range
        rxs = rxsite(...
            'Name', 'internal.coveragesite', ... % Specify to avoid default site naming
            'Latitude', latitude, ...
            'Longitude', longitude, ...
            'AntennaHeight', rxAntennaHeight);
        
        % Launch waitbar phase of dialog
        if ~outputRequested && showWaitbar
            viewer.showProgressDialog('Message', computingDataMsg, ...
                'Value', 0, ...
                'ShowPercentage', true, ...
                'Indeterminate', false, ...
                'Cancelable', true);
        end
        
        % Compute signal strength at each rxsite in the grid.
        data = sigstrength(rxs, txs, pm, ...
            'Type', type, ...
            'ReceiverGain', rxGain, ...
            'Map', mapStruct, ...
            'TransmitterAntennaSiteCoordinates', txsCoords);
        
    else
        % Compute data and corresponding locations using radial layout
        [latitude, longitude, ~, data] = radialReceiverLocationsLayoutData(txs, txsCoords, type, ...
            generatingMapMsg, computingDataMsg, pm, mapStruct, res, radialResFactor, datarange, ...
            lonWest, lonEast, latSouth, latNorth, rxGain, rxAntennaHeight);
    end
catch e
    if strcmp(e.identifier, "shared_channel:rfprop:ProgressDialogCancelled")
        return
    else
        rethrow(e)
    end
end

% Launch final phase of status dialog 
if ~outputRequested
    viewer.showProgressDialog('Message', generatingMapMsg, ...
        'Indeterminate', true, ...
        'Cancelable', false);
end

% Merge data of rx sites for use in single image, where signal strength at
% an rx site is the max strength due to any tx.
if numTx > 1
    data = max(data, [], 1);
end

% Remove NaN-data
isnn = ~isnan(data);
data = data(isnn);
latitude = latitude(isnn);
longitude = longitude(isnn);

if ~outputRequested 
    % Show rxsite
    if siteViewerWasClosed(viewer)
        return
    end
    if isa(rx,'rxsite')
        show(rx,'Map',viewer,'Animation','none','EnableWindowLaunch',false);
    end
end

% Create propagation data container
if strcmpi(type,'power')
    dataVariableName = 'Power';
else
    dataVariableName = 'Efield'; 
end
pd = propagationData(latitude, longitude, ...
    'Name', message('shared_channel:rfprop:CoveragePropagationDataTitle').getString, ...
    dataVariableName, data(:));

if outputRequested
    % Return propagationData with contour properties set except for ID,
    % which is not set since there is no plot to associate with
    pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize);
    varargout = {pd};
    return
end

% It is now safe to remove old contours associated with sites
viewer.remove(graphicsToRemove);

% Return early if no coverage map data meets minimum signal strength,
% including if no rays were found (all data = -Inf)
minstrength = min(strengths,[],"all");
[mindata, maxdata] = bounds(data,"all");
if maxdata < minstrength
    if maxdata == -Inf % No rays found
        warning(message("shared_channel:rfprop:NoCoverageMapAreaHasFiniteSignalStrength"))
    else
        warning(message("shared_channel:rfprop:NoCoverageMapAreaMeetsMinSignalStrength", ...
            num2str(minstrength),num2str(mindata),num2str(maxdata)))
    end
    
    return
end

% Create contour map
pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize,contourID);
contourArgs = {'Type', type, ...
    'Map', viewer, ...
    'Levels', strengths, ...
    'Transparency', transparency, ...
    'ShowLegend', showLegend, ...
    'ImageSize', imageSize, ...
    'LegendTitle', legendTitle, ...
    'ValidateColorConflicts', false};
if useColors
    contourArgs = [contourArgs, ...
        'Colors', colors];
else
    contourArgs = [contourArgs, ...
        'ColorLimits', clim, ...
        'Colormap', cmap];
end
contour(pd, contourArgs{:});

% Set the IDs of the appropriate Site Graphics
for i = 1:numel(txs)
    viewer.SiteGraphics.(txs(i).UID).contour = contourID;
    if (showLegend)
        viewer.SiteGraphics.(txs(i).UID).legend = ['legend' contourID];
    end
end

% Hide progress dialog
viewer.hideProgressDialog;

% Cancel forced exit cleanup, since normal execution has completed
cancel(resetProgressDialog);
if ~outputRequested
    cancel(resetColorGraphics);
end
end

function wasClosed = siteViewerWasClosed(viewer)

wasClosed = viewer.LaunchWebWindow && ~viewer.Visible;
end

function strengths = validateSignalStrengths(p, defaultStrengths)

try
    if ismember('SignalStrengths', p.UsingDefaults)
        strengths = defaultStrengths;
    else
        strengths = p.Results.SignalStrengths;
        validateattributes(strengths, {'numeric'}, ...
            {'real','finite','nonnan','nonsparse','vector','nonempty'}, 'coverage', 'SignalStrengths');
        if ~iscolumn(strengths)
            strengths = strengths(:);
        end
        strengths = double(strengths);
    end
catch e
    throwAsCaller(e);
end
end

function maxrange = validateMaxRange(p, type, pm, txs, txslatlon, strengths, rxGain, map)

try
    maxrange = p.Results.MaxRange;
    numSites = numel(txs);
    if ischar(maxrange) || isstring(maxrange)
        validatestring(maxrange, {'auto'}, 'coverage', 'MaxRange');
        
        % Get max range for each tx
        if pm.requiresTerrain || pm.isMultipathModel
            % Use default value for terrain prop model
            maxrange = rfprop.internal.Validators.defaultMaxRange(txslatlon, pm, map);
        else
            maxrange = zeros(1,numSites);
            ss = min(strengths);
            
            % Turn off warning on calling range since check is made below
            warnState = warning('off','shared_channel:rfprop:RangeGreaterThanMax');
            warnCleanup = onCleanup(@()warning(warnState));
            for k = 1:numSites
                maxrange(k) = range(txs(k), ss, rxGain(k), type, pm);
            end
            
            % Get propagation limit
            if isa(map,'siteviewer') 
                useTerrain = map.UseTerrain;
            else
                useTerrain = ~strcmp(map,'none');
            end
            if useTerrain
                maxRangeLimit = rfprop.Constants.MaxPropagationDistanceUsingTerrain;
                warnID = 'shared_channel:rfprop:MaxRangeAssignedTerrainLimit';
            else
                maxRangeLimit = rfprop.Constants.MaxPropagationDistance;
                warnID = 'shared_channel:rfprop:MaxRangeAssignedPropagationLimit';
            end
            
            % Saturate and throw warning if range is at limit. Use limit
            % equality instead of strict greater than check since "range"
            % above may have already saturated to the limit.
            exceedsMax = (maxrange >= maxRangeLimit);
            if any(exceedsMax)
                maxrange(exceedsMax) = maxRangeLimit;
                warning(message(warnID,round(maxRangeLimit/1000)));
            end
        end
    else
        maxrange = rfprop.internal.Validators.validateNumericMaxRange(maxrange, pm, numSites, map, 'coverage');
    end
catch e
    throwAsCaller(e);
end
end

function isWithinBldg = isWithinAnyBuilding(bldgs, lats, lons)

% Check if locations are within footprint of or touching any buildings
numBldg = numel(bldgs);
inp = false(numel(lons),numBldg);
onp = inp;
for bldgInd = 1:numBldg
    thisVertices = bldgs(bldgInd).Footprint.Vertices;
    [inp(:,bldgInd), onp(:,bldgInd)] = inpolygon(lats,lons,thisVertices(:,1),thisVertices(:,2));
end
isWithinBldg = any((inp | onp), 2);
end

function setColorGraphics(viewer,imageGraphics)

viewer.ColorGraphics = imageGraphics;
end