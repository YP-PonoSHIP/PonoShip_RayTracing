function varargout = sinr(txs, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate number of output arguments
nargoutchk(0,1)

% Validate site
validateattributes(txs,{'txsite'},{'nonempty'},'sinr','',1);

% Add parameters
p = inputParser;
if mod(numel(varargin),2) % Odd number of inputs
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('SignalSource', 'strongest');
p.addParameter('Values', -5:20);
p.addParameter('Resolution', 'auto');
p.addParameter('ReceiverGain', 2.1);
p.addParameter('ReceiverAntennaHeight', 1);
p.addParameter('ReceiverNoisePower', -107);
p.addParameter('Animation', '');
p.addParameter('MaxRange', []);
p.addParameter('Colormap', 'jet');
p.addParameter('ColorLimits', []);
p.addParameter('Transparency', 0.4);
p.addParameter('ShowLegend', true);
p.addParameter('ReceiverLocationsLayout', []);
p.addParameter('MaxImageResolutionFactor', 5);
p.addParameter('RadialResolutionFactor', 2);
p.addParameter('Map', []);
p.parse(varargin{:});

% Validate CoordinateSystem
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(txs);

% Get Site Viewer and validate web graphics
outputRequested = nargout > 0;
if ~outputRequested
    viewer = rfprop.internal.Validators.validateMap(p, 'sinr', usingCartesian);
    isViewerInitiallyVisible = viewer.Visible;
else
    isViewerInitiallyVisible = false;
end

% Create vector array of all txs
sigSource = validateSignalSource(p);
txs = txs(:);
if isa(sigSource,'txsite')
    if ~ismember(sigSource,txs)
        txs = [txs; sigSource];
    end
end
numTx = numel(txs);

if usingCartesian
    error(message('shared_channel:rfprop:CoordinateSystemCartesianSINRNotSupported'));
end

% Validate parameters
[map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'sinr');
terrainSource = rfprop.internal.Validators.validateTerrainSource(map, 'sinr');
animation = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'sinr');
values = validateValues(p);
clim = rfprop.internal.Validators.validateColorLimits(p, [-5 20], 'sinr');
cmap = rfprop.internal.Validators.validateColorMap(p, 'sinr');
transparency = rfprop.internal.Validators.validateTransparency(p, 'sinr');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'sinr');
pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'sinr');
rfprop.internal.Validators.validateMaxNumReflections(pm, 'sinr');
rfprop.internal.Validators.validateMaxNumDiffractions(pm, 'sinr'); % Disable 2 order diffractions
rxGain = rfprop.internal.Validators.validateReceiverGain(p, 'sinr');
rxAntennaHeight = rfprop.internal.Validators.validateReceiverAntennaHeight(p, 'sinr');
noisePower = validateReceiverNoisePower(p);
maxImageResFactor = p.Results.MaxImageResolutionFactor;
radialResFactor = p.Results.RadialResolutionFactor;

% Get site coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, map);
txslatlon = txsCoords.LatitudeLongitude;

% Validate dependent parameters
if ismember('MaxRange', p.UsingDefaults)
    maxrange = rfprop.internal.Validators.defaultMaxRange(txslatlon, pm, map);
else
    maxrange = rfprop.internal.Validators.validateNumericMaxRange(p.Results.MaxRange, pm, numTx, map, 'sinr');
end
[res, isAutoRes] = rfprop.internal.Validators.validateResolution(p, maxrange, 'sinr');
datarange = rfprop.internal.Validators.validateDataRange(txslatlon, maxrange, res, ~strcmp(terrainSource,'none'));
rxLocationsLayout = rfprop.internal.Validators.validateReceiverLocationsLayout(p, pm, txslatlon, 'sinr');

if outputRequested
    % Do not show progress dialog
    generatingMapMsg = '';
    computingDataMsg = '';
else
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
    colorData = struct('Colormap',cmap,'ColorLimits',clim);
    viewer.checkForGraphicsConflict('sinr', contourID, colorData);
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
    generatingMapMsg = message('shared_channel:rfprop:ProgressDialogGeneratingSINRMap').getString;
    computingDataMsg = message('shared_channel:rfprop:ProgressDialogComputingSINR').getString;
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
    latNorth, latSouth, lonEast, lonWest, res, isAutoRes, maxrange, maxImageSize, 'sinr');
gridSize = size(gridlats);

% Compute and validate image size for SINR map
imageSize = rfprop.internal.Validators.validateImageSize(...
    gridSize, maxImageResFactor, maxImageSize, res, 'sinr');

% Trim grid locations to those which are within data range
[datalats, datalons] = rfprop.internal.MapUtils.georange(...
    txs, gridlats(:), gridlons(:), datarange, terrainSource);

if ~outputRequested && siteViewerWasClosed(viewer)
    return % Abort if Site Viewer has been closed
end
if isempty(datalats)
    warning(message('shared_channel:rfprop:NoSINRMapArea'))
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
    type = 'power';
    if strcmp(rxLocationsLayout,'grid')
        % Define rxsites at grid locations within data range
        rxs = rxsite(...
            'Name', 'internal.sinrsite', ... % Specify to avoid default site naming
            'Latitude', datalats, ...
            'Longitude', datalons, ...
            'AntennaHeight', rxAntennaHeight);
        
        % Launch waitbar phase of dialog
        if ~outputRequested && showWaitbar
            viewer.showProgressDialog('Message', computingDataMsg, ...
                'Value', 0, ...
                'ShowPercentage', true, ...
                'Indeterminate', false, ...
                'Cancelable', true);
        end
        
        % Compute signal strength at each rxsite in the grid
        ss = sigstrength(rxs, txs, pm, ...
            'Type', type, ...
            'ReceiverGain', rxGain, ...
            'Map', mapStruct, ...
            'TransmitterAntennaSiteCoordinates', txsCoords);
    else
        % Compute data and corresponding locations using radial layout
        [datalats, datalons, rxs, ss] = radialReceiverLocationsLayoutData(txs, txsCoords, type, ...
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

% Cache signal strength on site coordinates
txsCoords.addCustomData('SignalStrength', ss);

if (numTx > 1) && strcmp(rxLocationsLayout,'radial')
    % Compute SINR at each rxsite in data grid, which is the same size as
    % the image grid. Create a dummy rxsite to call the method, but the
    % data used to compute values is passed as SignalStrength.
    dummyrx = rxsite(...
        'Name', 'internal.sinrsite', ...
        'Latitude', datalats(1), ...
        'Longitude', datalons(1));
    
    data = sinr(dummyrx, txs, ...
        'SignalSource', sigSource, ...
        'ReceiverNoisePower', noisePower, ...
        'PropagationModel', pm, ...
        'ReceiverGain', rxGain, ...
        'TransmitterAntennaSiteCoordinates', txsCoords, ...
        'Map', map);
else
    % Compute SINR at each rxsite in the grid
    data = sinr(rxs, txs, ...
        'SignalSource', sigSource, ...
        'ReceiverNoisePower', noisePower, ...
        'PropagationModel', pm, ...
        'ReceiverGain', rxGain, ...
        'TransmitterAntennaSiteCoordinates', txsCoords, ...
        'Map', map);
end

if ~outputRequested && siteViewerWasClosed(viewer)
    return
end

% Remove NaN-data
isnn = ~isnan(data);
data = data(isnn);
datalats = datalats(isnn);
datalons = datalons(isnn);

% Create propagation data container
pd = propagationData(datalats, datalons, ...
    'Name', message('shared_channel:rfprop:SINRPropagationDataTitle').getString, ...
    'SINR', data(:));

if outputRequested
    % Return propagationData with contour properties set except for ID,
    % which is not set since there is no plot to associate with
    pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize);
    varargout = {pd};
    return
end

% It is now safe to remove old contours associated with sites
viewer.remove(graphicsToRemove);

% Return early if no SINR map data meets minimum SINR value
minvalue = min(values,[],"all");
[mindata, maxdata] = bounds(data,"all");
if maxdata < minvalue
    warning(message("shared_channel:rfprop:NoSINRMapAreaMeetsMinValue", ...
        num2str(minvalue),num2str(mindata),num2str(maxdata)))
    return
end

% Create contour map
pd = pd.setContourProperties(txslatlon(:,1),txslatlon(:,2),maxrange,imageSize,contourID);
contour(pd, ...
    'Type','sinr', ...
    'Animation', animation, ...
    'Map', viewer, ...
    'Levels', values, ...
    'ColorLimits', clim, ...
    'Colormap', cmap, ...
    'Transparency', transparency, ...
    'ShowLegend', showLegend, ...
    'ImageSize', imageSize, ...
    'ValidateColorConflicts', false);

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

function sigsource = validateSignalSource(p)

try
    sigsource = p.Results.SignalSource;
    if ischar(sigsource) || isstring(sigsource)
        sigsource = validatestring(sigsource, {'strongest'}, ...
            'sinr','SignalSource');
    else
        validateattributes(sigsource,{'txsite'}, {'scalar'}, ...
            'sinr','SignalSource');
    end
catch e
    throwAsCaller(e);
end
end

function values = validateValues(p)

try
    values = p.Results.Values;
    validateattributes(values, {'numeric'}, ...
        {'real','nonnan','nonsparse','vector','nonempty'}, 'sinr', 'Values');
    if ~iscolumn(values)
        values = values(:);
    end
    values = double(values);
catch e
    throwAsCaller(e);
end
end

function noisePower =  validateReceiverNoisePower(p)

try
    noisePower = p.Results.ReceiverNoisePower;
    validateattributes(noisePower, {'numeric'}, {'real','finite','nonnan','nonsparse','scalar'}, ...
        'sinr', 'ReceiverNoisePower');
catch e
    throwAsCaller(e);
end
end

function setColorGraphics(viewer,imageGraphics)

viewer.ColorGraphics = imageGraphics;
end