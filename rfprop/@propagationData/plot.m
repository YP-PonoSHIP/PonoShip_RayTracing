function plot(pd,varargin)
%

% Copyright 2019-2024 The MathWorks, Inc.

validateattributes(pd,{'propagationData'},{'scalar'},'propagationData/plot','',1);

% Process optional name/value pairs
p = inputParser;
p.addParameter('Animation', '');
p.addParameter('EnableWindowLaunch', true);
p.addParameter('DataVariableName', pd.DataVariableName);
p.addParameter('Type', 'custom');
p.addParameter('Levels', []);
p.addParameter('Colors', []);
p.addParameter('ColorLimits', []);
p.addParameter('Colormap', 'jet');
p.addParameter('MarkerSize', 10);
p.addParameter('ShowLegend', true);
p.addParameter('LegendTitle', '');
p.addParameter('Map', []);
p.parse(varargin{:});

% Get Site Viewer and control parameters
viewer = rfprop.internal.Validators.validateMap(p, 'plot');
isViewerInitiallyVisible = viewer.Visible;
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(...
    p, isViewerInitiallyVisible, 'plot');

% Get data variable name and value
if ismember('DataVariableName',p.UsingDefaults)
    dataVariableName = pd.DataVariableName;
else
    dataVariableName = validateDataVariableName(pd, p.Results.DataVariableName, 'plot');
end
[data,lats,lons,alts] = pd.getDataVariable(dataVariableName);

% Validate Type and and get corresponding default parameter values
allowedTypes = {'power','efield','sinr','pathloss','custom'};
[type, levels, defaultColorLimits, defaultLegendTitle] = ...
    rfprop.internal.Validators.validateType(p, allowedTypes, 'plot', data);

% Validate graphical parameters
colors = rfprop.internal.Validators.validateColors(p, 'plot');
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'plot');
cmap = rfprop.internal.Validators.validateColorMap(p, 'plot');
isPathlossPlot = strcmp(type,'pathloss');
if isPathlossPlot
    % Flip colormap to maintain consistency of color meaning (red end of
    % jet means good signal), since path loss is inversely proportional to
    % received power
    cmap = flipud(cmap);
end
markerSize = validateMarkerSize(p);
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'plot');
legendTitle = rfprop.internal.Validators.validateLegendTitle(p, defaultLegendTitle, 'plot');

% Compute data/legend colors, removing data below the lowest level
cdata = pd.discretizeDataToLevels(data, levels);
isBelowLevels = isnan(cdata);
cdata(isBelowLevels) = [];
data(isBelowLevels) = [];
lats(isBelowLevels) = [];
lons(isBelowLevels) = [];
alts(isBelowLevels) = [];
useColors = ~isempty(colors);
[markerRGB, legendColors, legendColorValues] = rfprop.internal.ColorUtils.dataColors(...
    cdata, useColors, colors, cmap, clim, levels, showLegend);
if isPathlossPlot
    legendColors = fliplr(legendColors);
    legendColorValues = fliplr(legendColorValues);
end

if useColors
    colorData = struct('Levels',levels,'Colors',colors);
else
    colorData = struct('Colormap',cmap,'ColorLimits',clim);
end

% Compute plot locations for data, where elevation (Z) must match terrain
% source height reference
lats = double(lats);
lons = wrapTo180(double(lons));
switch pd.AltitudeReference
    case 'geoid'
        Z = alts;
        if viewer.UseTerrain && strcmpi(viewer.TerrainSource.HeightReference,'ellipsoid')
            Z = terrain.internal.HeightTransformation.orthometricToEllipsoidal(Z, lats, lons);
        end
    case 'ellipsoid'
        Z = alts;
        if viewer.UseTerrain && strcmpi(viewer.TerrainSource.HeightReference,'geoid')
            Z = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(Z, lats, lons);
        end
    case 'ground'
        Z = alts + rfprop.internal.AntennaSiteCoordinates.queryGroundHeightAboveTerrainReference(...
            lats, lons, viewer);
    case 'surface'
        Z = alts + rfprop.internal.AntennaSiteCoordinates.querySurfaceHeightAboveTerrainReference(...
            lats, lons, viewer);
end
plotLocations = [lats lons Z];
pdID = pd.pUID;
% Clear old plot color info for the pd. This is because a propData can
% never conflict with itself (since a propData will replot it's own points if
% plot is called on the same propData multiple times)
% Also mark the legend and old plot graphics to be removed
oldPlotID = viewer.getGraphic(pdID, 'plot');
graphicsToRemove = {};
if (~isempty(oldPlotID))
    viewer.removeColorData(oldPlotID);
    graphicsToRemove = [graphicsToRemove; oldPlotID];
end
if ~isempty(viewer.LegendID)
    graphicsToRemove = [graphicsToRemove; viewer.LegendID];
    viewer.LegendID = '';
    viewer.LegendEntities = strings(0);
end

% Validate the plot against other plots and contours
viewer.checkForGraphicsConflict(type, pd.plotID, colorData);

viewer.remove(graphicsToRemove);

% Add legend
composite = globe.internal.CompositeModel;
if showLegend
    lv = globe.internal.LegendViewer;
    legendID = ['legend' pdID];
    [~, legendDescriptor] = lv.buildPlotDescriptors(legendTitle, legendColors, legendColorValues, "ID", legendID);
    viewer.LegendID = legendID;
    composite.addGraphic("colorLegend", legendDescriptor);
end

% Compute descriptions
descriptionTitle = pd.Name;
numData = numel(data);
descriptions = cell(numData,1);
dec = rfprop.Constants.MaxLatLonInfoDecimalPlaces;
for k=1:numData
    descriptions{k} = [mat2str(round(lats(k), dec)) ', ' mat2str(round(lons(k), dec)) ...
        '<br />' dataVariableName ': ' mat2str(round(data(k), dec))];
end

% Add point markers
[~, pointDescriptor] = viewer.Instance.GlobeViewer.getViewer('Point').buildPlotDescriptors(plotLocations, ...
    'Color', num2cell(markerRGB(:,:),2), ...
    'Size', markerSize, ...
    'Animation', animation, ...
    'Name', descriptionTitle, ...
    'ID', {pd.plotID}, ...
    'Description', descriptions, ...
    'EnableWindowLaunch', enableWindowLaunch);
composite.addGraphic('point',pointDescriptor);
% Set the sync flag to true for the point plot to tell the JS code that
% when it clears the graphic, it can avoid manual clearing of entities to
% work around a 3p performance bug
composite.PlotDescriptors{end}.WaitForResponse = true;
compositeController = globe.internal.CompositeController(viewer.Instance.GlobeViewer.Controller);

% Zoom if the viewer isn't open yet
compositeArgs = composite.buildPlotDescriptors;
if ~isViewerInitiallyVisible && (isempty(animation) || strcmp(animation, 'zoom'))
    compositeArgs.Animation = 'zoom';
end

% Plot composite
compositeController.composite(compositeArgs)
viewer.SiteGraphics.(pdID).plot = pd.plotID;

% Show the viewer and bring it into focus
if ~isViewerInitiallyVisible && enableWindowLaunch
    viewer.Visible = true;
end
viewer.bringToFront;
end

function markerSize = validateMarkerSize(p)

try
    markerSize = p.Results.MarkerSize;
    validateattributes(markerSize, {'numeric'}, ...
        {'real','nonsparse','finite','scalar','positive','integer','<=',1000}, 'plot', 'MarkerSize');
catch e
    throwAsCaller(e);
end
end