function contour(pd, varargin)
%

% Copyright 2019-2024 The MathWorks, Inc.

validateattributes(pd,{'propagationData'},{'scalar'},'propagationData/contour','',1);

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
p.addParameter('ImageSize', []);
p.addParameter('Transparency', 0.4);
p.addParameter('ShowLegend', true);
p.addParameter('LegendTitle', '');
p.addParameter('Map', []);
p.addParameter('ValidateColorConflicts', true); % Allows us to skip color-conflict validation if coverage or sinr already did it
p.parse(varargin{:});

% Get Site Viewer and control parameters
viewer = rfprop.internal.Validators.validateMap(p, 'contour');
isViewerInitiallyVisible = viewer.Visible;
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(...
    p, isViewerInitiallyVisible, 'contour');

% Get data variable name and value
if ismember('DataVariableName',p.UsingDefaults)
    dataVariableName = pd.DataVariableName;
else
    dataVariableName = validateDataVariableName(pd, p.Results.DataVariableName, 'contour');
end
[data,lats,lons] = pd.getDataVariable(dataVariableName);

% Validate Type and and get corresponding default parameter values
allowedTypes = {'power','efield','sinr','pathloss','custom'};
[type, levels, defaultColorLimits, defaultLegendTitle] = ...
    rfprop.internal.Validators.validateType(p, allowedTypes, 'contour', data);

% Validate graphical parameters
colors = rfprop.internal.Validators.validateColors(p, 'contour');
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'contour');
cmap = rfprop.internal.Validators.validateColorMap(p, 'contour');
isPathlossPlot = strcmp(type,'pathloss');
if isPathlossPlot
    % Flip colormap to maintain consistency of color meaning (red end of
    % jet means good signal), since path loss is inversely proportional to
    % received power
    cmap = flipud(cmap);
end
transparency = rfprop.internal.Validators.validateTransparency(p, 'contour');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'contour');
legendTitle = rfprop.internal.Validators.validateLegendTitle(p, defaultLegendTitle, 'contour');
validateColorConflicts = p.Results.ValidateColorConflicts;

% Get color info
useColors = ~isempty(colors);
if useColors
    colorData = struct('Levels',levels,'Colors',colors);
else
    colorData = struct('Colormap',cmap,'ColorLimits',clim);
end

if (validateColorConflicts)
    oldContourID = viewer.getGraphic(pd.pUID, 'contour');
    if (~isempty(oldContourID))
        % Clear the old contour image info since a site can never
        % color-conflict with itself
        viewer.removeColorData(oldContourID);
    end
    
    % Validate image for color conflicts
    oldColorGraphics = viewer.ColorGraphics;
    viewer.checkForGraphicsConflict(type, pd.contourID, colorData);
    resetColorGraphics = rfprop.internal.onExit(@()setColorGraphics(viewer,oldColorGraphics));
end

% Get ImageSize from stored value or else from parameter
imageSize = pd.pImageSize;
if isempty(imageSize)
    imageSize = validateImageSize(p, viewer, data);
end

% Generate lat/lon grid for contour map image
[lonmin,lonmax] = bounds(lons);
[latmin,latmax] = bounds(lats);
imlonsv = linspace(lonmin,lonmax,imageSize(2));
imlatsv = linspace(latmin,latmax,imageSize(1));
[imlons,imlats] = meshgrid(imlonsv,imlatsv);
imlons = imlons(:);
imlats = imlats(:);

% Trim data to that within range (if stored on object)
maxrange = pd.pTransmitterMaxRange;
if ~isempty(maxrange)
    % Get tx locations
    txlat = pd.pTransmitterLatitude;
    txlon = pd.pTransmitterLongitude;

    % Find data locations that are within range of any transmitter site
    numTx = size(txlat,1);
    gc = nan(numel(imlons),numTx);
    for txInd = 1:numTx
        gc(:,txInd) = rfprop.internal.MapUtils.greatCircleDistance(...
            txlat(txInd), txlon(txInd), imlats, imlons);
    end
    isInRange = any(gc <= maxrange,2);
else
    isInRange = true(size(imlats));
end

% Create image data grid using interpolation. Query for locations within
% range of any transmitter site, or otherwise use NaN.
cdata = nan(imageSize);
cdata(isInRange) = interp(pd,imlats(isInRange),imlons(isInRange), ...
    'DataVariableName',dataVariableName);
cdata = flipud(cdata); % Start data columns from north instead of south

% Discretize image data so that each value is replaced by the corresponding
% contour level value or else NaN if it is below the minimum value
cdata = pd.discretizeDataToLevels(cdata, levels);

% Return early if no data to show, which likely means the minimum data
% value specified cannot be met
if isempty(cdata) || all(isnan(cdata(:)))
    minlevel = min(levels,[],"all");
    [mindata, maxdata] = bounds(data,"all");
    warning(message("shared_channel:rfprop:ContourmapNoDataArea", ...
        num2str(minlevel),num2str(mindata),num2str(maxdata)))
    viewer.removeColorData(pd.contourID);
    return
end

% Create image data matrix, mapping color data to RGB values
[imageRGB, legendColors, legendColorValues] = rfprop.internal.ColorUtils.dataColors(...
    cdata, useColors, colors, cmap, clim, levels, showLegend);
if isPathlossPlot
    legendColors = fliplr(legendColors);
    legendColorValues = fliplr(legendColorValues);
end

% Create transparency matrix (make clear where data is NaN)
imageAlpha = ones(size(imageRGB,1),size(imageRGB,2));
imageAlpha(isnan(cdata)) = 0;

% Create temp image file and mark for deletion on close of Site Viewer
fileLoc = [tempname, '.png'];
imwrite(imageRGB, fileLoc, 'Alpha', imageAlpha);
viewer.addTempFile(fileLoc);

% Remove old plot and legend
pdID = pd.pUID;
graphicsToRemove = {};
graphicsIDs = getGraphic(viewer, pdID, 'contour');
if ~isempty(graphicsIDs)
    graphicsToRemove = [graphicsToRemove; graphicsIDs];
end
if ~isempty(viewer.LegendID)
    graphicsToRemove = [graphicsToRemove; viewer.LegendID];
    viewer.LegendID = '';
    viewer.LegendEntities = strings(0);
end
viewer.remove(graphicsToRemove);

% Add legend
composite = globe.internal.CompositeModel;
if showLegend
    lv = globe.internal.LegendViewer;
    legendID = ['legend' pd.contourID];
    [~, legendDescriptor] = lv.buildPlotDescriptors(legendTitle, legendColors, legendColorValues, "ID", legendID);
    viewer.LegendID = legendID;
    viewer.addToLegendEntities(pd.contourID);
    composite.addGraphic("colorLegend", legendDescriptor);
end

% Add image
iv = viewer.Instance.GlobeViewer.getViewer("Image");
[~, imageDescriptor] = iv.buildPlotDescriptors(fileLoc, [latmin, lonmin, latmax, lonmax], ...
    "Transparency", transparency, ...
    "ID", {pd.contourID}, ...
    "Animation", animation);
composite.addGraphic("image", imageDescriptor);

% Plot composite
compositeController = globe.internal.CompositeController(viewer.Instance.GlobeViewer.Controller);
compositeController.composite(composite.buildPlotDescriptors)
viewer.SiteGraphics.(pdID).contour = pd.contourID;

% Show the viewer and bring it into focus
if ~isViewerInitiallyVisible && enableWindowLaunch
    viewer.Visible = true;
end
viewer.bringToFront;
if (validateColorConflicts)
    cancel(resetColorGraphics);
end
end

function imageSize = validateImageSize(p, viewer, data)

try
    if ismember('ImageSize',p.UsingDefaults)
        % Scale default ImageSize with square root of data size, using a
        % minimum of [500 500]
        imageLength = min(max(500, round(sqrt(numel(data)))),viewer.MaxImageSize);
        imageSize = [imageLength imageLength];
    else
        imageSize = p.Results.ImageSize;
        validateattributes(imageSize,{'numeric'}, ...
            {'real','finite','nonnan','nonsparse','positive','row','ncols',2}, ...
            'contour', 'ImageSize');
    end
catch e
    throwAsCaller(e);
end
end


function setColorGraphics(viewer,imageGraphics)

viewer.ColorGraphics = imageGraphics;
end