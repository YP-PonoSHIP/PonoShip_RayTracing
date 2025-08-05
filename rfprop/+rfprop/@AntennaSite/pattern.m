function pattern (site, varargin)
%

% Copyright 2018-2024 The MathWorks, Inc.

% Ensure that the input site is scalar
validateattributes(site, {'txsite', 'rxsite'}, {'scalar'}, 'pattern', '', 1);

% Special case for Rxsites: The second input must be a frequency
p = inputParser;
if (isa(site, 'rxsite'))
    p.addRequired('Frequency');
    % Explicitly check for a frequency value
    if isempty(varargin)
        error(message('MATLAB:InputParser:notEnoughInputs'));
    elseif ~isnumeric(varargin{1}) && numel(varargin) > 1
        % If more than 1 additional input is specified, then varargin
        % corresponds to name-value pairs in the absence of the required
        % frequency input.
        error(message('shared_channel:rfprop:ReceiverPatternFrequencyRequired'));
        % If exactly 1 additional input is specified, then it must be
        % checked for validity as a frequency.
    end
else
    p.addOptional('Frequency', site.TransmitterFrequency);
end
% Parse Inputs
p.addParameter('Animation', '');
p.addParameter('Size', "auto");
p.addParameter('Transparency', .4);
p.addParameter('Colormap', jet);
p.addParameter('ColorLimits', []);
p.addParameter('Map', []);
p.addParameter('Resolution', "high");
p.parse(varargin{:});

% Validate and assign parsed inputs
isCartesian = strcmp(site.CoordinateSystem, 'cartesian');
patternSize = validateSize(p, isCartesian);
transparency = rfprop.internal.Validators.validateTransparency(p, 'pattern');
cmap = rfprop.internal.Validators.validateColorMap(p, 'pattern');
colorLimits = rfprop.internal.Validators.validateColorLimits(p, [-1,1], 'pattern');
frequency = validateFrequency(p);

% Get Site Viewer visibility state and validate web graphics
viewer = rfprop.internal.Validators.validateMap(p, 'pattern', isCartesian);
if isCartesian
    siteCoordSys = 'cartesian';
else
    siteCoordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(viewer, siteCoordSys);
isViewerInitiallyVisible = viewer.Visible;

% Validate viewer-dependent inputs
coords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(site, viewer);
if isCartesian
    % The size of the pattern model will always be scaled based on 1 unit
    % for cartesian cases. The JavaScript code will take care of scaling it
    % based on the model or input size.
    patternSizeValue = 1;
else
    % For geographic cases, the size of the model is just the specified
    % size.
    patternSizeValue = patternSize;
end
[animation, ~] = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'pattern');
ant = site.Antenna;

isElectromagneticAntenna = rfprop.AntennaSite.isElectromagneticAntenna(ant);
resolution = validateResolution(p);

if (isElectromagneticAntenna)
    % The antenna is from Antenna Toolbox
    resolutionValue = rfprop.Constants.AntennaPatternResolution.(resolution);
    interpRes = rfprop.Constants.AntennaPatternInterpolation.(resolution);
elseif ischar(ant)
    % The antenna is the default 'isotropic' string
    resolutionValue = rfprop.Constants.IsotropicPatternResolution.(resolution);
    interpRes = resolutionValue;
else
    % The antenna is a phased antenna or an arrayConfig object (which
    % behaves like a phased antenna)
    resolutionValue = rfprop.Constants.PhasedPatternResolution.(resolution);
    interpRes = resolutionValue;
end

% Convert pattern coordinates for mesh X, Y, Z
% Call pattern
if ischar(ant)
    % Implement isotropic antenna
    pAz = -180:resolutionValue:180;
    pEl = -90:resolutionValue:90;
    r = zeros(numel(pEl),numel(pAz));
else
    % For antenna toolbox antennas using "medium" or "high" resolution, we
    % use 'auto' resolution which defaults to the pattern function's
    % resolution
    if (strcmp(resolutionValue, 'auto'))
        % Use default angular resolution of pattern, which is antenna-dependent
        [r,pAz,pEl] = pattern(ant, frequency);
    else
        % Use custom resolution
        pAz = -180:resolutionValue:180;
        pEl = -90:resolutionValue:90;
        [r, pAz, pEl] = pattern(ant, frequency, pAz, pEl);
    end
end

%Phased array uses a different way to scale the directivities so
%we must account for that.
if (rfprop.AntennaSite.isPhasedAntenna(ant)) || ...
   (rfprop.AntennaSite.isCommAntenna(ant))
    dBRange = 50; % find out how to calculate dBRange
    r =  limitDynamicdBRange(r,dBRange);
end

[mAz, mEl] = meshgrid(pAz, pEl);
%Find the resolution of the figure pattern function
patternRes = 360 / (length(pAz) - 1);

%If the resolution of the figure pattern is worse than 1, then we will
%interpolate the data to achieve 1 degree resolution
if (interpRes < patternRes)
    % Interpolate the azimuth, elevation, and calculated directivities.
    queryAz = -180:interpRes:180;
    queryEl = -90:interpRes:90;
    [queryAzMesh, queryElMesh] = meshgrid(queryAz, queryEl);
    r = interp2(mAz, mEl, r, queryAzMesh, queryElMesh);
    az = queryAz;
    el = queryEl;
else
    % No need to interpolate
    az = pAz;
    el = pEl;
end


[X, Y, Z, r2] = formatDataForPlot(az, el, r);

% Calculate the scale of the pattern given patternSize
[~, maxInd] = max(r2(:));
maxXYZ = [X(maxInd), Y(maxInd), Z(maxInd)];
distToBoresight = sqrt(maxXYZ(1).^2 + maxXYZ(2).^2 + maxXYZ(3).^2);
scaleFactor = patternSizeValue / distToBoresight;

% Triangulate the data using a triangle-based patch
[tri, xyzData] = surf2patch(X,Y,Z,'triangles');

% Map directivity values to color values
rVec = r2(:);
minr = min(rVec);
maxr = max(rVec);
% Normalize the directivity
normR = rVec + abs(minr);
normR = normR/maxr;
minColor = colorLimits(1);
maxColor = colorLimits(2);
% Set the values to correspond to the min/max colors
diff = maxColor - minColor;
normR = normR * diff;
normR = normR + minColor;
% Index into the directivities using the triangulation
% To get the color of each triangle.
rColorsTriangles = normR(tri);
rColorsTriangles = rColorsTriangles';
CData = rColorsTriangles;

% Convert color data to RGB triplets 
cmin = min(CData(:));
cmax = max(CData(:));
% If the antenna's gain values are all the same
% Then we must make the range go from -1 to 1 around the value.
% This allows us to generate a valid set of colors for our pattern.
if (cmin == cmax)
    cmin = cmin - 1;
    cmax = cmax + 1;
end
m = size(cmap,1);
colormapInd = fix((CData-cmin)/(cmax-cmin)*m) + 1;
colormapInd(colormapInd>m) = m;
colormapInd(colormapInd<1) = 1;
imgRGB = ind2rgb(colormapInd, cmap);

rmin = min(min(r));
rmax = max(max(r));
[legendColors, legendColorValues] = imageRGBData(cmap, [rmin, rmax]);
legendTitle = message('shared_channel:rfprop:RadiationPatternLegendTitle').getString;

if isViewerInitiallyVisible && ~viewer.Visible
    return % Abort if Site Viewer has been closed (test before show)
end

% Queue plots before calling show to synchronize graphics
viewer.Visualizer.queuePlots();
resetUnqueuePlots = rfprop.internal.onExit(@()unqueuePlots(viewer));
show(site,'AntennaSiteCoordinates',coords,'Map',viewer,'EnableWindowLaunch', false)

% Remove graphics first in our plot queue to ensure globe is clean before
% we plot
viewer.remove(viewer.getGraphic(site.UID, 'pattern'));
viewer.remove(viewer.getGraphic(site.UID, 'infoboxLegend'));

% Transpose triangles for easier use in JS
tri = tri';
indices = tri(:);

% Scale pattern
xyzData = xyzData * scaleFactor;

% reshape the colors to a M*N matrix, where M is the number of
% vertices, N is 3 (RGB channels)
numIndices = numel(indices);
numVertices = size(xyzData, 1);
imgRGB = reshape(imgRGB, numIndices, 3); 
vColors = zeros(numVertices, 3);
for i = 1:numIndices
    idx = indices(i);
    color = imgRGB(i,:);
    vColors(idx,:) = color;
end
% Convert color to linear RGB space
vColors = rfprop.internal.ColorUtils.srgb2lin(vColors);

% Get rid of NaN values in xyzData
tempxyzData = fillmissing(xyzData, 'constant', 0);

% Create a triangulation object from the temp data which has NaNs removed
T = reshape(indices, 3, [])';
TR = triangulation(T, tempxyzData);

% Convert the triangulation to a generic struct using the original data
% containing NaN values. This allows NaN values in patterns to appear as
% empty spaces instead of being filled in.
TR = struct("Points", xyzData, "ConnectivityList", TR.ConnectivityList, ...
    "faceNormal", TR.faceNormal);

patternID = viewer.getId(1);
patternID = ['pattern' num2str(patternID{1})];
viewer.SiteGraphics.(site.UID).pattern = patternID;

if size(site.AntennaAngle, 1) == 2
    modelAngle = [site.AntennaAngle;0];
else
    modelAngle = [site.AntennaAngle;0;0];
end

% Create and show 3D model for pattern
patternModel = globe.internal.Geographic3DModel(TR, ...
    'VertexColor', vColors, ...
    'YUpCoordinate', ~viewer.IsCartesian, ...
    'EnableLighting', false);

% Turn clipping on
viewer.turnClippingOn();

% Create an infobox color legend
idNum = viewer.getId(1);
legendID = ['infoboxLegend' num2str(idNum{1})];
viewer.Visualizer.legend(legendTitle, legendColors, legendColorValues, ...
    "ParentGraphicID", viewer.SiteGraphics.(site.UID).marker.ID, ...
    "ID", legendID, ...
    "InfoboxLegend", true);

viewer.SiteGraphics.(site.UID).infoboxLegend = legendID;

if (~isViewerInitiallyVisible && (isempty(animation) || strcmp(animation, 'zoom')))
    submitAnim = 'zoom';
else
    submitAnim = 'fly';
end
if viewer.IsCartesian
    location = coords.AntennaPosition;
    viewer.Visualizer.model3D(patternModel, location, ...
        'ID', patternID, ...
        'Transparency', transparency, ...
        'Rotation', [modelAngle(3) -modelAngle(2), modelAngle(1)], ...
        'Size', patternSize, ...
        'Animation', submitAnim);
else
    location = [coords.LatitudeLongitude(1), coords.LatitudeLongitude(2), coords.AntennaHeightAboveTerrainReference];
    viewer.Visualizer.geoModel3D(patternModel, location, ...
        'Animation', submitAnim, ...
        'Persistent', false, ...
        'Transparency', transparency,...
        'BoundingSphereRadius', patternSize,...
        'ID', patternID, ...
        'ShowIn2D', true, ...
        'Rotation', [modelAngle(1), modelAngle(2) modelAngle(3)], ...
        'ForcePlot', true);
end


viewer.Visualizer.submitPlots('Animation', 'none');

viewer.Visible = true;
viewer.bringToFront;
cancel(resetUnqueuePlots);

% Delete the generated 3D model file
try
    % Also disable permissions-related warnings (g3052830)
    warnState = warning('off', 'MATLAB:DELETE:Permission');
    warnCleanup = onCleanup(@()warning(warnState));
    delete(patternModel.File);
catch
    % Silently proceed if permissions are denied
end
end

% Validation Functions
function patternSize = validateSize(p, isCartesian)
try
    patternSize = p.Results.Size;
    if strcmp(patternSize, "auto")
        if ~isCartesian
            patternSize = 50;
        end
    else
        validateattributes(patternSize, {'numeric'}, {'nonnan','scalar','positive', ...
            'real','finite','<=',rfprop.Constants.MaxPropagationDistance}, 'pattern','Size');
    end
catch e
    throwAsCaller(e);
end
end

function resolution = validateResolution(p)
try
    resolution = validatestring(p.Results.Resolution, ["low", "medium", "high"]);
catch e
    throwAsCaller(e);
end
end

function fq = validateFrequency(p)
try
    % Validate frequency type and shape but allow antenna pattern function
    % to validate range
    fq = p.Results.Frequency;
    validateattributes(fq, {'numeric'}, {'nonnan','scalar','positive', ...
        'real','finite'}, 'pattern','Frequency');
catch e
    throwAsCaller(e);
end
end

function [legendColors, legendColorValues] = imageRGBData(cmap, clim)
% Return image RGB and alpha matrices

% Use Colors if user specified or if single signal level and no Colormap
% specified
% Convert image to indices in colormap using color limits. Code adapted
% from Tip on caxis reference page
cmin = clim(1);
cmax = clim(2);
m = size(cmap, 1);

% Compute legend values. The levels in the legend are generated using
% ColorLimits. The strategy is to show the color limits and
% intermediate values at a fixed step size, where the step size is
% chosen to be intuitive and generate at least 10 levels. A default
% step size of 10 is used. If this fails to produce at least 10 levels,
% then a step size of 5 is used. If that also fails to produce at least
% 10 levels, a step size of 2 is used. If that also fails, a final step
% size of 1 is used.

% If the antenna's gain values are all the same, then we will show a range
% of +- 1 around that gain value in the legend.
if (cmin == cmax)
    % Specify non-integer colors for the legend
    % Since the legend will be from -1 to 1 around cmin or cmax
    integerColors = false;
    roundedColor = round(cmin);
    cmin = roundedColor - 1;
    cmax = roundedColor + 1;
    colorMax = cmax;
    colorMin = cmin;
else
    % Get rounded color limits
    % Specify integer colors for the legend
    integerColors = true;
    colorMax = ceil(cmax);
    colorMin = floor(cmin);
end
% Compute color strengths, ensuring that max is included
% Interpolate to ensure that at least 64 boxes are used.
colorStrengths = linspace(colorMin, colorMax, 64);
if colorStrengths(end) ~= colorMax
    colorStrengths = [colorStrengths, colorMax];
end

% Grow legend values
numColorStrengths = numel(colorStrengths);
legendColors = strings(numColorStrengths, 1);
legendColorValues = strings(numColorStrengths, 1);
for strengthInd = 1:numColorStrengths
    colorStrength = colorStrengths(strengthInd);
    legendInd = fix((colorStrength-cmin)/(cmax - cmin)*m)+1;
    rgb = ind2rgb(legendInd, cmap);
    legendColors(strengthInd) = rfprop.internal.ColorUtils.rgb2css(rgb(:));
    if (integerColors)
        legendColorValues(strengthInd) = mat2str(round(colorStrength));
    else
        legendColorValues(strengthInd) = mat2str(round(colorStrength,1));
    end
end
legendColors = fliplr(legendColors);
legendColorValues = fliplr(legendColorValues);

end

function [X,Y,Z,r2] = formatDataForPlot(az, el, r)
phi1 = az;
theta1 = 90-el;
[theta, phi] = meshgrid(theta1, phi1);
MagE1 = r';
[offset, maxE1] = bounds(MagE1(:));
rRange = maxE1 - offset;
if (rRange == 0)
    % If all the gain values are the same, then it means we have a
    % perfect sphere pattern. This means the size of the pattern
    % doesn't matter, it just has to have a constant radius.
    % Phased pattern sets the offset to 1, so we will do the same.
    offset = 1;
end
MagE = reshape(MagE1,length(phi1),length(theta1));
r2 = MagE - offset;
[X, Y, Z] = psph2cart(phi, theta, r2./max(max(r2)));
end

function [X, Y, Z]= psph2cart(phi, theta, r)

Z  = r.*cosd(theta);
X  = r.*sind(theta).*cosd(phi);
Y  = r.*sind(theta).*sind(phi);
end

function dBRespLimited = limitDynamicdBRange(dBResp,dRange)
% From the
% phased.internal.AbstractRadiationPattern3D.limitDynamicdBRange function 
% when dRange is not a cell

respmax = max(dBResp, [], 'all');  % maximum response value in dB
respmin = respmax - dRange;  % minimum response value in dB
dBRespLimited = dBResp;
dBRespLimited(dBResp<respmin) = respmin;
end

function unqueuePlots(viewer)
    viewer.Visualizer.unqueuePlots();
end