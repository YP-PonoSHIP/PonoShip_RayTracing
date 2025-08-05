function varargout = los(observerSite, targetSites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate sites
validateattributes(observerSite,{'txsite' 'rxsite'},{'scalar'},'los','',1);
validateattributes(targetSites,{'txsite' 'rxsite'},{'nonempty'},'los','',2);
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(observerSite, targetSites);

% Check if too many sites to compute or show
outputRequested = nargout > 0;
numTargetSites = numel(targetSites);
if ~outputRequested && (numTargetSites >= rfprop.Constants.MaxNumSitesShow)
    error(message('shared_channel:rfprop:ShowTooManySites', rfprop.Constants.MaxNumSitesShow));
end

% Validate observer is not also a target
if strcmp(class(observerSite),class(targetSites)) && ismember(observerSite,targetSites)
    error(message('shared_channel:rfprop:LOSRepeatedSite', observerSite.Name));
end

% Process optional name/value pairs
p = inputParser;
p.addParameter('Animation', '');
p.addParameter('EnableWindowLaunch', true);
p.addParameter('Resolution','auto');
p.addParameter('VisibleColor', 'green');
p.addParameter('ObstructedColor', 'red');
p.addParameter('Map', []);
p.parse(varargin{:});

% Get Site Viewer visibility state and validate web graphics
if outputRequested
    nargoutchk(1,1)
    if usingCartesian
        map = rfprop.internal.Validators.validateCartesianMap(p);
    else
        map = rfprop.internal.Validators.validateMapTerrainSource(p, 'los');
    end
    isViewerInitiallyVisible = false;
else
    map = rfprop.internal.Validators.validateMap(p, 'los', usingCartesian);
    isViewerInitiallyVisible = map.Visible;
end
if usingCartesian
    siteCoordSys = 'cartesian';
else
    siteCoordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, siteCoordSys);
% Get site antenna coordinates
observerCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(observerSite, map);
targetsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(targetSites, map);

% Allocate output matrix
status = zeros(numTargetSites, 1, 'logical');

% Validate and get parameters
[animation, ~] = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'los');
visibleColor = rfprop.internal.Validators.validateLineColor(p, 'VisibleColor', 'los');
obstructedColor = rfprop.internal.Validators.validateLineColor(p, 'ObstructedColor', 'los');

% Cartesian line-of-sight, which only supports simple logical output
if usingCartesian
    if outputRequested
        triMap = map;
    else
        if ~isempty(map.ModelMesh)
            warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
            warnCleanup = onCleanup(@()warning(warnState));
            triMap = Geometry.meshToTriangulation(map.ModelMesh);
        else
            triMap = 'none';
        end
    end
    if strcmp(triMap,'none')
        % Trivial case: all sites are visible
        [status(:)] = true;
    else
        % Use ray tracing to analyze los
        env = struct( ...
            "Triangulation", triMap, ...
            "RayTracer", matlabshared.internal.StaticSceneRayTracer(triMap), ...
            "SharpEdgeFlags", comm.internal.geometry.getSharpEdges(triMap));    
        [firstPt,~,isobstructed] = comm.internal.geometry.firstIntersect( ...
            env, observerCoords.AntennaPosition, targetsCoords.AntennaPosition, 'segment');
        status = ~isobstructed;
    end
    if outputRequested
        varargout = {status};
        return
    else
        % Plot in cartesian siteviewer
        map.Visualizer.queuePlots;
        % Show the sites
        show(observerSite,'AntennaSiteCoordinates',observerCoords,'Map',map, ...
            'Animation','none','EnableWindowLaunch',false);
        show(targetSites,'AntennaSiteCoordinates',targetsCoords,'Map',map, ...
            'Animation','none','EnableWindowLaunch',false);
        for k = 1:numTargetSites
            % Associate and dissasociate graphics IDs with the sites
            observerID = observerSite.UID;
            targetID = targetSites(k).UID;
            graphicsToRemove = {};
            oldLOSLines = map.getGraphic(observerID, 'los');
            if (isfield(oldLOSLines, targetID))
                graphicsToRemove = [graphicsToRemove; oldLOSLines.(targetID)];
            end
            map.disassociateSiteGraphics(targetID, observerID, 'los');

            % Remove existing rays between the sites
            oldRayLines = map.getGraphic(observerID, 'rays');
            if (isfield(oldRayLines, targetID))
                map.removeFromLegendEntities(oldRayLines.(targetID));
                graphicsToRemove = [graphicsToRemove; oldRayLines.(targetID)];
            end
            map.disassociateSiteGraphics(targetID, observerID, 'rays');
            
            % Remove stale graphics
            map.remove(graphicsToRemove);
            
            % Get new IDs for the LOS lines
            losID = "los" + map.getId(1);
            
            % Plot the LOS lines
            targetCoords = targetsCoords.extract(k);
            if status(k)
                map.Visualizer.line([observerCoords.AntennaPosition; targetCoords.AntennaPosition], ...
                    "Color", visibleColor, "ID", losID, ...
                    "Name", message('shared_channel:rfprop:LOSName').getString, ...
                    "Description", message('shared_channel:rfprop:LOSDescriptionVisibleStatus').getString);
                map.associateSiteGraphics(targetID, observerID, 'los', losID);
            else
                map.Visualizer.line([observerCoords.AntennaPosition; firstPt(k, :)], ...
                    "Color", visibleColor, "DepthTest", true, "ID", losID, ...
                    "Name", message('shared_channel:rfprop:LOSName').getString, ...
                    "Description", message('shared_channel:rfprop:LOSDescriptionNotVisibleStatus').getString);
                losID2 = "los" + losID;
                map.Visualizer.line([firstPt(k, :); targetCoords.AntennaPosition], ...
                    "Color", obstructedColor, "DepthTest", true, "ID", losID2, ...
                    "Name", message('shared_channel:rfprop:LOSName').getString, ...
                    "Description", message('shared_channel:rfprop:LOSDescriptionNotVisibleStatus').getString);
                map.associateSiteGraphics(targetID, observerID, 'los', {losID; losID2});
            end
        end
        if (~isViewerInitiallyVisible && isempty(animation))
            animation = 'zoom';
        else
            animation = 'fly';
        end
        map.Visualizer.submitPlots("Animation", animation);
        if (~map.Visible)
            map.Visible = true;
        end
        map.bringToFront;
        return;
    end
end

res = validateResolution(p);
%Defining another variable for auto_resolution case handling
res_type = strcmpi(res,'auto');

% Compute visibility from observer to all targets
startIDs = cell(1,numTargetSites);
endIDs = startIDs;
visprofile = startIDs;
longitudes = startIDs;
latitudes = startIDs;
elevations = startIDs;
lineInfo = startIDs;
showObstructionPoint = startIDs;
dsGc = distance(observerSite, targetSites, 'greatcircle', 'Map', map, ...
    'SourceAntennaSiteCoordinates', observerCoords, ...
    'TargetAntennaSiteCoordinates', targetsCoords);
[dsEuc, ~, criticalAngles] = distanceangle(observerSite, targetSites, 'Map', map, ...
    'SourceAntennaSiteCoordinates', observerCoords, ...
    'TargetAntennaSiteCoordinates', targetsCoords);
sceneTri = [];
for k = 1:numTargetSites
    targetSite = targetSites(k);
    targetCoords = targetsCoords.extract(k);
    
    % Get antenna elevation and distances
    targetAntennaElevation = targetCoords.AntennaHeightAboveTerrainReference;
    dGc = dsGc(k);
    dEuc = dsEuc(k);
    
    % Track whther to show obstruction at all (i.e. line with green and red
    % segments, and also whether obstruction is from buildings analysis
    % (in which case a point is plotted)
    canShowObstruction = false;
    useBuildingsObstruction = false;
    
    maxTerrainSampleDist = rfprop.Constants.MaxLOSTerrainSampleDistance;
    if isnan(dGc) || (dGc > maxTerrainSampleDist)
        % Use fast calculation if distance between sites is very large.
        % This approach is not preferred because it is approximate and does
        % not calculate an intersection point.
        targetIsVisible = computeLOSUsingSpheroid(observerCoords, targetCoords);
    else
        % Check if sites are within buildings region
        observerLatLon = observerCoords.LatitudeLongitude;
        observerLat = observerLatLon(1);
        observerLon = observerLatLon(2);
        targetLatLon = targetCoords.LatitudeLongitude;
        endlats = [observerLat; targetLatLon(1)];
        endlons = [observerLon; targetLatLon(2)];
        withinBuildingsLimits = all(targetCoords.withinBuildingsLimits(endlats, endlons));
        
        if dGc == 0
            % Distance between sites is 0, so there must be visibility
            targetIsVisible = true;
        elseif withinBuildingsLimits
            % Compute LOS using ray tracing
            % Pass in the sceneTri in case multiple sites need to compute
            % los using raytracing. This avoids recomputing sceneTri every
            % iteration.
            [targetIsVisible, obstructingLocation, obstructingElevation, sceneTri] = ...
                computeLOSUsingRayTracing(observerCoords, targetCoords, map, sceneTri);
            canShowObstruction = true;
            useBuildingsObstruction = true;
        else
            % Compute resolution
            if res_type
                res = autoResolution(dGc);
            elseif res > dGc
                dGeoKmFormat = mat2str(round(dGc/1000, rfprop.Constants.MaxDistanceInfoDecimalPlaces));
                error(message('shared_channel:rfprop:LOSResolutionTooGreat', ...
                    observerSite.Name, targetSite.Name, dGeoKmFormat));
            end
            
            % Create array of sites along terrain profile, removing first
            % and last locations returned from great circle because those
            % correspond to observer and target sites
            [lats,lons] = sampleGreatCircle(observerSite, targetSite, res, 'Map', map, ...
                'SourceAntennaSiteCoordinates', observerCoords, ...
                'TargetAntennaSiteCoordinates', targetCoords);
            
            % Compute LOS using terrain profile. Only include profile
            % points which are outside of the buildings region, since a
            % separate LOS check will be made for the buildings region.
            profileWithinBuildingsLimits = targetCoords.withinBuildingsLimits(lats, lons);
            profileInd = ~profileWithinBuildingsLimits;
            profileInd([1 end]) = true; % Always include endpoints
            [targetIsVisible, canShowObstruction, obstructingLocation, obstructingElevation] = ...
                computeLOSUsingTerrainProfile(observerSite, observerCoords, criticalAngles(k), lats(profileInd), lons(profileInd), map);
            
            % If terrain profile intersects buildings region, also test LOS
            % against buildings
            profileIntersectsBuildings = any(profileWithinBuildingsLimits);
            if profileIntersectsBuildings
                % Pass in the sceneTri in case multiple sites need to compute
                % los using raytracing. This avoids recomputing sceneTri every
                % iteration.
                [targetIsVisibleBldgs, obstructingLocationBldgs, obstructingElevationBldgs, sceneTri] = ...
                    computeLOSUsingRayTracing(observerCoords, targetCoords, map, sceneTri);
                
                % Use buildings obstruction if it exists and terrain
                % obstruction does not, or if it is closer than terrain
                % obstruction
                if ~targetIsVisibleBldgs && ~targetIsVisible
                    bldgsDistToObstruction = rfprop.internal.MapUtils.greatCircleDistance(observerLat, observerLon, ...
                        obstructingLocationBldgs(1), obstructingLocationBldgs(2));
                    terrainDistToObstruction = rfprop.internal.MapUtils.greatCircleDistance(observerLat, observerLon, ...
                        obstructingLocation(1), obstructingLocation(2));
                    useBuildingsObstruction = (bldgsDistToObstruction <= terrainDistToObstruction);
                else
                    useBuildingsObstruction = ~targetIsVisibleBldgs && targetIsVisible;
                end
                
                if useBuildingsObstruction
                    targetIsVisible = targetIsVisibleBldgs;
                    canShowObstruction = true;
                    obstructingLocation = obstructingLocationBldgs;
                    obstructingElevation = obstructingElevationBldgs;
                end
            end
        end
    end
    
    % Populate output or line data for plot
    if outputRequested
        status(k) = targetIsVisible;
    else
        startIDs{k} = observerSite.UID;
        endIDs{k} = targetSite.UID;
        dEucKmFormat = mat2str(round(dEuc/1000, rfprop.Constants.MaxDistanceInfoDecimalPlaces));
        lineInfoDist = [message('shared_channel:rfprop:LOSDescriptionDistance', dEucKmFormat).getString, '<br>'];
        observerLocation = observerCoords.DisplayLatitudeLongitude;
        targetLocation = targetCoords.DisplayLatitudeLongitude;
        if targetIsVisible || ~canShowObstruction
            visprofile{k} = [targetIsVisible targetIsVisible];
            latitudes{k} = [observerLocation(1) targetLocation(1)];
            longitudes{k} = [observerLocation(2) targetLocation(2)];
            elevations{k} = [observerCoords.AntennaHeightAboveTerrainReference targetAntennaElevation];
            showObstructionPoint{k} = false;
        else % Partially visible line
            visprofile{k} = [true true false];
            latitudes{k} = [observerLocation(1) obstructingLocation(1) targetLocation(1)];
            longitudes{k} = [observerLocation(2) obstructingLocation(2) targetLocation(2)];
            elevations{k} = [observerCoords.AntennaHeightAboveTerrainReference obstructingElevation targetAntennaElevation];
            showObstructionPoint{k} = useBuildingsObstruction;
        end
        
        % Set description
        if targetIsVisible
            lineInfo{k} = [lineInfoDist, message('shared_channel:rfprop:LOSDescriptionVisibleStatus').getString];
        else
            lineInfo{k} = [lineInfoDist, message('shared_channel:rfprop:LOSDescriptionNotVisibleStatus').getString];
        end
    end
end

if outputRequested
    varargout = {status};
else
    % Show sites (but no animation and do not launch window)
    if isViewerInitiallyVisible && ~map.Visible
        return % Abort if Site Viewer has been closed (test before show)
    end
    show(observerSite,'AntennaSiteCoordinates',observerCoords,'Map',map, ...
        'Animation','none','EnableWindowLaunch',false);
    if isViewerInitiallyVisible && ~map.Visible
        return % Abort if Site Viewer has been closed (test before show)
    end
    show(targetSites,'AntennaSiteCoordinates',targetsCoords,'Map',map, ...
        'Animation','none','EnableWindowLaunch',false);
    
    % Get viewer and plot visibility line (with blocking until complete)    
    numLOS = numel(longitudes);
    losLocations = cell(1,numLOS);
    obstructionPointLocations = cell(1,numLOS);
    lineNVs = cell(numLOS);
    pointNVs = cell(numLOS);
    graphicsToRemove = {};
    for i = 1:numLOS
        numLines = numel(longitudes{i}) - 1;
        lineLocations = zeros(numLines + 1, 3);
        currentVisProfile = visprofile{i};
        losID = map.getId(1);
        losID = ['los' num2str(losID{1})];
        targetID = targetSites(i).UID;
        observerID = observerSite.UID;
        % Remove existing LOS lines between the sites
        oldLOSLines = map.getGraphic(observerID, 'los');
        if (isfield(oldLOSLines, targetID))
            graphicsToRemove = [graphicsToRemove; oldLOSLines.(targetID)];
        end
        map.disassociateSiteGraphics(targetID, observerID, 'los');
        
        % Remove existing rays between the sites
        oldRayLines = map.getGraphic(observerID, 'rays');
        if (isfield(oldRayLines, targetID))
            map.removeFromLegendEntities(oldRayLines.(targetID));
            graphicsToRemove = [graphicsToRemove; oldRayLines.(targetID)];
        end
        map.disassociateSiteGraphics(targetID, observerID, 'rays');
        
        % Initialize an empty point name-value pairs, to be filled if there
        % is an obstruction.
        pointNV = {};
        % There can only ever be 1 or 2 lines per LOS total (2 or 3 points total per LOS)
        if (numLines == 1)
            lineLocations(1,:) = [latitudes{i}(1), longitudes{i}(1), elevations{i}(1)];
            lineLocations(2,:) = [latitudes{i}(2), longitudes{i}(2), elevations{i}(2)];
            if (all(currentVisProfile(:)))
                color = visibleColor;
            else
                color = obstructedColor;
            end
            lineNV = {"Color", color, "ShowArrow", true, "Width", 9, "Description", lineInfo{i},...
                "Name", message('shared_channel:rfprop:LOSName').getString,"ID", {losID}};
            map.associateSiteGraphics(targetID, observerID, 'los', losID);
        else
            % Only a double-line can be obstructed.
            lineLocations(1,:) = [latitudes{i}(1), longitudes{i}(1), elevations{i}(1)];
            lineLocations(2,:) = [latitudes{i}(2), longitudes{i}(2), elevations{i}(2)];
            lineLocations(3,:) = [latitudes{i}(3), longitudes{i}(3), elevations{i}(3)];
            colors = {visibleColor, obstructedColor};
            
            % The line without the arrow must have a different width
            fullLOSID = {losID; ['los' losID]};
            lineNV = {"Color", colors, "ShowArrow", [false, true], "Width", [2, 9], "Description", lineInfo{i}, ...
                "Name", message('shared_channel:rfprop:LOSName').getString,"ID", fullLOSID};
            map.associateSiteGraphics(targetID, observerID, 'los', fullLOSID);
            
            if showObstructionPoint{i}
                % Get obstruction point location and elevation
                obstructionPointLocation = lineLocations(2,:);
                obstructionPointLat = obstructionPointLocation(1);
                obstructionPointLon = obstructionPointLocation(2);
                obstructionPointElevation = obstructionPointLocation(3);
                
                % Add obstruction point to list for plotting
                obstructionPointLocations{i} = obstructionPointLocation;
                
                % Obstruction elevation is with reference to TerrainSource reference,
                % but we need geoid reference for description
                if map.UseTerrain && strcmpi(map.TerrainSource.HeightReference,'ellipsoid')
                    obstructionPointElevation = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(...
                        obstructionPointElevation,obstructionPointLat,obstructionPointLon);
                end
                
                % Add a point to where the red and green line touch
                pointDesc = [num2str(obstructionPointLat) ',' num2str(obstructionPointLon) '<br>' ...
                    message('shared_channel:rfprop:LOSObstructionElevationDescription', round(obstructionPointElevation)).getString];
                pointID = ['point' losID];
                pointNV = {"Color", obstructedColor, "Description", {pointDesc}, ...
                    "Name", message('shared_channel:rfprop:LOSObstructionName').getString, ...
                    "ID", {pointID}};
                % Add the point to both site's graphics
                map.associateSiteGraphics(targetID, observerID, 'los', pointID);
            end
            
        end
        losLocations{i} = lineLocations;
        lineNV = [lineNV {"LinkedGraphics", {{targetID; observerID}}}]; %#ok<AGROW>
        lineNVs{i} = lineNV;
        pointNVs{i} = pointNV;
    end
    if isViewerInitiallyVisible && ~map.Visible
        return % Abort if Site Viewer has been closed (test before line plot)
    end
    map.remove(graphicsToRemove);
    losviewer = map.Instance.GlobeViewer.getViewer("LOS");
    if (~isViewerInitiallyVisible && isempty(animation))
        animation = 'zoom';
    end
    losviewer.los(losLocations, obstructionPointLocations, map.Instance.GlobeViewer.getViewer("Line"), ...
        map.Instance.GlobeViewer.getViewer("Point"), lineNVs, pointNVs, [showObstructionPoint{:}], animation);
    % Turn clipping on
    map.turnClippingOn();
    
    if (~map.Visible)
        map.Visible = true;
    end
    map.bringToFront;
end
end

function [isvis, obsLoc, obsEl, tri] = computeLOSUsingRayTracing(observerCoords, targetCoords, map, tri)

% Initialize outputs to correspond to no-obstruction case
obsLoc = [];
obsEl = [];

% Use ray tracing to compute obstruction point
observerPos = observerCoords.enuFromRegionCenter;
targetPos = targetCoords.enuFromRegionCenter;
% If the triangulation wasn't passed in, compute it based on the sceneModel
if isempty(tri)
    gsc = map.GeographicSceneModel;
    tri = Geometry.meshToTriangulation(gsc.sceneMesh(gsc.BuildingsModel.BuildingsLimits(1:2), ...
            gsc.BuildingsModel.BuildingsLimits(3:4), ...
            gsc.BuildingsModel.BuildingsCenter));
end
env = struct( ...
    "Triangulation", tri, ...
    "RayTracer", matlabshared.internal.StaticSceneRayTracer(tri), ...
    "SharpEdgeFlags", comm.internal.geometry.getSharpEdges(tri));    
[intPt,~,isobstructed] = comm.internal.geometry.firstIntersect( ...
    env, observerPos, targetPos, 'segment');
isvis = ~isobstructed;

% Compute obstruction point
if isobstructed
    [obsLat,obsLon,obsEl] = targetCoords.geodeticFromRegionCenter(intPt);
    obsLoc = [obsLat obsLon];
end
end

function [isvis, canShowObstruction, obsLoc, obsEl] = computeLOSUsingTerrainProfile(observerSite, observerCoords, criticalAngle, lats, lons, map)

% Initialize outputs to correspond to no-obstruction case
isvis = true;
canShowObstruction = false;
obsLoc = [];
obsEl = [];

% Return early if not enough data points. Analysis requires at least three
% points: two end points and one middle point.
if numel(lats) <= 2
    return
end

% Compute visibility status by looking for first obstruction along terrain
% between observer and target, which is located where elevation angle from
% observer is greater than angle between observer and target.
rxprofile = rxsite('Name', 'internal.terrainsamplesite', ... % Specify to avoid default site naming
    'AntennaHeight', 0, ...
    'Latitude', lats(2:end-1), ...
    'Longitude', lons(2:end-1));
[~, el] = angle(observerSite, rxprofile, 'euclidean', 'Map', map, ...
    'SourceAntennaSiteCoordinates', observerCoords);
obstructionInd = find(el >= criticalAngle, 1);
isvis = isempty(obstructionInd);

% If first sample is obstruction, assume there is complete obstruction
% (i.e. no point to show as visible)
canShowObstruction = obstructionInd ~= 1;

% Compute obstruction point
if ~isvis
    obstructingSite = rxprofile(obstructionInd);
    obstruction = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(...
        obstructingSite, map);
    obsLoc = obstruction.DisplayLatitudeLongitude;
    obsEl = obstruction.SurfaceHeightAboveTerrainReference;
end
end

function isvis = computeLOSUsingSpheroid(observerCoords, targetCoords)

    % Create spheroid and geodetic coordinates for Earth center
    w = map.geodesy.internal.wgs84Spheroid;
    lat0 = 0;
    lon0 = 0;
    h0 = -w.SemimajorAxis;

    % Get Cartesian coordinates for observer and target from Earth center
    obsLatLon = observerCoords.LatitudeLongitude;
    [Xobs,Yobs,Zobs] = rfprop.internal.MapUtils.geodetic2enu(lat0, lon0, h0, ...
        obsLatLon(1), obsLatLon(2), observerCoords.AntennaHeightAboveEllipsoid, w);
    tgtLatLon = targetCoords.LatitudeLongitude;
    [Xtgt,Ytgt,Ztgt] = rfprop.internal.MapUtils.geodetic2enu(lat0, lon0, h0, ...
        tgtLatLon(1), tgtLatLon(2), targetCoords.AntennaHeightAboveEllipsoid, w);

    % Compute shortest distance from line segment connecting observer and
    % target to Earth center
    dsquared = lineSegmentToOrigin(Xobs,Yobs,Zobs,Xtgt,Ytgt,Ztgt);

    % Visibility exists if there is no intersection between the line-of-sight
    % and the Earth, which occurs if the shortest distance from the line to the
    % Earth center is greater than the Earth radius. Use spherical Earth radius
    % for approximation.
    isvis = dsquared > rfprop.Constants.SphericalEarthRadius^2;
end

function [s,t] = lineSegmentToOrigin(x1,y1,z1,x2,y2,z2)
% Return the square, s, of 3-D Euclidean distance to the origin from the
% line segment connecting point 1 (x1,y1,z1) to point 2 (x2,y2,z2).
% Optionally return parameter t, which indicates the fraction of the
% distance along the segment connecting point 1 to point 2 at which the
% point closest to the origin is located. The value of t is constrained to
% the interval [0, 1]. When t == 0, the point 1 is closer to the origin
% than any other point along the segment. When t == 1, point 2 is closer to
% the origin than any other point along the segment. Otherwise, the closest
% point to the origin falls somewhere along the segment between point 1 and
% point 2.
%
% Using a parametric form for a 3-D line:
%
%     x(t) = x1 .* t + x2 .* (1 - t)
%     y(t) = y1 .* t + y2 .* (1 - t)
%     z(t) = z1 .* t + z2 .* (1 - t)
%
% this function finds the value of t that minimizes:
%
%     s(t) = x(t).^2 + y(t).^2 + z(t).^2
%
% subject to the contraints that 0 <= t <= 1. Using a, b, and c, defined
% below,
%
%     s(t) = a*t.^2 - 2*b*t + c
%
% The first derivative is s'(t) = 2*a*t - 2*b. Setting s'(t) = 0 indicates
% that the unconstrained minimum occurs at t = b ./ a, which can then be
% clamped to [0, 1] using a min-max pattern.
%
% This function is elementwise and fully vectorized: The inputs can have
% arbitrary size, as long as their sizes all match, and the output, s and
% t, will match the inputs in size. (An arbitrary number of point-pairs
% can be processed in a single call).

    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;
    a = dx .* dx + dy .* dy + dz .* dz;
    b = x2 .* dx + y2 .* dy + z2 .* dz;
    c = x2 .* x2 + y2 .* y2 + z2 .* z2;
    t = min(max(b ./ a, 0), 1);
    s = t .* (a .* t - 2 * b) + c;  % a*t.^2 - 2*b*t + c via Horner's rule
end

function res = validateResolution(p)

try
    %Resolution validation
    res = p.Results.Resolution;
    
    if ischar(res) || isstring(res)
        res = validatestring(res, {'auto'}, 'los', 'Resolution');
    else
        validateattributes(res,{'numeric'}, ...
            {'positive','real','finite','nonsparse','scalar'}, ...
            'los','Resolution');
        res = double(res);
    end
catch e
    throwAsCaller(e);
end
end

function  res = autoResolution(dGeo)
res = rfprop.Constants.DefaultLOSTerrainResolution;
maxSamples = rfprop.Constants.MaxLOSNumTerrainSamples;
if (dGeo/res > maxSamples) % Saturate at resolution that generates max samples
    res = dGeo/maxSamples;
elseif res >= dGeo % Guarantee at least one sample
    res = dGeo/2;
end
end