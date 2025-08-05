function varargout = raytrace(txs, rxs, varargin)
%

% Copyright 2019-2024 The MathWorks, Inc.

% Validate sites
validateattributes(txs,{'txsite'},{'nonempty'},'raytrace','',1);
validateattributes(rxs,{'rxsite'},{'nonempty'},'raytrace','',2);
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(txs, rxs);

% Process optional name/value pairs
p = inputParser;
if mod(numel(varargin),2) % Odd number of inputs
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Animation', '');
p.addParameter('EnableWindowLaunch', true);
p.addParameter('Type', 'power');
p.addParameter('ColorLimits', []);
p.addParameter('Colormap', 'jet');
p.addParameter('ShowLegend', true);
p.addParameter('Map', []);
p.parse(varargin{:});

% Get Site Viewer visibility state and validate web graphics
outputRequested = nargout > 0;
if ~outputRequested
    viewer = rfprop.internal.Validators.validateMap(p, 'raytrace', usingCartesian);
    isViewerInitiallyVisible = viewer.Visible;
else
    isViewerInitiallyVisible = false;
end

% Validate and get parameters
if usingCartesian
    [map, sceneMesh, sceneMaterials] = rfprop.internal.Validators.validateCartesianMap(p);
    siteCoordSys = 'cartesian';
else
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'raytrace');
    % The GeographicSceneModel will create the sceneMesh and sceneMaterials
    sceneMaterials = "";
    siteCoordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, siteCoordSys);
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(...
    p, isViewerInitiallyVisible, 'raytrace');
[pm, usingDefaultPropModel] = validateRaytracingPropagationModel(p,txs,rxs);
if isa(pm, 'rfprop.RayTracing')
    rtpm = pm;
else
    for cpm = pm.PropagationModels
        if isa(cpm, 'rfprop.RayTracing')
            rtpm = cpm;
            break
        end
    end
end

numReflections = validateNumReflections(rtpm);

% Verify CoordinateSystem consistency
if ~strcmp(rtpm.CoordinateSystem,txs(1).CoordinateSystem)
    error(message('shared_channel:rfprop:RaytraceCoordinateSystemInconsistent'))
end

% Validate Type and and get corresponding default parameter values
allowedTypes = {'power','pathloss'};
[type, ~, defaultColorLimits] = ...
    rfprop.internal.Validators.validateType(p, allowedTypes, 'raytrace');

% Validate and get graphical parameters
cmap = rfprop.internal.Validators.validateColorMap(p, 'raytrace');
clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'raytrace');
showLegend = rfprop.internal.Validators.validateShowLegend(p, 'raytrace');

% Get site antenna coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, map);
rxsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(rxs, map);

warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
warnCleanup = onCleanup(@()warning(warnState));

if usingCartesian
    if strcmp(map,'none')
        % With no map data, there is always a single line-of-sight path.
        % Get this result using a trivial triangulation map and forcing
        % numReflections to 0.
        txsrxsPos = [txsCoords.AntennaPosition; rxsCoords.AntennaPosition];
        minZ = min(txsrxsPos(:,3));
        zpt = minZ - 10;
        pts = [0 0 zpt;
            1 0 zpt;
            0 1 zpt];
        tri = triangulation(1:3, pts);
        numReflections = 0;
        sceneMesh = Geometry.mesh(tri);
    else
        tri = map;
    end
    
    % Package material information
    mtl = struct('Material',  "", 'MaterialName', "");
    for k = 1:numel(sceneMaterials)
        mtl(k) = struct('Material', string(sceneMaterials(k)), 'MaterialName', string(sceneMaterials(k)));
    end

    srfMaterial = string(rtpm.SurfaceMaterial);
    srfMaterialStr = srfMaterial;

    % If the sceneMesh has no materials defined, then it is an STL file
    % that needs the materials from the propagationModel.
    if sceneMesh.GetMaterials.Dimensions == 1 % the first material is always empty in the geometry library
        AddMaterial(sceneMesh, srfMaterialStr);
        Geometry.setMaterial(sceneMesh, srfMaterialStr);
    end

    if srfMaterial == "custom"
        srfMaterial = [rtpm.SurfaceMaterialPermittivity; rtpm.SurfaceMaterialConductivity];
        % If the propagation-model override is being used, then this custom
        % material will override the entire "mtl" struct in the next block.
        % If 'scene' is being used, then all empty materials
        % (unspecified/unrecognized) will be replaced with the custom
        % material "srfMaterial".
    end
    
    % When the material is 'auto', the materials will default to
    % those found in the scene. Otherwise the materials of the prop model
    % will override those of the scene.
    if ~strcmp(srfMaterialStr, 'auto')
        sceneTri = Geometry.meshToTriangulation(sceneMesh);
        % Convert the triangulation back to a mesh but with a uniform
        % material
        sceneMesh = Geometry.mesh(sceneTri, srfMaterialStr);
        % Overwrite the current instance of the materials table with a
        % single entry for our prop model material.
        mtl = struct('Material', srfMaterial, 'MaterialName', srfMaterialStr);
    else
        % Replace all unspecified and unrecognized materials (denoted with a
        % "") with the default material specified in rfprop.Constants. This
        % will be concrete unless rfprop.Constants is modified.
        for k = 1:numel(mtl)
            if strcmp(mtl(k).Material, "")
                mtl(k).Material = rfprop.Constants.DefaultRaytracingSurfaceMaterial;
                mtl(k).MaterialName = rfprop.Constants.DefaultRaytracingSurfaceMaterial;
            end
        end
    end
    
    % Formulate internal site objects for ray tracing
    cTxs = comm.internal.Site(txs);
    cRxs = comm.internal.Site(rxs);
else
    % Site Viewer constructs its own GeographicSceneModel, but if raytrace
    % is used without a siteviewer then a GeographicSceneModel must be
    % constructed here for mesh-related operations.
    if isa(map, "siteviewer")
        gsc = map.GeographicSceneModel;
        % The only geographic case in which we have "matched" materials
        % right now is when siteviewer has buildings. If we create a
        % dictionary mapping materials to their matched material, the
        % geoscene can easily determine what material each building in the
        % mesh should be assigned. If a key isn't found, it'll simply use
        % the default material found in the "" entry.
        matchedMaterials = map.Materials.MatchedCatalogMaterial;
        matchedMaterials(matchedMaterials == "") = rtpm.BuildingsMaterial;
        materialsDict = dictionary([map.Materials.Material], matchedMaterials);
    else
        gsc = rfprop.internal.GeographicSceneModel;
        if isa(map, "terrain.internal.TerrainSource")
            gsc.Terrain = map.Name;
            gsc.TerrainSource = map;
        end
        % Without a siteviewer, there is no materials matching in the
        % geographic case. So we can simply create an empty dictionary that
        % will be populated by the materials specified by the prop model.
        materialsDict = dictionary("", rtpm.BuildingsMaterial);
    end
    dgc = distance(txs,rxs,'geodesic');
    validatePropagationDistance(dgc);

    % Package material information.

    % Initialize the materials struct array that will be used by the
    % comm.internal.raytrace function
    % For most materials, the Material and MaterialName field will be the
    % same. For custom materials, the MaterialName will contain "custom"
    % and the Material will be the parameters provided in the prop model.
    mtl = struct("Material", "fillerMaterial", "MaterialName", "fillerName");
    mtl = mtl([]); % Initialize the struct as truly empty so it can be concatenated with later
    materialKeys = materialsDict.keys;
    for k = 1:materialsDict.numEntries
        mtl(k) = struct("Material", materialsDict(materialKeys(k)), "MaterialName", materialsDict(materialKeys(k)));
    end

    % Handle materials prior to generating the scene to avoid having to
    % regenerate the scene with new materials.
    bldgsMaterial = string(rtpm.BuildingsMaterial);
    bldgsMaterialStr = bldgsMaterial;
    if bldgsMaterial == "custom"
        % Replace the custom material entry with the custom material being
        % used
        bldgsMaterial = [rtpm.BuildingsMaterialPermittivity; rtpm.BuildingsMaterialConductivity];
        % To differentiate between custom buildings material and custom
        % terrain materials, add the index of the material to the material
        % string.
        bldgsMaterialStr = "custom" + (materialsDict.numEntries+1);
        % Map the custom material to itself so the geoscene can set the
        % materials of the buildings accordingly
        materialsDict(bldgsMaterialStr) = bldgsMaterialStr;
    end
    mtl(end+1) = struct("Material", bldgsMaterial, "MaterialName", bldgsMaterialStr);

    terrainMaterial = string(rtpm.TerrainMaterial);
    if terrainMaterial == "custom"
        terrainMaterial = [rtpm.TerrainMaterialPermittivity; rtpm.TerrainMaterialConductivity];
        % To differentiate between custom buildings material and custom
        % terrain materials, add the index of the material to the material
        % string.
        gsc.setTerrainMaterial("custom" + (materialsDict.numEntries+1));
        % We don't need to add the terrain material to the materials
        % dictionary, because terrain will only have a single material
        % always. The materials dictionary is to allow us to specify
        % buildings materials at runtime.
    else
        gsc.setTerrainMaterial(terrainMaterial);
    end
    mtl(end+1) = struct("Material", terrainMaterial, "MaterialName", gsc.TerrainMaterial);
    % Set the materials of the geoscene based on the prop model. If the
    % prop model does not use the 'auto' materials, then the materials
    % will be overridden with those of the prop model.
    if ~strcmp(bldgsMaterialStr, 'auto')
        % If the buildings materials are to be overridden, then all
        % materials in the buildings should map to the specified material.
        % Override our materials dictionary with a single default value of
        % "" mapping to our specified material.
        materialsDict = dictionary("", bldgsMaterialStr);
    else
        % If the scene data's materials are to be used, then only map empty
        % values to the default value for buildings materials found in
        % rfprop.Constants.
        materialsDict("") = rfprop.Constants.DefaultRaytracingBuildingsMaterial;
    end

    % Generate triangulation of 3D environment for tx/rx locations
    txsLatLon = txsCoords.LatitudeLongitude;
    rxsLatLon = rxsCoords.LatitudeLongitude;
    lats = [txsLatLon(:,1); rxsLatLon(:,1)];
    lons = [txsLatLon(:,2); rxsLatLon(:,2)];
    [latmin, latmax] = bounds(lats(:));
    [lonmin, lonmax] = bounds(lons(:));
    
    % If within building limits, then the boundaries of the scene are
    % simply the buildings limits.
    if all(txsCoords.withinBuildingsLimits(lats, lons))
        sceneMesh = gsc.sceneMesh(gsc.BuildingsModel.BuildingsLimits(1:2), ...
            gsc.BuildingsModel.BuildingsLimits(3:4), ...
            gsc.BuildingsModel.BuildingsCenter, ...
            "TerrainResolution", terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution, ...
            "Materials", materialsDict);
        regionCenter = gsc.BuildingsModel.BuildingsCenter;
    else
        % Since not all sites are within the limits, a buffer region must
        % be constructed around them. Then, we must determine if any of
        % this buffer region intersects the buildings region.
        [latlim, lonlim] = bufferRegion([latmin latmax],[lonmin lonmax]);
        
        % First we need to determine if the region contains the buildings
        terrainOnlyTri = isempty(gsc.BuildingsModel);
        if ~terrainOnlyTri
            bldgsLimits = gsc.BuildingsModel.BuildingsLimits;
            terrainOnlyTri = (bldgsLimits(1) > latlim(2)) || (latlim(1) > bldgsLimits(2)) || ...
                (bldgsLimits(3) > lonlim(2)) || (lonlim(1) > bldgsLimits(4));
        end
        % If only terrain is used, then the region should be centered
        % around itself
        if terrainOnlyTri
            if strcmp(gsc.TerrainSource, 'none')
                terrainResolution = 250;
            else
                terrainResolution = gsc.TerrainSource.IntrinsicResolutionArcLength;
            end
            sceneMesh = gsc.sceneMesh(latlim, lonlim, "TerrainResolution", terrainResolution, "Materials", materialsDict);
            regionCenter = [sum(latlim)/2, sum(lonlim)/2, 0];
        else
            % If buildings are used as well, then the region should be
            % centered at the center of the buildings
            [sceneMesh, regionCenter] = gsc.sceneMesh(latlim, lonlim, gsc.BuildingsModel.BuildingsCenter, ...
                "TerrainResolution", terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution, ...
                "Materials", materialsDict);
        end
    end
    tri = Geometry.meshToTriangulation(sceneMesh);

    % Remove any materials that aren't used by the sceneMesh and reorder them
    % such that they match the indices used by the mesh. We can't rely on
    % the ordering specified by the buildings in the mesh because we don't
    % know what order the materials may have appeared in our buildings
    % list, and this order determines their index in the geometry library.
    [~, sceneMaterials] = Geometry.extractMesh(sceneMesh);
    uniqueMaterialIndices = unique(sceneMaterials(:,4));
    uniqueMaterials = struct("Material", "fillerMaterial", "MaterialName", "fillerName");
    allMaterialNames = [mtl.MaterialName];
    for k = 1:numel(uniqueMaterialIndices)
        currentMaterial = sceneMesh.GetMaterial(uniqueMaterialIndices(k)).GetName;
        currentMaterialIndex = strcmpi(currentMaterial, allMaterialNames);
        % In cases where the terrain material matches the buildings
        % material, there will be multiple entries matching. Use only the
        % first of these matches.
        materialEntry = mtl(currentMaterialIndex);
        uniqueMaterials(k) = materialEntry(1);
    end
    mtl = uniqueMaterials;
    
    % Compute txs coordinate info
    txsCoords.RegionCenter = regionCenter;
    rxsCoords.RegionCenter = regionCenter;
    txsLatLon = txsCoords.LatitudeLongitude;
    txsAntennaElevation = txsCoords.AntennaHeightAboveTerrainReference;
    txsPos = txsCoords.enuFromRegionCenter';
    txsAxes = antennaAxes(txs,txsLatLon,txsAntennaElevation,regionCenter);
    
    % Compute rxs coordinate info
    rxsLatLon = rxsCoords.LatitudeLongitude;
    rxsAntennaElevation = rxsCoords.AntennaHeightAboveTerrainReference;
    rxsPos = rxsCoords.enuFromRegionCenter';
    rxsAxes = antennaAxes(rxs,rxsLatLon,rxsAntennaElevation,regionCenter);
    
    % Formulate internal tx site objects for ray tracing
    numTx = numel(txs);
    cTxs = repmat(comm.internal.Site, 1, numTx);
    txPos = mat2cell(txsPos, 3, ones(1, numTx));
    ant = {txs.Antenna};
    freq = {txs.TransmitterFrequency};
    [cTxs.Position] = txPos{:};
    [cTxs.Antenna] = ant{:};
    [cTxs.OrientationAxes] = txsAxes{:};
    [cTxs.Frequency] = freq{:};
    
    % Formulate internal rx site objects for ray tracing
    numRx = numel(rxs);
    cRxs = repmat(comm.internal.Site, 1, numRx);
    rxPos = mat2cell(rxsPos, 3, ones(1, numRx));
    ant = {rxs.Antenna};
    [cRxs.Position] = rxPos{:};
    [cRxs.Antenna] = ant{:};
    [cRxs.OrientationAxes] = rxsAxes{:};
end

% Perform ray tracing analysis for all txs to all rxs
env = struct( ...
    "Triangulation", tri, ...
    "RayTracer", matlabshared.internal.StaticSceneRayTracer(tri), ...
    "SharpEdgeFlags", comm.internal.geometry.getSharpEdges(tri));

isImageMethod = strcmp(rtpm.Method, 'image');
if isImageMethod
    cfg = struct( ...
        'Method', 'image', ...
        'NumReflections', numReflections, ...
        'MaxAbsolutePathLoss', rtpm.MaxAbsolutePathLoss, ...
        'MaxRelativePathLoss', rtpm.MaxRelativePathLoss);
else
    cfg = struct( ...
        'Method', 'sbr', ...
        'AngularSeparation',  rtpm.AngularSeparation, ...
        'MaxNumReflections',  rtpm.MaxNumReflections, ...
        'MaxNumDiffractions', rtpm.MaxNumDiffractions, ...
        'MaxAbsolutePathLoss', rtpm.MaxAbsolutePathLoss, ...
        'MaxRelativePathLoss', rtpm.MaxRelativePathLoss, ...
        'UseGPU', rtpm.UseGPU);
end

raytraceArgs = {env, mtl, cTxs, cRxs, cfg, sceneMesh};

if isa(map,'siteviewer') && ~isempty(map.ProgressDialogHandle) && isvalid(map.ProgressDialogHandle)
    raytraceArgs{end+1} = map.ProgressDialogHandle;
end
rays = comm.internal.raytrace(raytraceArgs{:});

% Update rays to have geographic location
if ~outputRequested
    raysToPlot = {};
    txsPlot = txsite.empty;
    rxsPlot = rxsite.empty;
end

% Path loss thresholds
PathLossThresholds = [rtpm.MaxAbsolutePathLoss rtpm.MaxRelativePathLoss];
usePathLossThresholds = any(PathLossThresholds);

% Process rays
if isImageMethod
    tempRayInt = repmat(struct('Type', 'Reflection', 'Location', [10; 10; 0], ...
        'MaterialName', strings(0)), 1, rtpm.MaxNumReflections + rtpm.MaxNumDiffractions);
end
updateRays = isImageMethod || ~usingCartesian;
isCompositePropagationModel = isa(pm,'rfprop.CompositePropagationModel');
for txInd = 1:numel(txs)
    tx = txs(txInd);
    if ~usingCartesian
        txLatLon = txsLatLon(txInd,:);
        txElevation = txsAntennaElevation(txInd);
    end
    
    for rxInd = 1:numel(rxs)
        rx = rxs(rxInd);
        if ~usingCartesian
           rxLatLon = rxsLatLon(rxInd,:);
           rxElevation = rxsAntennaElevation(rxInd);
        end
        
        % Get rays, where each ray object corresponds to a propagation
        % path between the tx and rx (whether LOS or reflected)
        txrxRays = rays{txInd,rxInd};
        if isempty(txrxRays)
            continue
        end
        
        % Update rays with specific path losses (e.g. rain, gas)
        if isCompositePropagationModel
            for rayInd = 1:numel(txrxRays)
                ray = txrxRays(rayInd);
                txrxRays(rayInd).PathLoss = ray.PathLoss + ...
                    pm.specificPathLoss(rx, tx, ...
                    ray.PropagationDistance, ray.AngleOfDeparture);
            end
        end

        % Filter rays using path loss thresholds
        if usePathLossThresholds
            PL = [txrxRays.PathLoss];
            keepRays = PL <= min(PathLossThresholds(1), min(PL)+PathLossThresholds(2));
            txrxRays = txrxRays(keepRays);
        end
        
        % Update rays for geographic locations and interactions
        if updateRays
            for rayInd = 1:numel(txrxRays)
                ray = txrxRays(rayInd);

                % Update ray with geographic locations
                if ~usingCartesian
                    [rayLats,rayLons,rayElevs] = rayGeographicLocations(...
                        ray, txsCoords, txLatLon, txElevation, rxLatLon, rxElevation);
                    ray.CoordinateSystem = 'Geographic';
                    ray.TransmitterLocation = [rayLats(1); rayLons(1); rayElevs(1)];
                    ray.ReceiverLocation = [rayLats(end); rayLons(end); rayElevs(end)];
                end

                % Update interactions
                if ~ray.LineOfSight
                    numRayInt = length(ray.Interactions);
                    if isImageMethod
                        rayInt = tempRayInt(1:numRayInt);
                        for intactIdx = 1:numRayInt
                            rayInt(intactIdx).Type = ray.Interactions(intactIdx).Type;
                            if ~usingCartesian
                                rayInt(intactIdx).Location = [rayLats(intactIdx+1); 
                                    rayLons(intactIdx+1); rayElevs(intactIdx+1)];
                            else
                                rayInt(intactIdx).Location = ray.Interactions(intactIdx).Location;
                            end
                            rayInt(intactIdx).MaterialName = ray.Interactions(intactIdx).MaterialName;
                        end
                        ray.Interactions = rayInt;
                    elseif ~usingCartesian % SBR
                        for intactIdx = 1:numRayInt
                            ray.Interactions(intactIdx).Location = [rayLats(intactIdx+1); 
                                rayLons(intactIdx+1); rayElevs(intactIdx+1)];
                        end
                    end
                end

                txrxRays(rayInd) = ray;
            end
        end

        % Arrays for display
        if ~outputRequested
            for rayInd = 1:numel(txrxRays)
                raysToPlot{end + 1} = txrxRays(rayInd); %#ok<*AGROW>
                % Use builtin to concatenate to avoid expensive validation
                txsPlot = builtin('horzcat',txsPlot,tx);
                rxsPlot = builtin('horzcat',rxsPlot,rx);
            end
        end
        
        % Update rays cell array
        rays{txInd,rxInd} = txrxRays;
    end
end

% Set ray UIDs so they can be plotted later
for k = 1:numel(rays) % ray cells
    for ii = 1:numel(rays{k}) % individual rays
        rays{k}(ii).UID = rays{k}(ii).createUIDNumber;
    end
end

% Return rays if output requested
if outputRequested
    varargout = {rays};
    return
end

% Show sites
if isViewerInitiallyVisible && ~viewer.Visible
    return % Abort if Site Viewer has been closed
end
show(txs,'AntennaSiteCoordinates',txsCoords,'Map',viewer, 'Animation','none','EnableWindowLaunch',false);
show(rxs,'AntennaSiteCoordinates',rxsCoords,'Map',viewer, 'Animation','none','EnableWindowLaunch',false);

% Warn if no paths to plot and return early
allrays = [raysToPlot{:}];
if isempty(allrays)
    if usingDefaultPropModel
        warning(message('shared_channel:rfprop:NoPathsFoundUsingDefaultPropagationModel'));
    else
        warning(message('shared_channel:rfprop:NoPathsFoundUsingCustomPropagationModel'));
    end
    return
end

% Plot rays
plot(allrays,'Type', type, ...
    'TransmitterSite', txsPlot, ...
    'ReceiverSite', rxsPlot, ...
    'AssociateWithSites', true, ...
    'ColorLimits', clim, ...
    'ColorMap', cmap, ...
    'ShowLegend', showLegend, ...
    'Map', viewer, ...
    'Animation', animation, ...
    'EnableWindowLaunch', enableWindowLaunch);
end

function [pm, usingDefaultPropModel] = validateRaytracingPropagationModel(p,txs,rxs)

try
    if ismember('PropagationModel', p.UsingDefaults)
        pm = propagationModel('raytracing',...
            'CoordinateSystem',char(txs(1).CoordinateSystem));
        usingDefaultPropModel = true;
    else
        pm = p.Results.PropagationModel;
        if ischar(pm) || isstring(pm)
            % Normalize name and validate it corresponds to ray tracing
            % model. Allow hyphenated and spaced forms of "raytracing".
            pmName = pm;
            raytracematches = ["ray tracing" "ray-tracing"];
            if startsWith(pmName,"ray") && any(startsWith(raytracematches, pmName))
                pmName = 'raytracing';
            else
                pmName = validatestring(pmName, {'raytracing'}, 'raytrace', 'PropagationModel');
            end
            pm = propagationModel(pmName);
        else
            % Validate model supports multi-path
            if ~isa(pm,'rfprop.PropagationModel') || ~isscalar(pm) || ~pm.isMultipathModel
                error(message('shared_channel:rfprop:RaytracePropagationModelValue'))
            end
        end
        usingDefaultPropModel = false;
    end
    
    % Validate antenna site properties
    pm.validateTransmitterSites(txs);
    pm.validateReceiverSites(rxs);
catch e
    throwAsCaller(e)
end
end

function validatePropagationDistance(dgc)
%validatePropagationDistance   Validate propagation great circle distance

try
    % Verify that geodesic, or great circle, distance from tx to rx does
    % not exceed maximum propagation distance. Re-use limit for terrain.
    dgc = dgc(:);
    maxPropDistance = rfprop.Constants.MaxPropagationDistanceUsingTerrain;
    if any(isnan(dgc)) || any(dgc > maxPropDistance)
        error(message('shared_channel:rfprop:RaytraceDistanceGreaterThanMax',round(maxPropDistance/1000)));
    end
catch e
    throwAsCaller(e);
end
end

function numReflections = validateNumReflections(pm)

try
    % Compute NumReflections from MaxNumReflections
    rfprop.internal.Validators.validateMaxNumReflections(pm,"raytrace");
    numReflections = 0:pm.MaxNumReflections;
catch e
    throwAsCaller(e);
end
end

function [rayLats,rayLons,rayElevs] = rayGeographicLocations(...
    ray, txsCoords, txLatLon, txElevation, rxLatLon, rxElevation)
%rayGeographicLocations   Compute geographic locations of ray vertices

% Initialize location variables, which includes the tx and rx locations
% along with any reflections
numLocations = 2 + ray.NumInteractions;
rayLats = zeros(numLocations,1);
rayLons = zeros(numLocations,1);
rayElevs = zeros(numLocations,1);

% Ray path starts with tx location
rayLats(1) = txLatLon(1);
rayLons(1) = txLatLon(2);
rayElevs(1) = txElevation;

% Add reflection locations
tolerance = 1e-3;
if ~ray.LineOfSight
    for intactIdx = 1:ray.NumInteractions
        rayIntactPos = ray.Interactions(intactIdx).Location;
        [rayLats(intactIdx+1),rayLons(intactIdx+1),rayElev] = ...
            txsCoords.geodeticFromRegionCenter(rayIntactPos');
        rayElev(abs(rayElev) < tolerance) = 0; % Snap very small values to 0
        rayElevs(intactIdx+1) = rayElev;
    end
end

% Ray path ends with rx location
rayLats(end) = rxLatLon(1);
rayLons(end) = rxLatLon(2);
rayElevs(end) = rxElevation;
end

function axes = antennaAxes(sites, latlon, h, regionCenter)
%antennaAxes   Return antenna axes for geographic antenna sites

axes = cell(1, numel(sites)); 
lat0 = regionCenter(1);
lon0 = regionCenter(2);
h0 = regionCenter(3);

rotMtx = rfprop.internal.MapUtils.enuRotation( ...
    lat0, lon0, h0, latlon(:,1), latlon(:,2), h);

for siteInd = 1:numel(sites)
    site = sites(siteInd);
    antAngle = site.AntennaAngle;
    
    % Populate antenna axes for site
    phi = antAngle(1);
    if numel(antAngle) > 1
        theta = antAngle(2);
    else
        theta = 0;
    end
    localAntennaAxes = ([1 0 0; 0 cosd(theta) sind(theta); 0 -sind(theta) cosd(theta)] * ...
        [cosd(phi) sind(phi) 0; -sind(phi) cosd(phi) 0; 0 0 1])';

    axes{siteInd} = rotMtx{siteInd} * localAntennaAxes;
end
end

function [latlim, lonlim] = bufferRegion(latlim, lonlim)

% Generate four corner locations from limits
lats = [latlim latlim];
lons = [lonlim(1) lonlim(1) lonlim(2) lonlim(2)];

% Get locations defined by buffer distance from corners
buffer = terrain.internal.TerrainSource.TerrainRegionBuffer;
[latfwd, lonfwd] = rfprop.internal.MapUtils.greatCircleForward(...
    lats, lons, buffer, [0 90 180 270]);

% Get new corners using buffer distance
[latmin, latmax] = bounds(latfwd(:));
[lonmin, lonmax] = bounds(lonfwd(:));
latlim = [latmin latmax];
lonlim = [lonmin, lonmax];
end