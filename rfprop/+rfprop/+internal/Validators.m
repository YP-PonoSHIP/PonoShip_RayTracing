classdef (Sealed, Hidden) Validators
    %Validators   Reserved for MathWorks internal use only
    
    %   Copyright 2018-2024 The MathWorks, Inc. 
    
    properties(Constant)        
        ScreenSize = getScreenSize
    end
    
    methods(Static)
        function [type, levels, defaultColorLimits, defaultLegendTitle] = validateType(p, allowedTypes, fcnName, data)
            
            % Validate Type
            try
                % Default value must be defined by caller
                type = p.Results.Type;
                type = validatestring(type, allowedTypes, fcnName, 'Type');
            catch e
                throwAsCaller(e);
            end
            
            % Get default values according to plot type
            isCoverage = ismember(fcnName,{'coverage'});
            switch(type)
                case 'power'
                    if isCoverage
                        % Coverage map thresholds at a single value by default
                        defaultLevels = -100; % Unit: dBm
                    end
                    defaultColorLimits = [-120 -5];
                    defaultLegendTitle = message('shared_channel:rfprop:PowerLegendTitle').getString;
                case 'efield'
                    if isCoverage
                        % Coverage map thresholds at a single value by default
                        defaultLevels = 40; % Unit: dBuV/m
                    end
                    defaultColorLimits = [20 135];
                    defaultLegendTitle = message('shared_channel:rfprop:EfieldLegendTitle').getString;
                case 'sinr'
                    defaultColorLimits = [-5 20];
                    defaultLegendTitle = message('shared_channel:rfprop:SINRLegendTitle').getString;
                case 'pathloss'
                    defaultColorLimits = [45 160];
                    defaultLegendTitle = message('shared_channel:rfprop:PathlossLegendTitle').getString;
                case 'custom'
                    data = data(isfinite(data)); % Remove Inf, -Inf
                    [minData, maxData] = bounds(data);
                    defaultColorLimits = [minData maxData];
                    defaultLegendTitle = '';
            end
            
            % Compute levels from color limits
            if ~isCoverage
                defaultLevels = linspace(defaultColorLimits(1), defaultColorLimits(2), 64);
            end
            
            % If "custom" data is allowed, validate "Levels" parameter 
            levels = defaultLevels;
            if ismember('custom',allowedTypes)
                try
                    if ~ismember('Levels', p.UsingDefaults)
                        levels = p.Results.Levels;
                        validateattributes(levels, {'numeric'}, ...
                            {'real','nonnan','nonsparse','vector','nonempty'}, fcnName, 'Levels');
                        if ~iscolumn(levels)
                            levels = levels(:);
                        end
                        levels = double(levels);
                    end
                catch e
                    throwAsCaller(e);
                end
            end
        end
        
        function levels = validateLevels(p, defaultLevels, fcnName)
            try
                if ismember('Levels', p.UsingDefaults)
                    levels = defaultLevels;
                else
                    levels = p.Results.Levels;
                    validateattributes(levels, {'numeric'}, ...
                        {'real','nonnan','nonsparse','vector','nonempty'}, fcnName, 'Levels');
                    if ~iscolumn(levels)
                        levels = levels(:);
                    end
                    levels = double(levels);
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function [res, isAutoRes] = validateResolution(p, maxrange, fcnName)
            
            try
                res = p.Results.Resolution;
                isAutoRes = ischar(res) || isstring(res);
                if isAutoRes
                    validatestring(res, {'auto'}, fcnName, 'Resolution');
                    
                    % Compute resolution to scale with maxrange. The smaller the resolution,
                    % the higher the quality but at the cost of performance. Select value
                    % that will produce a grid that is approximately 100-by-100, which
                    % provides good quality without too much performance cost. Note that
                    % when there are multiple tx, maxrange is a vector. Choose the largest
                    % value since it corresponds to largest coverage map that will be shown. 
                    % If all tx have the same maxrange, then quality of the coverage map is 
                    % the same whether a single or multiple tx is specified.
                    d = 2*max(maxrange);
                    gridLength = 99;
                    res = d / gridLength;
                else
                    validateattributes(res,{'numeric'}, ...
                        {'positive','real','finite','nonnan','nonsparse','scalar','<',min(maxrange)}, ...
                        fcnName,'Resolution');
                    res = double(res);
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function datarange = validateDataRange(txslatlon, maxrange, res, useTerrain)
            
            try
                % Define a range for data collection that goes a little beyond maxrange.
                % This is necessary so that interpolation can be performed at the maxrange
                % boundary to provide a smoother outer contour (see contourmap). In
                % addition, this provides a little padding so the image fits within the
                % screen on zoom.
                datarange = maxrange + 2*res;
                
                % Saturate data range at max prop distance. This may have already been done
                % in rfprop.PropagationModel/range but is needed here to account for buffer
                % added beyond maxrange. This cannot only be done here because range can be
                % called as its own method.
                if useTerrain
                    maxDataRange = rfprop.Constants.MaxPropagationDistanceUsingTerrain;
                    maxDataRangeErrorID = 'shared_channel:rfprop:PathlossDistanceGreaterThanTerrainMax';
                else
                    maxDataRange = rfprop.Constants.MaxPropagationDistance;
                    maxDataRangeErrorID = 'shared_channel:rfprop:PathlossDistanceGreaterThanMax';
                end
                datarange(datarange > maxDataRange) = maxDataRange;
                
                % Validate that distance between txsites does not exceed propagation limit,
                % since all txsites are included in signal strength computations of region
                % around all other txsites
                numTxs = size(txslatlon,1);
                if numTxs > 1
                    lats = txslatlon(:,1);
                    lons = txslatlon(:,2);
                    for txInd = 1:numTxs
                        gc = rfprop.internal.MapUtils.greatCircleDistance(lats(txInd), lons(txInd), lats, lons);
                        gcmax = max(gc(:));
                        if isnan(gcmax) || gcmax > maxDataRange
                            error(message(maxDataRangeErrorID, round(maxDataRange/1000)));
                        end
                    end
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function maxrange = defaultMaxRange(txslatlon, pm, map)
            
            % Initialize output
            if isMultipathModel(pm)
                defaultMaxRange = rfprop.Constants.DefaultMultipathModelMaxRange;
            else
                defaultMaxRange = rfprop.Constants.DefaultMaxRange;
            end
            
            % If using terrain propagation model and Site Viewer has
            % buildings, calculate maxrange to extend to furthest building
            numTx = size(txslatlon,1);
            if requiresTerrain(pm) && isa(map,'siteviewer') && map.HasBuildings
                % Get lats/lons for buildings region corners
                bldgLimits = map.BuildingsLimits;
                bldgLimitsLats = [bldgLimits(1); bldgLimits(1); bldgLimits(2); bldgLimits(2)];
                bldgLimitsLons = [bldgLimits(3); bldgLimits(4); bldgLimits(3); bldgLimits(4)];
                
                % Calculate distance from sites to building limits
                maxrange = zeros(1,numTx);
                for k = 1:numTx
                    distToBldgLimits = rfprop.internal.MapUtils.greatCircleDistance(...
                        txslatlon(k,1), txslatlon(k,2), bldgLimitsLats, bldgLimitsLons);
                    maxrange(k) = min(defaultMaxRange,max(distToBldgLimits));
                end
            else
                maxrange = repmat(defaultMaxRange, 1, numTx);
            end
        end
        
        function maxrange = validateNumericMaxRange(maxrange, pm, numSites, map, fcnName)
            
            try
                maxPropDistance = rfprop.Constants.MaxPropagationDistance;
                validateattributes(maxrange,{'numeric'}, ...
                    {'positive','real','finite','nonnan','nonsparse','scalar', ...
                    '<=',maxPropDistance},fcnName,'MaxRange');
                maxrange = double(maxrange);

                % Validate MaxRange is greater than terrain resolution
                if requiresTerrain(pm)
                    minMaxRange = terrainProfileResolution(pm, map);
                    if maxrange < minMaxRange
                        error(message('shared_channel:rfprop:MinMaxRangeTerrain',num2str(minMaxRange)));
                    end
                end
                
                % Enforce limit if using terrain
                if isa(map,'siteviewer')
                    useTerrain = map.UseTerrain;
                else
                    useTerrain = ~strcmp(map,'none');
                end
                if useTerrain
                    maxPropDistance = rfprop.Constants.MaxPropagationDistanceUsingTerrain;
                    if maxrange > maxPropDistance
                        error(message('shared_channel:rfprop:MaxRangeTerrain',round(maxPropDistance/1000)));
                    end
                end
                
                % Repeat range for each tx
                maxrange = repmat(maxrange, 1, numSites);
            catch e
                throwAsCaller(e);
            end
        end
        
        function rxAntennaHeight = validateReceiverAntennaHeight(p, fcnName)
            
            try
                rxAntennaHeight = p.Results.ReceiverAntennaHeight;
                validateattributes(rxAntennaHeight, {'numeric'}, ...
                    {'nonnegative','real','finite','nonnan','nonsparse','scalar'}, ...
                    fcnName, 'ReceiverAntennaHeight');
            catch e
                throwAsCaller(e);
            end
        end
        
        function colors = validateColors(p, fcnName)
            
            try
                colors = p.Results.Colors;
                if ~ismember('Colors', p.UsingDefaults)
                    % Validate that neither Colormap nor ColorLimits specified
                    usingDefaultColormap = ismember('Colormap', p.UsingDefaults);
                    usingDefaultColorLimits = ismember('ColorLimits', p.UsingDefaults);
                    if ~usingDefaultColormap || ~usingDefaultColorLimits
                        error(message('shared_channel:rfprop:InvalidColorParameters', ...
                            'Colors', 'Colormap', 'ColorLimits'));
                    end
                    
                    validColors = rfprop.internal.ColorUtils.colors;
                    if ischar(colors)
                        colors = validatestring(colors, validColors, fcnName, 'Colors');
                        colors = rfprop.internal.ColorUtils.str2rgb(colors);
                    elseif iscellstr(colors)
                        rgbColors = [];
                        for k = 1:numel(colors)
                            color = validatestring(colors{k}, validColors, fcnName, 'Colors');
                            rgbColors = [rgbColors; rfprop.internal.ColorUtils.str2rgb(color)]; %#ok<*AGROW>
                        end
                        colors = rgbColors;
                    elseif isstring(colors)
                        rgbColors = [];
                        for k = 1:numel(colors)
                            color = validatestring(colors(k), validColors, fcnName, 'Colors');
                            rgbColors = [rgbColors; rfprop.internal.ColorUtils.str2rgb(color)];
                        end
                        colors = rgbColors;
                    else
                        validateattributes(colors, {'numeric'}, ...
                            {'real','finite','nonnan','nonsparse','ncols',3,'>=',0,'<=',1}, ...
                            fcnName, 'Colors');
                    end
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function clim = validateColorLimits(p, defaultColorLimits, fcnName)
            
            try
                clim = p.Results.ColorLimits;
                if ismember('ColorLimits', p.UsingDefaults)
                    clim = defaultColorLimits;
                else
                    validateattributes(clim,{'numeric'}, ...
                        {'real','finite','nonnan','nonsparse','increasing','row','ncols',2}, ...
                        fcnName, 'ColorLimits');
                    clim = double(clim);
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function showLegend = validateShowLegend(p, fcnName)
            
            try
                showLegend = p.Results.ShowLegend;
                validateattributes(showLegend, {'logical'}, {'nonsparse','scalar'}, ...
                    fcnName, 'ShowLegend');
            catch e
                throwAsCaller(e);
            end
        end
        
        function legendTitle = validateLegendTitle(p, defaultLegendTitle, fcnName)
            
            try
                if ismember('LegendTitle', p.UsingDefaults)
                    legendTitle = defaultLegendTitle;
                else
                    legendTitle = p.Results.LegendTitle;
                    validateattributes(legendTitle, {'char','string'}, {'scalartext'}, ...
                        fcnName, 'LegendTitle');
                    legendTitle = char(legendTitle);
                    
                    % Replace newlines with <br>
                    legendTitle = strrep(legendTitle,newline,'<br>');
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function cmap = validateColorMap(p, fcnName)
            
            try
                cmap = p.Results.Colormap;
                if ischar(cmap) || isstring(cmap)
                    validateattributes(cmap,{'char','string'}, {'scalartext'}, ...
                        fcnName, 'Colormap');
                    
                    % Use colormap approach of performing feval on text
                    cmap = char(lower(cmap));
                    k = min(strfind(cmap,'('));
                    if ~isempty(k)
                        cmap = feval(cmap(1:k-1),str2double(cmap(k+1:end-1)));
                    else
                        cmap = feval(cmap);
                    end
                else
                    % If string, let colormap function validate
                    if ~ischar(cmap) && ~isstring(cmap)
                        validateattributes(cmap,{'numeric'}, ...
                            {'real','finite','nonnan','nonsparse','ncols',3,'>=',0,'<=',1}, ...
                            fcnName, 'Colormap');
                    end
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function transparency = validateTransparency(p, fcnName)
            
            try
                transparency = p.Results.Transparency;
                validateattributes(transparency, {'numeric'}, ...
                    {'real','finite','nonnan','nonsparse','scalar','>=',0,'<=',1}, ...
                    fcnName, 'Transparency');
                transparency = double(transparency);
            catch e
                throwAsCaller(e)
            end
        end
        
        function rxGain = validateReceiverGain(p, fcnName)
            
            try
                rxGain = p.Results.ReceiverGain;
                if ismember('ReceiverGain', p.UsingDefaults)
                    expShape = 'vector'; % Allow vector if computing gain for multiple tx given rx
                else
                    expShape = 'scalar';
                end
                validateattributes(rxGain,{'numeric'}, ...
                    {'real','finite','nonnan','nonsparse',expShape},fcnName,'ReceiverGain');
            catch e
                throwAsCaller(e);
            end
        end
        
        function imageSize = validateImageSize(dataSize, maxSizeFactor, maxImageSize, res, visualName)
            
            % Compute size factor limited by max image size and error if
            % data size exceeds max
            maxImageSizeFactor = floor(maxImageSize./dataSize);
            if maxImageSizeFactor < 1
                resStr = mat2str(round(res,rfprop.Constants.MaxDistanceInfoDecimalPlaces));
                error(message('shared_channel:rfprop:ContourmapTooMuchData', visualName, resStr))
            end
            
            % Compute size factor limited by screen size. If data size
            % exceeds screen size, factor is 1 and therefore no
            % interpolation is performed to generate image.
            screenSize = rfprop.internal.Validators.ScreenSize;
            screenSizeFactor = ceil(screenSize./dataSize);
            
            % Compute size factor
            sizeFactor = min([maxSizeFactor, screenSizeFactor, maxImageSizeFactor]);
            
            % Compute image size so that resolution (length/numSteps) is
            % inversely proportional to size factor
            numSteps = dataSize-1;
            imageSize = numSteps*sizeFactor+1;
        end
        
        function [map, mesh, materials, mapStruct] = validateCartesianMap(p)
            
            try
                % Suppress possible Triangulation warning if input
                % triangulation or STL file has un-used vertices
                warnState = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
                warnCleanup = onCleanup(@()warning(warnState));

                % Check if input is a pre-validated struct
                mapInput = p.Results.Map;
                isDefaultSource = ismember('Map', p.UsingDefaults);
                if ~isDefaultSource && isstruct(mapInput)
                    if isequal(sort(string(fieldnames(mapInput))), ...
                            ["Map"; "Materials"; "Mesh"])
                        map = mapInput.Map;
                        mesh = mapInput.Mesh;
                        materials = mapInput.Materials;
                        mapStruct = mapInput;
                        
                        return
                    else
                        error(message('shared_channel:rfprop:CartesianSiteInvalidMap'))
                    end
                end

                % Parse Map input or current SiteViewer (if any)
                materials = "";
                mesh = Geometry.mesh;
                if isDefaultSource
                    % Check current SiteViewer
                    if ~isempty(siteviewer.all) && siteviewer.current.IsCartesian
                        if ~isempty(siteviewer.current.ModelMesh)
                            map = Geometry.meshToTriangulation(siteviewer.current.ModelMesh);
                            mesh = siteviewer.current.ModelMesh;
                        else
                            map = 'none';
                        end
                        materials = siteviewer.current.Materials.MatchedCatalogMaterial;
                    else
                        map = 'none';
                    end
                else
                    map = mapInput;
                    istri = isa(map,'triangulation');
                    isviewer = isa(map, 'siteviewer');

                    % Parse if siteviewer
                    if isviewer
                        map = rfprop.internal.Validators.validateMap(p, '', map.IsCartesian);

                        % Must be cartesian
                        if ~map.IsCartesian
                            error(message('shared_channel:rfprop:CartesianSiteInvalidMap'))
                        end

                        % Parse model mesh
                        if isempty(map.ModelMesh)
                            map = 'none';
                        else
                            map = Geometry.meshToTriangulation(map.ModelMesh);
                            istri = isa(map,'triangulation');
                            mesh = siteviewer.current.ModelMesh;
                        end
                        materials = siteviewer.current.Materials.MatchedCatalogMaterial;
                    elseif ~istri && ~matlab.internal.datatypes.isScalarText(map)
                        error(message('shared_channel:rfprop:CartesianSiteInvalidMap'))
                    end
                    
                    % Parse triangulation or STL/glTF file
                    if istri
                        mesh = Geometry.mesh(map);
                    elseif ~strcmp(map,'none')
                        % Validate file name has STL or glb/glTF extension
                        [~,~,ext] = fileparts(map);
                        isSTL = strcmpi(ext, '.stl');
                        isGLTF = any(strcmpi(ext, {'.glb', '.gltf'}));
                        if ~(isSTL || isGLTF)
                            error(message("shared_channel:rfprop:CartesianSiteInvalidMap"));
                        end
                        
                        % Rely on stlread and the geometry library for
                        % remaining validation. Note that relative names
                        % and names on the path are accepted.
                        if isSTL
                            map = stlread(map);
                            mesh = Geometry.mesh(map);
                        else % glTF
                            mesh = Geometry.meshRead(map);
                            map = Geometry.meshToTriangulation(mesh);

                            % Materials from file
                            numMaterials = mesh.GetMaterials.Dimensions - 1;
                            materials = arrayfun(@(x) mesh.GetMaterial(x).GetName, ...
                                1:numMaterials);

                            % Match the materials against the catalog
                            materials = rfprop.internal.matchMaterial(materials);
                        end
                    end
                end
            catch e
                throwAsCaller(e);
            end

            % Validated map struct
            if nargout >= 4
                mapStruct = struct("Map", map, "Materials", materials, ...
                    "Mesh", mesh);
            end
        end
        
        function [map, mapStruct] = validateMapTerrainSource(p, fcnName)
            try
                % Check if input is a pre-validated struct
                mapInput = p.Results.Map;
                isDefaultSource = ismember('Map', p.UsingDefaults);
                if ~isDefaultSource && isstruct(mapInput)
                    if isequal(fieldnames(mapInput), {'Map'})
                        map = mapInput.Map;
                        mapStruct = mapInput;
                        
                        return
                    else
                        error(message('shared_channel:rfprop:GeographicSiteInvalidMap'))
                    end
                end

                % Parse Map input or siteviewer (if any)
                if isDefaultSource
                    % Get the current Site Viewer if any exist. Otherwise,
                    % just get the default terrain source. Only use the
                    % current viewer if it is geographic.
                    if ~isempty(siteviewer.all) && strcmp('geographic', ...
                            siteviewer.current.CoordinateSystem)
                        map = siteviewer.current;
                    else
                        map = siteviewer.defaultTerrainName;
                        if ~strcmp(map, 'none')
                            map = terrain.internal.TerrainSource.createFromSettings(map);
                        end
                    end
                else
                    if isa(mapInput,'siteviewer')
                        map = rfprop.internal.Validators.validateMap(p, fcnName);
                    elseif isa(mapInput,'terrain.internal.TerrainSource') || ...
                            strcmp(mapInput,'none')
                        map = mapInput;
                    else
                        terrainChoices = terrain.internal.TerrainSource.terrainchoices;
                        terrainName = validatestring(mapInput, terrainChoices, fcnName, 'Map');
                        map = terrain.internal.TerrainSource.createFromSettings(terrainName);
                    end
                end
            catch e
                throwAsCaller(e);
            end

            % Validated map struct
            if nargout >= 2
                mapStruct = struct("Map", map);
            end
        end

        function validateMapCoordinateSystemMatch(map, coordinateSystem)
            if strcmp(map, 'none')
                return
            end

            if (isa(map, 'siteviewer') && ~strcmp(map.CoordinateSystem, coordinateSystem))
                % Siteviewer type does not match desired coordinate system
                % Throw error mismatch
                if strcmp(map.CoordinateSystem, 'geographic')
                    error(message('shared_channel:rfprop:CartesianSiteInvalidMap'));
                else
                    error(message('shared_channel:rfprop:GeographicSiteInvalidMap'));
                end
            end

            if strcmp(coordinateSystem, 'geographic') && isa(map, 'triangulation')
                % Triangulation incompatible with geographic scenes
                % Throw error mismatch
                error(message('shared_channel:rfprop:GeographicSiteInvalidMap'));
            elseif strcmp(coordinateSystem, 'cartesian') && ischar(map)
                % Terrain name incompatible with cartesian scenes
                % Throw error mismatch
                error(message('shared_channel:rfprop:CartesianSiteInvalidMap'));
            end
        end
        
        function terrainSource = validateTerrainSource(mapOrTerrain, fcnName)
            
            if nargin < 2
                fcnName = '';
            end
            
            try
                if isa(mapOrTerrain,'siteviewer')
                    terrainSource = mapOrTerrain.Visualizer.Controller.TerrainSource;
                elseif isa(mapOrTerrain,'terrain.internal.TerrainSource') || strcmp(mapOrTerrain,'none')
                    terrainSource = mapOrTerrain; 
                else
                    terrainChoices = terrain.internal.TerrainSource.terrainchoices;
                    terrainName = validatestring(mapOrTerrain, terrainChoices, fcnName, 'Map');
                    terrainSource = terrain.internal.TerrainSource.createFromSettings(terrainName);
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function usingCartesian = validateCoordinateSystem(varargin)
            
            try
                % Validate that all sites have the same CoordinateSystem value
                expCoordinateSystem = '';
                usingCartesian = false;
                for sitesInd = 1:numel(varargin)
                    sites = varargin{sitesInd};
                    if isempty(sites)
                        % Empty is possible with concatenation (e.g. [txsite txsite.empty])
                        continue
                    elseif isempty(expCoordinateSystem)
                        usingCartesian = strcmp(sites(1).CoordinateSystem,'cartesian');
                        if usingCartesian
                            expCoordinateSystem = 'cartesian';
                        else
                            expCoordinateSystem = 'geographic';
                        end
                    end
                    
                    siteCoords = cell(numel(sites),1);
                    for k = 1:numel(sites)
                        siteCoords{k} = sites(k).CoordinateSystem;
                    end
                    if ~all(strcmp(siteCoords, expCoordinateSystem))
                        error(message('shared_channel:rfprop:CoordinateSystemInconsistent'));
                    end
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function validateGeographicSites(sites, fcnName)
            
            try
                if any(strcmp([sites.CoordinateSystem],'cartesian'))
                    error(message('shared_channel:rfprop:CoordinateSystemCartesianNotSupported', fcnName));
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function viewer = validateMap(p, fcnName, usingCartesian)
            if nargin < 3
                usingCartesian = false;
            end
            try
                isDefaultMap = ismember('Map', p.UsingDefaults);
                if isDefaultMap
                    % Get current Site Viewer
                    if usingCartesian
                        viewer = siteviewer.current('cartesian');
                    else
                        viewer = siteviewer.current('geographic');
                    end
                else
                    map = p.Results.Map;
                    validateattributes(map, {'siteviewer'}, {'scalar'}, ...
                        fcnName, 'Map');
                    viewer = map;
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function coords = validateAntennaSiteCoordinates(coords, sites, viewer, fcnName)
            
            try
                if isempty(coords)
                    coords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(sites, viewer);
                else
                    validateattributes(coords, {'rfprop.internal.AntennaSiteCoordinates'}, ...
                        {'scalar'}, fcnName, 'AntennaSiteCoordinates');
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function pm = validatePropagationModel(p, fcnName)
            %validatePropagationModel   Validate propagation model
            
            pm = p.Results.PropagationModel;
            if ischar(pm) || isstring(pm)
                pmName = rfprop.PropagationModel.validateName(...
                    pm, fcnName, 'PropagationModel');
                pm = propagationModel(pmName);
            else
                validateattributes(pm, {'rfprop.PropagationModel'}, ...
                    {'scalar'}, fcnName, 'PropagationModel');
            end
        end
        
        function pm = validateGeographicPropagationModel(p, map, fcnName)
            %validateGeographicPropagationModel   Validate propagation model for geographic sites
            
            try
                terrainSource = rfprop.internal.Validators.validateTerrainSource(map, fcnName);
                if ismember('PropagationModel', p.UsingDefaults)
                    % Choose default propagation model based on whether terrain
                    % is enabled (Longley-Rice) or disabled (free space)
                    if strcmp(terrainSource,'none')
                        pm = rfprop.FreeSpace;
                    else
                        pm = rfprop.LongleyRice;
                    end
                else
                    pm = rfprop.internal.Validators.validatePropagationModel(p, fcnName);
                end
                % Early validation for TIREM that library can be accessed
                if isa(pm,'rfprop.TIREM')
                    tirem.internal.Validators.validateSetup
                end
                rfprop.internal.Validators.validateTerrainAvailability(terrainSource, pm);
            catch e
                throwAsCaller(e)
            end
        end
        
        function pm = validateCartesianPropagationModel(p, map, fcnName)
            %validateCartesianPropagationModel   Validate propagation model for Cartesian sites
            
            try
                if ismember('PropagationModel', p.UsingDefaults)
                    % Choose default propagation model based on whether map
                    % is 'none' or a 3-D environment (STL or triangulation)
                    if strcmp(map,'none')
                        pm = rfprop.FreeSpace;
                    else
                        pm = rfprop.RayTracing;
                    end
                else
                    pm = rfprop.internal.Validators.validatePropagationModel(p, fcnName);
                end
                
                % Disallow terrain propagation models
                if requiresTerrain(pm)
                    error(message('shared_channel:rfprop:CartesianSiteInvalidPropagationModel'));
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function pm = validateMaxNumReflections(pm, fcnName)
            %validateMaxNumReflections   Validate MaxNumReflections for ray tracing model
            
            try
                % Validate MaxNumReflections limit for ray tracing model
                if isa(pm,'rfprop.RayTracing')
                    if strcmp(pm.Method, 'image')
                        numReflectionsLimit = 2;
                    else % sbr
                        numReflectionsLimit = 100;
                    end
                    validateattributes(pm.MaxNumReflections, {'numeric'}, ...
                        {'scalar','integer','<=',numReflectionsLimit}, fcnName, 'MaxNumReflections');
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function pm = validateMaxNumDiffractions(pm, fcnName)
            %validateMaxNumDiffractions   Validate MaxNumDiffractions for
            %SBR ray tracing model in point-to-area functions. 
            
            try
                % Validate MaxNumDiffractions limit for point-to-area functions
                if isa(pm,'rfprop.RayTracing') && strcmp(pm.Method, 'sbr')
                    validateattributes(pm.MaxNumDiffractions, {'numeric'}, ...
                        {'scalar','integer','<=',1}, fcnName, 'MaxNumDiffractions');
                end
            catch e
                throwAsCaller(e)
            end
        end
                
        function [animation, enableWindowLaunch] = validateGraphicsControls(p, visible, fcnName)
            %validateGraphicsControls   Validate graphic control parameters
            
            try
                % Validate Animation. Note that empty is valid default to defer
                % determination of which animation to use.
                animation = p.Results.Animation;
                if ~isempty(animation)
                    animation = validatestring(animation, ...
                        {'none','fly','zoom'}, fcnName, 'Animation');
                else
                    % If the viewer is not visible, then default to zoom
                    % animation
                    if (~visible)
                        animation = 'zoom';
                    end
                end
                
                % Validate EnableWindowLaunch
                if isfield(p.Results,'EnableWindowLaunch')
                    enableWindowLaunch = p.Results.EnableWindowLaunch;
                    validateattributes(enableWindowLaunch, {'logical'}, {'nonsparse','scalar'}, ...
                        fcnName, 'EnableWindowLaunch');
                else
                    enableWindowLaunch = true;
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function color = validateLineColor(p, colorName, fcnName)
            
            try
                color = p.Results.(colorName);
                if ischar(color) || isstring(color)
                    % Get list of valid colors.
                    validColors = rfprop.internal.ColorUtils.colors;
                    color = validatestring(color, validColors, fcnName, colorName);
                    color = rfprop.internal.ColorUtils.str2rgb(color); % Convert to RGB
                else
                    validateattributes(color, {'numeric'}, ...
                        {'real','finite','nonnan','nonsparse','row','ncols',3,'>=',0,'<=',1}, ...
                        fcnName, colorName);
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function validateTerrainAvailability(mapOrTerrainSource, pm)
            %validateTerrainAvailability   Validate terrain availability
            
            try
                if (nargin < 2) || requiresTerrain(pm)
                    % Get terrain source
                    terrainSource = rfprop.internal.Validators.validateTerrainSource(mapOrTerrainSource);
                    
                    % Error if no terrain source
                    if strcmp(terrainSource,'none')
                        error(message('shared_channel:rfprop:TerrainRequired'));
                    end
                    
                    % Error if terrain is unavailable
                    if ~terrainSource.isLocationAvailable
                        if terrainSource.IsURLLocation
                            error(message('shared_terrain:terrain:TerrainNoInternet', ...
                                terrainSource.Name));
                        else
                            error(message("shared_terrain:terrain:TerrainFolderNotFound", ...
                                terrainSource.Name, terrainSource.Location))
                        end
                    end
                end
            catch e
                throwAsCaller(e);
            end
        end
        
        function rxLocationsLayout = validateReceiverLocationsLayout(p, pm, txslatlon, fcnName)
            
            try
                rxLocationsLayout = p.Results.ReceiverLocationsLayout;
                
                if ismember('ReceiverLocationsLayout', p.UsingDefaults)
                    % Use radial for terrain propagation models and grid
                    % for non-terrain models
                    if pm.requiresTerrain
                        rxLocationsLayout = 'radial';
                    else
                        rxLocationsLayout = 'grid';
                    end
                end
                
                % Force grid layout if any txsite is near pole, since
                % radial uses great-circle forward algorithm which breaks
                % down near poles
                txslat = txslatlon(:,1);
                if any(txslat > 89.9) || any(txslat < -89.9)
                    rxLocationsLayout = 'grid';
                end
                
                rxLocationsLayout = validatestring(rxLocationsLayout, ...
                    {'grid','radial'}, fcnName, 'ReceiverLocationsLayout');
            catch e
                throwAsCaller(e);
            end
        end
        
        function validateModel(model)
            if ~isempty(model) % Allow default empty value
                validateattributes(model.Points, {'numeric'}, {'ncols',3}, '', 'Model');
            end
        end
    end
end

function screenSize = getScreenSize

% Get screen position and put in MATLAB size format
screenPosition = get(groot, 'ScreenSize');
screenSize = screenPosition([4 3]);
end
