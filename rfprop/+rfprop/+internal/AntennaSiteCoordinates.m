classdef (Sealed, Hidden) AntennaSiteCoordinates < handle
    %AntennaSiteCoordinates   Reserved for MathWorks internal use only
    
    %   Copyright 2019-2021 The MathWorks, Inc.
    
    properties(SetAccess = private)
        CoordinateSystem
        AntennaPosition
        LatitudeLongitude
        AntennaHeight
        CustomData
    end
    
    properties
        RegionCenter = [0 0 0]
    end
    
    properties(Dependent, SetAccess = private)
        DisplayLatitudeLongitude
        GroundHeightAboveGeoid
        GroundHeightAboveEllipsoid
        GroundHeightAboveTerrainReference
        SurfaceHeightAboveGeoid
        SurfaceHeightAboveEllipsoid
        SurfaceHeightAboveTerrainReference
        AntennaHeightAboveGeoid
        AntennaHeightAboveEllipsoid
        AntennaHeightAboveTerrainReference
    end
    
    properties(Access = private)
        TerrainSource
        BuildingsLimits
        BuildingsTerrainTriangulation
        pGroundHeightAboveGeoid
        pGroundHeightAboveEllipsoid
        pSurfaceHeightAboveGeoid
        pSurfaceHeightAboveEllipsoid
        pDisplayLatitudeLongitude
    end
    
    methods(Static)
        function coords = createFromAntennaSites(sites, map)
            %createFromAntennaSites   Create coordinates object from antenna sites
            
            numSites = numel(sites);
            if strcmp(sites(1).CoordinateSystem, 'geographic')
                anthts = zeros(numSites,1);
                lat = zeros(numSites,1);
                lon = zeros(numSites,1);
                
                for k = 1:numSites
                    site = sites(k);
                    lat(k) = site.Latitude;
                    lon(k) = site.Longitude;
                    anthts(k) = site.AntennaHeight;
                end
                
                % Convert to double before concatenation so no precision is
                % lost if Lat or Lon are different numeric types. Wrap
                % Longitude to [-180,180] range to make subsequent operations
                % safe.
                latlon = [double(lat) wrapTo180(double(lon))];
                
                % Create object to contain site coordinates
                coords = rfprop.internal.AntennaSiteCoordinates(latlon, anthts, map);
            else
                pos = zeros(numSites,3);
                
                for k = 1:numSites
                    pos(k,:) = sites(k).AntennaPosition';
                end
                
                % Create object to contain site coordinates
                coords = rfprop.internal.AntennaSiteCoordinates(pos);
            end
        end
        
        function Z = queryGroundHeightAboveTerrainReference(lats, lons, map)
            %queryGroundHeightAboveTerrainReference   Query ground height above terrain reference at locations
            
            coords = rfprop.internal.AntennaSiteCoordinates([lats, lons], 0, map);
            Z = coords.GroundHeightAboveTerrainReference;
        end
        
        function Z = querySurfaceHeightAboveTerrainReference(lats, lons, map)
            %querySurfaceHeightAboveTerrainReference   Query surface height above terrain reference at locations
            
            coords = rfprop.internal.AntennaSiteCoordinates([lats, lons], 0, map);
            Z = coords.SurfaceHeightAboveTerrainReference;
        end
        
        function Z = querySurfaceHeightAboveGeoid(lats, lons, map)
            %querySurfaceHeightAboveGeoid   Query surface height above geoid at locations
            
            coords = rfprop.internal.AntennaSiteCoordinates([lats, lons], 0, map);
            Z = coords.SurfaceHeightAboveGeoid;
        end
    end
    
    methods
        function coords = AntennaSiteCoordinates(latlonOrPos, antht, map)
            
            if nargin < 2
                coords.CoordinateSystem = 'cartesian';
                coords.AntennaPosition = latlonOrPos;
                return
            else
                coords.CoordinateSystem = 'geographic';
                latlon = latlonOrPos;
            end
            
            % If third column provided in latlon, it corresponds to surface
            % height above geoid
            if size(latlon,2) == 3
                coords.pSurfaceHeightAboveGeoid = latlon(:,3);
                latlon = latlon(:,[1 2]);
            end
            
            % Populate geographic coordinate properties
            coords.LatitudeLongitude = latlon;
            if isscalar(antht)
                antht = antht*ones(size(latlon,1),1);
            end
            coords.AntennaHeight = antht;
            
            % Get terrain source and Site Viewer from input
            viewer = [];
            if isa(map,'siteviewer')
                viewer = map;
                terrainSource = viewer.TerrainSource;
            else
                terrainSource = map;
            end
            
            % Set terrain source for later usage
            if ~strcmpi(terrainSource,'none')
                coords.TerrainSource = terrainSource;
            end
            
            % Set region and buildings data
            if ~isempty(viewer) && viewer.HasBuildings
                coords.BuildingsTerrainTriangulation = viewer.BuildingsTerrainTriangulation;
                coords.BuildingsLimits = viewer.BuildingsLimits;
                coords.RegionCenter = viewer.BuildingsCenter;
            end
        end
        
        function inLimits = withinBuildingsLimits(viewer, lats, lons)
            %withinBuildingsLimits   Return whether locations are within buildings limits
            
            buildingsLimits = viewer.BuildingsLimits;
            if isempty(buildingsLimits)
                inLimits = false(numel(lats),1);
            else
                % Do not include locations that are exactly on the limits
                latlim = buildingsLimits(1:2);
                lonlim = buildingsLimits(3:4);
                inLimits = lats > latlim(1) & lats < latlim(2) & ...
                    lons > lonlim(1) & lons < lonlim(2);
            end
        end
        
        function varargout = enuFromRegionCenter(coords, lats, lons, ht)
            %enuFromRegionCenter   Return Cartesian position at locations
            
            % Use geographic location if specified or else use this
            % object's location
            if nargin == 1
                latlon = coords.LatitudeLongitude;
                lats = latlon(:,1);
                lons = latlon(:,2);
                ht = coords.AntennaHeightAboveTerrainReference;
            else
                lats = lats(:);
                lons = lons(:);
                ht = ht(:);
            end
            
            % Get region center location to use as reference origin
            center = coords.RegionCenter;
            lat0 = center(1);
            lon0 = center(2);
            h0 = center(3);
            
            % Convert geographic location to ENU
            [x,y,z] = rfprop.internal.MapUtils.geodetic2enu(...
                lat0, lon0, h0, lats, lons, ht);
            
            % Return position as single variable or three outputs
            if nargout == 1
                pos = [x y z];
                varargout = {pos};
            else
                varargout = {x,y,z};
            end
        end
        
        function [lat,lon,ht] = geodeticFromRegionCenter(coords, pos)
            %geodeticFromRegionCenter   Return geographic locations at positions
            
            % Get region center location to use as reference origin
            center = coords.RegionCenter;
            lat0 = center(1);
            lon0 = center(2);
            h0 = center(3);
            
            % Convert ENU to geographic location
            [lat,lon,ht] = rfprop.internal.MapUtils.enu2geodetic(...
                lat0, lon0, h0, pos');
        end
        
        function coordsout = extract(coords, ind)
            %extract   Create coordinates object for site index
            
            % Create output coordinates with values corresponding to ind
            if strcmp(coords.CoordinateSystem,'cartesian')
                coordsout = rfprop.internal.AntennaSiteCoordinates(...
                    coords.AntennaPosition(ind,:));
            else
                coordsout = rfprop.internal.AntennaSiteCoordinates(...
                    coords.LatitudeLongitude(ind,:), ...
                    coords.AntennaHeight(ind), ...
                    coords.TerrainSource);
                
                % Set private properties
                coordsout.BuildingsTerrainTriangulation = coords.BuildingsTerrainTriangulation;
                coordsout.BuildingsLimits = coords.BuildingsLimits;
                coordsout.RegionCenter = coords.RegionCenter;
                if ~isempty(coords.pSurfaceHeightAboveGeoid)
                    coordsout.pSurfaceHeightAboveGeoid = coords.pSurfaceHeightAboveGeoid(ind);
                end
                if ~isempty(coords.pSurfaceHeightAboveEllipsoid)
                    coordsout.pSurfaceHeightAboveEllipsoid = coords.pSurfaceHeightAboveEllipsoid(ind);
                end
                if ~isempty(coords.pGroundHeightAboveGeoid)
                    coordsout.pGroundHeightAboveGeoid = coords.pGroundHeightAboveGeoid(ind);
                end
                if ~isempty(coords.pGroundHeightAboveEllipsoid)
                    coordsout.pGroundHeightAboveEllipsoid = coords.pGroundHeightAboveEllipsoid(ind);
                end
                displatlon = coords.DisplayLatitudeLongitude;
                coordsout.pDisplayLatitudeLongitude = displatlon(ind,:);
            end
            
            % Set custom data
            customData = coords.CustomData;
            if ~isempty(customData)
                extractedCustomData = struct();
                fields = fieldnames(customData);
                for fieldInd = 1:numel(fields)
                    field = fields{fieldInd};
                    value = customData.(field);
                    extractedCustomData.(field) = value(ind,:);
                end
                coordsout.CustomData = extractedCustomData;
            end
        end
        
        function addCustomData(coords, varargin)
            %addCustomData   Add custom data to object
            
            if isempty(coords.CustomData)
                customData = struct();
            else
                customData = coords.CustomData;
            end
            
            if strcmp(coords.CoordinateSystem,'cartesian')
                expDataLength = size(coords.AntennaPosition,1);
            else
                expDataLength = size(coords.LatitudeLongitude,1);
            end
            for argInd = 1:2:numel(varargin)
                field = varargin{argInd};
                value = varargin{argInd+1};
                
                % Verify that custom data is same length as number of sites
                validateattributes(value,{'numeric'},{'nrows',expDataLength});
                
                % Add custom data field
                customData.(field) = value;
            end
            
            coords.CustomData = customData;
        end
        
        function latlon = get.DisplayLatitudeLongitude(coords)
            if isempty(coords.pDisplayLatitudeLongitude)
                latlon = coords.computeDisplayLatitudeLongitude;
                coords.pDisplayLatitudeLongitude = latlon;
            else
                latlon = coords.pDisplayLatitudeLongitude;
            end
        end
        
        function Z = get.GroundHeightAboveGeoid(coords)
            if isempty(coords.pGroundHeightAboveGeoid)
                if isempty(coords.pGroundHeightAboveEllipsoid) || isempty(coords.TerrainSource)
                    Z = coords.computeGroundHeight('geoid');
                else
                    latlon = coords.LatitudeLongitude;
                    lats = latlon(:,1);
                    lons = latlon(:,2);
                    Zellipsoid = coords.pGroundHeightAboveEllipsoid;
                    Z = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(Zellipsoid, lats, lons);
                    Z = terrain.internal.TerrainSource.snapMeanSeaLevel(Z);
                end
                coords.pGroundHeightAboveGeoid = Z;
            else
                Z = coords.pGroundHeightAboveGeoid;
            end
        end
        
        function Z = get.GroundHeightAboveEllipsoid(coords)
            if isempty(coords.pGroundHeightAboveEllipsoid)
                if isempty(coords.pGroundHeightAboveGeoid) || isempty(coords.TerrainSource)
                    Z = coords.computeGroundHeight('ellipsoid');
                else
                    latlon = coords.LatitudeLongitude;
                    lats = latlon(:,1);
                    lons = latlon(:,2);
                    Zgeoid = coords.pGroundHeightAboveGeoid;
                    Z = terrain.internal.HeightTransformation.orthometricToEllipsoidal(Zgeoid, lats, lons);
                end
                coords.pGroundHeightAboveEllipsoid = Z;
            else
                Z = coords.pGroundHeightAboveEllipsoid;
            end
        end
        
        function Z = get.GroundHeightAboveTerrainReference(coords)
            terrainSource = coords.TerrainSource;
            if isempty(terrainSource) || strcmpi(terrainSource.HeightReference,'ellipsoid')
                Z = coords.GroundHeightAboveEllipsoid;
            else
                Z = coords.GroundHeightAboveGeoid;
            end
        end
        
        function Z = get.SurfaceHeightAboveGeoid(coords)
            if isempty(coords.pSurfaceHeightAboveGeoid)
                if isempty(coords.pSurfaceHeightAboveEllipsoid) || isempty(coords.TerrainSource)
                    Z = coords.computeSurfaceHeight('geoid');
                else
                    latlon = coords.LatitudeLongitude;
                    lats = latlon(:,1);
                    lons = latlon(:,2);
                    Zellipsoid = coords.pSurfaceHeightAboveEllipsoid;
                    Z = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(Zellipsoid, lats, lons);
                    Z = terrain.internal.TerrainSource.snapMeanSeaLevel(Z);
                end
                coords.pSurfaceHeightAboveGeoid = Z;
            else
                Z = coords.pSurfaceHeightAboveGeoid;
            end
        end
        
        function Z = get.SurfaceHeightAboveEllipsoid(coords)
            if isempty(coords.pSurfaceHeightAboveEllipsoid)
                if isempty(coords.pSurfaceHeightAboveGeoid) || isempty(coords.TerrainSource)
                    Z = coords.computeSurfaceHeight('ellipsoid');
                else
                    latlon = coords.LatitudeLongitude;
                    lats = latlon(:,1);
                    lons = latlon(:,2);
                    Zgeoid = coords.pSurfaceHeightAboveGeoid;
                    Z = terrain.internal.HeightTransformation.orthometricToEllipsoidal(Zgeoid, lats, lons);
                end
                coords.pSurfaceHeightAboveEllipsoid = Z;
            else
                Z = coords.pSurfaceHeightAboveEllipsoid;
            end
        end
        
        function Z = get.SurfaceHeightAboveTerrainReference(coords)
            terrainSource = coords.TerrainSource;
            if isempty(terrainSource) || strcmpi(terrainSource.HeightReference,'ellipsoid')
                Z = coords.SurfaceHeightAboveEllipsoid;
            else
                Z = coords.SurfaceHeightAboveGeoid;
            end
        end
        
        function Z = get.AntennaHeightAboveGeoid(coords)
            Z = coords.SurfaceHeightAboveGeoid + coords.AntennaHeight;
        end
         
        function Z = get.AntennaHeightAboveEllipsoid(coords)
            Z = coords.SurfaceHeightAboveEllipsoid + coords.AntennaHeight;
        end
        
        function Z = get.AntennaHeightAboveTerrainReference(coords)
            Z = coords.SurfaceHeightAboveTerrainReference + coords.AntennaHeight;
        end
        
        function [latout, lonout, htout] = geodetic(coords, pos)
            %geodetic   Return geodetic coordinates
            
            h0 = coords.AntennaHeightAboveTerrainReference;
            latlon = coords.LatitudeLongitude;
            lat0 = latlon(:,1);
            lon0 = latlon(:,2);
            if nargin > 1
                % Return coordinates of input, which is defined in ENU
                % coordinates from this object (which must be scalar)
                validateattributes(coords,{'rfprop.internal.AntennaSiteCoordinates'},{'scalar'});
                [latout, lonout, htout] = ...
                    rfprop.internal.MapUtils.enu2geodetic(lat0, lon0, h0, pos);
            else
                % Return coordinates of this object
                latout = lat0;
                lonout = lon0;
                htout = h0;
            end
        end
        
        function [X,Y,Z] = position(coords, referenceCoords)
            %position   Return Cartesian position coordinates from reference
            
            if strcmp(coords.CoordinateSystem,'cartesian')
                [X,Y,Z] = positionCartesian(coords, referenceCoords);
            else
                [X,Y,Z] = positionGeographic(coords, referenceCoords);
            end
        end
    end
    
    methods(Hidden)
        function [X,Y,Z] = positionCartesian(coords, referenceCoords)
            %positionCartesian   Return position from Cartesian reference
            
            srcpos = referenceCoords.AntennaPosition;
            tgtpos = coords.AntennaPosition;
            
            % Initialize outputs
            numSourceSites = size(srcpos,1);
            X = zeros(size(tgtpos,1), numSourceSites);
            Y = X;
            Z = X;
            
            for srcInd = 1:numSourceSites
                relPos = tgtpos - srcpos(srcInd,:);
                X(:,srcInd) = relPos(:,1);
                Y(:,srcInd) = relPos(:,2);
                Z(:,srcInd) = relPos(:,3);
            end
        end
        
        function [X,Y,Z] = positionGeographic(coords, referenceCoords)
            %positionGeographic   Return ENU position from geographic reference
            
            % Get geodetic antenna coordinates
            [srclat, srclon, srcht] = referenceCoords.geodetic;
            [tgtlat, tgtlon, tgtht] = coords.geodetic;
            
            % Initialize outputs
            numSourceSites = numel(srclat);
            X = zeros(numel(tgtlat), numSourceSites);
            Y = X;
            Z = X;
            
            % Populate ENU coordinates from each source site to all target sites
            for srcInd = 1:numSourceSites
                [X(:,srcInd),Y(:,srcInd),Z(:,srcInd)] = rfprop.internal.MapUtils.geodetic2enu(...
                    srclat(srcInd), srclon(srcInd), srcht(srcInd), tgtlat, tgtlon, tgtht);
            end
            
            % Snap very small values to 0 to increase robustness of sensitive
            % calculations like angle
            tolerance = 1e-6; % One micron
            X(abs(X)<tolerance) = 0;
            Y(abs(Y)<tolerance) = 0;
            Z(abs(Z)<tolerance) = 0;
        end
        
        function Z = computeTriangulationElevation(coords, tri, lats, lons, outputHeightReference)
            %computeTriangulationElevation   Return elevation at triangulation surface
            
            % Get ENU points corresponding to locations
            [x,y,z] = coords.enuFromRegionCenter(lats, lons, 0);
            
            % Compute normal vectors at locations
            [xp,yp,zp] = coords.enuFromRegionCenter(lats, lons, 1);
            xnorm = xp - x;
            ynorm = yp - y;
            znorm = zp - z;
            
            % Query triangulation surface points
            pos = comm.internal.geometry.profileSurface(tri, ...
                [x y z], [xnorm ynorm znorm]);
            
            % Compute elevation of surface point
            [~,~,Z] = coords.geodeticFromRegionCenter(pos);
            
            % The elevation is with respect to the reference specified by
            % the terrain source. Transform it if different from the
            % requested output height reference.
            terrainSource = coords.TerrainSource;
            if ~isempty(terrainSource)
                if strcmpi(terrainSource.HeightReference,'ellipsoid') && strcmpi(outputHeightReference,'geoid')
                    Z = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(Z, lats, lons);
                    Z = terrain.internal.TerrainSource.snapMeanSeaLevel(Z);
                elseif strcmpi(terrainSource.HeightReference,'geoid') && strcmpi(outputHeightReference,'ellipsoid')
                    Z = terrain.internal.HeightTransformation.orthometricToEllipsoidal(Z, lats, lons);
                end
            end
        end
        
        function Z = computeSurfaceHeight(coords, outputHeightReference)
            %computeSurfaceHeight   Return surface height at locations
            
            latlon = coords.LatitudeLongitude;
            lats = latlon(:,1);
            lons = latlon(:,2);
            
            % Use ground elevation if there are no buildings
            bldgsLimits = coords.BuildingsLimits;
            isGeoidReference = strcmpi(outputHeightReference,'geoid');
            if isempty(bldgsLimits)
                if isGeoidReference
                    Z = coords.GroundHeightAboveGeoid;
                else
                    Z = coords.GroundHeightAboveEllipsoid;
                end
                return
            end
            
            % Initialize output
            Z = zeros(size(latlon,1),1);
            
            % Populate surface elevation of locations within buildings region
            inBldgsLimits = coords.withinBuildingsLimits(lats, lons);
            if any(inBldgsLimits)
                % Get locations
                inBldgsLats = lats(inBldgsLimits);
                inBldgsLons = lons(inBldgsLimits);
                
                % Compute elevation from buildings/terrain triangulation
                tri = coords.BuildingsTerrainTriangulation;
                Z(inBldgsLimits) = computeTriangulationElevation(coords, ...
                    tri, inBldgsLats, inBldgsLons, outputHeightReference);
            end
            
            % Populate surface elevation of locations outside of buildings
            % region (or if populated with NaN above), using previously
            % cached values if available. Do not get public properties
            % because that would query terrain at all points which may be
            % unneeded.
            outsideBldgsLimits = ~inBldgsLimits | isnan(Z);
            if isGeoidReference && ~isempty(coords.pGroundHeightAboveGeoid)
                Z(outsideBldgsLimits) = coords.pGroundHeightAboveGeoid(outsideBldgsLimits,:);
            elseif ~isGeoidReference && ~isempty(coords.pGroundHeightAboveEllipsoid)
                Z(outsideBldgsLimits) = coords.pGroundHeightAboveEllipsoid(outsideBldgsLimits,:);
            else
                Z(outsideBldgsLimits) = terrain.internal.TerrainSource.queryTerrain(...
                    coords.TerrainSource, latlon(outsideBldgsLimits,:), outputHeightReference);
            end
        end
        
        function Z = computeGroundHeight(coords, outputHeightReference)
            %computeGroundHeight   Return antenna ground elevation at locations
            
            Z = terrain.internal.TerrainSource.queryTerrain(...
                coords.TerrainSource, coords.LatitudeLongitude, outputHeightReference);
        end
        
        function displatlon = computeDisplayLatitudeLongitude(coords)
            %computeDisplayLatitudeLongitude   Return site locations for display
            
            % Use small offset to avoid graphical artifacts produced when
            % plotting at north or south pole
            latlon = coords.LatitudeLongitude;
            lat = latlon(:,1);
            lat(lat == 90) = 89.999999;
            lat(lat == -90) = -89.999999;
            displatlon = [lat latlon(:,2)];
        end
    end
end
