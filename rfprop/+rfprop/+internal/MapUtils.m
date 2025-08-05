classdef (Sealed, Hidden) MapUtils
    %MapUtils   Reserved for MathWorks internal use only
    
    %   Copyright 2018-2021 The MathWorks, Inc. 
        
    methods(Static)
        function [lats, lons, res] = geogrid(latNorth, latSouth, lonEast, lonWest, res, isAutoRes, maxrange, maxImageSize, fcnName)
            % Return rectangular grid locations
            
            % Create lat/lon grid. The grid needs to match the format of MATLAB image
            % data, where the first row of data is the top edge of the image and the
            % first column of data is the left side of the image. Since the image will
            % be put on a map, north needs to be the first row and west needs to be the
            % first column:
            %  [NW ... NE]
            %  |   ...   |
            %  [SW ... SE]
            [lonv, latv, numGridEl] = rfprop.internal.MapUtils.geovec(res, latNorth, latSouth, lonEast, lonWest);
            
            % Protect against case where auto-generated resolution is too small for
            % geographic boundaries and grid size blows up, which can occur for
            % multiple tx that are distant. Increase resolution until grid size does
            % not exceed maximum. This is needed here because automatic resolution is 
            % computed only from maxrange (without geobounds).
            if isAutoRes                
                while numGridEl > rfprop.Constants.MaxNumSitesAutoResolutionGrid || ...
                        ((length(lonv) > maxImageSize) || (length(latv) > maxImageSize))
                    res = res + .1*res; % Increase resolution by 10%
                    [lonv, latv, numGridEl] = rfprop.internal.MapUtils.geovec(res, latNorth, latSouth, lonEast, lonWest);
                end
            end
            
            % Validate that data grid does not exceed maximum image size
            if (length(lonv) > maxImageSize) || (length(latv) > maxImageSize)
                resStr = mat2str(round(res,rfprop.Constants.MaxDistanceInfoDecimalPlaces));
                error(message('shared_channel:rfprop:ContourmapTooMuchData', fcnName, resStr))
            end
            
            % Validate that resolution is less than max range, or else there will
            % be no data to show
            if res >= min(maxrange)
                maxRangeStr = mat2str(round(min(maxrange),rfprop.Constants.MaxDistanceInfoDecimalPlaces));
                error(message('shared_channel:rfprop:ContourmapMaxRangeTooSmall', fcnName, maxRangeStr))
            end
            
            % Create sites to fill up geographical boundaries according to lat/lon
            % grid. Signal strengths will be computed at the sites that are within
            % range of each site.
            [lons, lats] = meshgrid(lonv, latv);
        end
        
        function [latNorth, latSouth, lonEast, lonWest, animation] = geobounds(latlon, maxrange, animation)
            %geobounds   Geographic boundaries for site range
            %   [LAT_N,LAT_S,LON_E,LON_W] = geobounds(LATLON, RANGES) returns north and
            %   south latitudes and east and west longitudes defining the rectangular
            %   region that contains RANGES distance from antenna sites. The LATLON
            %   and RANGES inputs must be the same length.
            
            earthR = rfprop.Constants.MinorEarthRadius;
            latNorth = [];
            latSouth = [];
            lonEast = [];
            lonWest = [];
            
            % Boundaries define a single window around all sites, which will be used to
            % create grid to sample signal strengths and also define the rectangle for
            % image overlay on a map. A single image is necessary so that coverage maps
            % of all txs can be merged together.
            for k = 1:size(latlon,1)
                % Get site location, assuming longitude has already been
                % wrapped to 180
                lat = latlon(k,1);
                lon = latlon(k,2);
                
                % Get max range as angle of Earth, using minor axis radius of WGS-84,
                % which will slightly overestimate result and guarantee grid contains
                % all points in range
                maxrangeRad = maxrange(k) / earthR;
                maxrangeDeg = rad2deg(maxrangeRad);
                
                % Bounding latitudes are along north/south meridian, saturating at poles
                siteLatSouth = lat - maxrangeDeg;
                if siteLatSouth < -90
                    siteLatSouth = -90;
                end
                siteLatNorth = lat + maxrangeDeg;
                if siteLatNorth > 90
                    siteLatNorth = 90;
                end
                
                if (siteLatSouth <= -90) || (siteLatNorth >= 90)
                    % If bounds include a pole, use full range of longitude
                    siteLonWest = -180;
                    siteLonEast = 180;
                else
                    % Bounding longitudes are NOT simply along east/west latitude but need
                    % to be computed according to tangent longitudes to edge of circle
                    % defined by maxrange
                    % Reference: http://gis.stackexchange.com/questions/19221/find-tangent-point-on-circle-furthest-east-or-west
                    deltaLon = rad2deg(asin(sin(maxrangeRad)/cos(deg2rad(lat))));
                    siteLonWest = lon - deltaLon;
                    siteLonEast = lon + deltaLon;
                    
                    % If longitude outside of [-180,180] range, then date line has been
                    % crossed, in which case use full range of longitude. The values cannot
                    % simply be converted to wrapped values (e.g. [-190,-170] =>
                    % [170,-170]), since that would transform the west value to an east
                    % value and cause image overlay on a map to use the wrong area.
                    if siteLonWest < -180 || siteLonEast > 180 || ~isreal(siteLonWest) || ~isreal(siteLonEast)
                        siteLonWest = -180;
                        siteLonEast = 180;
                        animation = 'none'; % Turn off animation since image wrapped around Earth
                    end
                end
                
                % Assign outputs, where boundaries defined by outermost boundary for
                % any site
                if isempty(latNorth) || (siteLatNorth > latNorth)
                    latNorth = siteLatNorth;
                end
                if isempty(latSouth) || (siteLatSouth < latSouth)
                    latSouth = siteLatSouth;
                end
                if isempty(lonEast) || (siteLonEast > lonEast)
                    lonEast = siteLonEast;
                end
                if isempty(lonWest) || (siteLonWest < lonWest)
                    lonWest = siteLonWest;
                end
            end
        end
        
        function [lats, lons] = georange(sites, lats, lons, datarange, terrainSource)
            % Return locations within range of sites
            
            % Return early if no sites
            if isempty(lats) || isempty(lons)
                return
            end
            
            % Compute great circle range
            gc = nan(numel(lats),numel(sites));
            for siteInd = 1:numel(sites)
                site = sites(siteInd);
                lat = double(site.Latitude);
                lon = wrapTo180(double(site.Longitude));
                gc(:,siteInd) = rfprop.internal.MapUtils.greatCircleDistance(lat, lon, lats, lons);
            end
            
            % Apply data range
            isInDataRange = gc <= datarange;
            
            % Apply terrain range
            if ~strcmp(terrainSource,'none')
                latLimits = terrainSource.LatitudeLimits;
                lonLimits = terrainSource.LongitudeLimits;
                if ~isequal(latLimits,[-90 90]) || ~isequal(lonLimits, [-180 180])
                    outOfTerrainRange = terrainSource.isOutOfRange(lats,lons);
                    isInDataRange = isInDataRange & ~outOfTerrainRange;
                end
            end
            
            % Trim locations to those which are within range
            isInDataRange = any(isInDataRange,2);
            lats = lats(isInDataRange);
            lons = lons(isInDataRange);
        end
        
        function [lonv, latv, numGridEl] = geovec(res, latNorth, latSouth, lonEast, lonWest)
            % Return geographic vectors for data grid
            
            % Convert resolution from geodesic meters to degrees.
            
            % In Longitudinal direction (North-South, or y-direction), all distances
            % cover the same length in degrees since all lines of Longitude are great
            % circles. The length in this direction defines the angular resolution for
            % Latitudes.
            earthR = rfprop.Constants.SphericalEarthRadius;
            resDegLat = rad2deg(res/earthR);
            
            % In Latitudinal direction (East-West, or x-direction), length in degrees
            % varies with location since lines of Latitude are not great circles but
            % have varying radius. Use the average Latitude covered in the range for
            % the conversion radius. The length in this direction defines the angular
            % resolution for Longitudes.
            if (latNorth >= 0) && (latSouth <= 0) % Range goes through Equator
                latRadius = earthR;
            else
                largestLat = mean([abs(latNorth), abs(latSouth)]);
                latRadius = cosd(largestLat) * earthR;
            end
            resDegLon = rad2deg(res/latRadius);
            
            % Use linspace to guarantee endpoints of grid are the geo boundaries
            numlonv = floor(abs(diff([lonWest, lonEast])) / resDegLon);
            numlatv = floor(abs(diff([latNorth, latSouth])) / resDegLat);
            lonv = linspace(lonWest, lonEast, numlonv); % Longitude is angle in x-direction
            latv = linspace(latNorth, latSouth, numlatv); % Latitude is angle in y-direction
            numGridEl = numel(lonv)*numel(latv);
        end
        
        function [d, az] = greatCircleDistance(lat1, lon1, lat2, lon2)
            % Return great circle distance and angle between points using spherical Earth
            
            radius = rfprop.Constants.SphericalEarthRadius;
            [d, az] = map.geodesy.internal.greatCircleDistance(lat1, lon1, lat2, lon2, radius);
        end
        
        function [latfwd, lonfwd] = greatCircleForward(lat1, lon1, d, az)
            % Return location moving from a point on a spherical Earth
            
            az = mod(-az + 90, 360); % Convert angle to clockwise from north
            radius = rfprop.Constants.SphericalEarthRadius;
            [latfwd,lonfwd] = map.geodesy.internal.greatCircleTrace(lat1, lon1, d(:), az(:), radius);
            lonfwd = wrapTo180(lonfwd);
        end
        
        function [latout,lonout,actres,d] = sampleGreatCircle(lat1, lon1, lat2, lon2, res)
            % Return points sampled along a great circle
            
            radius = rfprop.Constants.SphericalEarthRadius;
            if isscalar(lat2)
                [latout,lonout,S] = map.geodesy.internal.sampleGreatCircleArc(lat1, lon1, lat2, lon2, res, radius);
                actres = S(2); % Resolution is distance to second sample
                d = S(end); % Total distance is cumulative distance to last sample
            else
                numLocations = numel(lat2);
                latout = cell(numLocations,1);
                lonout = cell(numLocations,1);
                actres = zeros(numLocations,1);
                d = zeros(numLocations,1);
                for k = 1:numLocations
                    [latout{k},lonout{k},S] = map.geodesy.internal.sampleGreatCircleArc(lat1, lon1, lat2(k), lon2(k), res, radius);
                    actres(k) = S(2);
                    d(k) = S(end);
                end
            end
        end
        
        function [X,Y,Z] = geodetic2enu(lat0, lon0, h0, lat, lon, h, spheroid)
            % Return relative ENU coordinates
            
            if nargin < 7
                spheroid = map.geodesy.internal.wgs84Spheroid;
            end
            [X,Y,Z] = geodetic2enu(lat, lon, h, lat0, lon0, h0, spheroid);
        end

        function [lat,lon,h] = enu2geodetic(lat0, lon0, h0, pos, spheroid)
            % Return geodetic coordinates from relative ENU coordinates
            
            if nargin < 5
                spheroid = map.geodesy.internal.wgs84Spheroid;
            end
            [lat,lon,h] = enu2geodetic(pos(1,:)',pos(2,:)',pos(3,:)', ...
                lat0, lon0, h0, spheroid);
        end
        
        function rotMtx = enuRotation(lat1, lon1, height1, lat2, lon2, height2)
            % Rotation matrix from ENU of the 2nd point into the ENU of the 1st point
            
            spheroid = map.geodesy.internal.wgs84Spheroid;
            
            % Calculate lat/lon/height of axes points ([1;0;0],[0;1;0],[0;0,1]) in ENU
            % of 2nd point
            [latXYZ, lonXYZ, heightXYZ] = arrayfun( ...
                @(x,y,z)rfprop.internal.MapUtils.enu2geodetic( ...
                x, y, z, eye(3), spheroid), ...
                lat2, lon2, height2, 'UniformOutput', false);
            
            % Calculate [X;Y;Z] of 2nd point in the ENU of 1st point
            [x0, y0, z0] = arrayfun( ...
                @(x,y,z)rfprop.internal.MapUtils.geodetic2enu( ...
                lat1, lon1, height1, x,y,z, spheroid), ...
                lat2, lon2, height2, 'UniformOutput', false);
            
            % Calculate [X;Y;Z] of axes points in the ENU of 1st point
            [x1, y1, z1] = cellfun( ...
                @(x,y,z)rfprop.internal.MapUtils.geodetic2enu( ...
                lat1, lon1, height1, x, y, z, spheroid), ...
                latXYZ, lonXYZ, heightXYZ, 'UniformOutput', false);
            
            % The difference between 2nd point's axes and 2nd point's
            % origin in the ENU of 1st point gives the rotation matrix.
            rotMtx = cellfun( ...
                @(x0,y0,z0,x1,y1,z1)([x1'-x0'; y1'-y0'; z1'-z0']), ...
                x0, y0, z0, x1,y1,z1, 'UniformOutput', false);            
        end
    end
end
