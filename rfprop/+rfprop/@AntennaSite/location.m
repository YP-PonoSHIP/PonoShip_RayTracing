function varargout = location(sites, dist, az)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate number of output arguments
nargoutchk(0,2);

% Validate first input
validateattributes(sites,{'txsite','rxsite'},{'nonempty'}, ...
    'location','',1);
coords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(...
    sites, 'none');
usingCartesian = strcmp(coords.CoordinateSystem,'cartesian');

% Validate second and third inputs. Note that if second input is specified,
% third must also be specified.
if nargin > 1
    if usingCartesian
        error(message('shared_channel:rfprop:LocationGreatCircleCartesianNotSupported'))
    end
    
    validateattributes(dist,{'numeric'},{'vector','real','finite','nonnegative','nonsparse'}, ...
        'location','',2);
    
    % Impose limit on distance using half-circumference of Earth
    maxDistance = rfprop.Constants.SphericalEarthRadius * pi;
    if dist > maxDistance
        error(message('shared_channel:rfprop:LocationDistanceGreaterThanMax', ...
            sprintf('%.0f',maxDistance/1000)));
    end
    
    % Perform scalar expansion
    if isscalar(dist) && ~isscalar(az)
        dist = dist + zeros(numel(az),1);
    elseif isscalar(az) && ~isscalar(dist)
        az = az + zeros(numel(dist),1);
    end
    
    validateattributes(az,{'numeric'},{'vector','real','finite','nonsparse','numel',numel(dist)}, ...
        'location','',3);
    
    % Validate scalar site if passing non-scalar distance/azimuths
    if numel(dist) > 1
        validateattributes(sites,{'txsite','rxsite'},{'scalar'}, ...
            'location','',1);
    end
end

varargout = cell(1,nargout);
if usingCartesian
    nargoutchk(0,1)
    varargout{1} = coords.AntennaPosition;
else
    % Get Latitude and Longitude values
    latlon = coords.LatitudeLongitude;
    lats = latlon(:,1);
    lons = latlon(:,2);
    
    % Compute output locations
    if nargin >= 2
        [lats,lons] = rfprop.internal.MapUtils.greatCircleForward(lats, lons, dist, az);
    end
    
    % Assign outputs
    if nargout <= 1 % Can be 0 to assign into ans
        varargout{1} = [lats lons];
    else
        varargout{1} = lats;
        varargout{2} = lons;
    end
end