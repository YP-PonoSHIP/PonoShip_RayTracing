function [az,el] = angle(sourceSites, targetSites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate sites
validateattributes(sourceSites,{'txsite','rxsite'},{'nonempty'}, ...
    'angle','',1);
validateattributes(targetSites,{'txsite','rxsite'},{'nonempty'}, ...
    'angle','',2);

% Process optional name/value pairs
p = inputParser;
addOptionalPath = nargin > 2 && mod(numel(varargin),2);
travelPath = 'euclidean';
if addOptionalPath
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('Path',travelPath,@(x)ischar(x)||isstring(x));
end
p.addParameter('Map', []);
p.addParameter('SourceAntennaSiteCoordinates', []);
p.addParameter('TargetAntennaSiteCoordinates', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.SourceAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(sourceSites, targetSites);
else
    usingCartesian = strcmp(p.Results.SourceAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

if addOptionalPath
    try
        travelPath = validatestring(p.Results.Path, {'euclidean','greatcircle'}, ...
            'angle','',3);
    catch me
        % Check if option is a match for deprecated "geodesic" option
        try
            travelPath = validatestring(p.Results.Path, {'geodesic'}, ...
                'angle','',3);
        catch
            rethrow(me)
        end
    end
end

% Get antenna site coordinates
if usingCartesian
    map = rfprop.internal.Validators.validateCartesianMap(p);
    coordSys = 'cartesian';
else
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'angle');
    coordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, coordSys);
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSites, map, 'angle');
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'angle');

% Compute angle
if strcmp(travelPath, 'euclidean')
    [X,Y,Z] = targetCoords.position(sourceCoords);
    [az,el] = cart2sph(X,Y,Z);
    az = rad2deg(az);
    el = rad2deg(el);
else
    % Validate sites are all geographic
    if usingCartesian
        error(message('shared_channel:rfprop:GreatCircleCartesianNotSupported'))
    end
    
    [d, heading] = greatCircleDistance(sourceSites, targetSites, ...
        'SourceAntennaSiteCoordinates', sourceCoords, ...
        'TargetAntennaSiteCoordinates', targetCoords);
    az = -heading + 90; % Convert clockwise-from-north to counter-clockwise-from-east
    
    % Detect where distance is zero, so angle should also be zero
    az((d == 0)) = 0;
    
    % Elevation is always zero for 2-D great circle
    el = zeros(size(az));
end

% Guarantee azimuth in range [-180,180]
az = wrapTo180(az);
end