function d = distance(sourceSites, targetSites, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate sites
validateattributes(sourceSites,{'txsite','rxsite'},{'nonempty'}, ...
    'distance','',1);
validateattributes(targetSites,{'txsite','rxsite'},{'nonempty'}, ...
    'distance','',2);

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
            'distance','',3);
    catch me
        % Check if option is a match for deprecated "geodesic" option
        try
            travelPath = validatestring(p.Results.Path, {'geodesic'}, ...
                'distance','',3);
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
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'distance');
    coordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, coordSys);
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSites, map, 'distance');
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'distance');

% Compute distance
if strcmp(travelPath, 'euclidean')    
    [X,Y,Z] = targetCoords.position(sourceCoords);
    [~,~,d] = cart2sph(X,Y,Z);
else
    % Validate sites are all geographic
    if usingCartesian
        error(message('shared_channel:rfprop:GreatCircleCartesianNotSupported'))
    end
    
    % test if any antipodal pairs
    d = greatCircleDistance(sourceSites, targetSites, ...
        'SourceAntennaSiteCoordinates', sourceCoords, ...
        'TargetAntennaSiteCoordinates', targetCoords);
end
end