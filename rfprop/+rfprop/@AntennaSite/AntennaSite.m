classdef (Abstract, Hidden) AntennaSite < handle & matlab.mixin.Copyable & matlab.mixin.CustomDisplay
    %AntennaSite   Abstract class for antenna sites
    
    %   Copyright 2017-2021 The MathWorks, Inc.
    
    properties(Dependent)
        %Name   Site name
        %   Site name, specified as a character vector or string.
        Name(1,:) {validateName}
        %CoordinateSystem   Coordinate system of site location
        %   Coordinate system of site location, specified as 'geographic' or
        %   'cartesian'. If specified as 'geographic', the site location is defined
        %   using Latitude, Longitude, and AntennaHeight. If specified as
        %   'cartesian', the site location is defined using AntennaPosition. The
        %   default values is 'geographic'.
        CoordinateSystem(1,:) char
        %Latitude   Site latitude coordinates
        %   Site latitude coordinates in decimal degrees, specified as a numeric
        %   scalar in range -90 to 90. Coordinates are defined with reference to
        %   Earth ellipsoid model WGS-84. This property applies when
        %   CoordinateSystem is set to 'geographic'.
        Latitude(1,1) {validateLatitude}
        %Longitude   Site longitude coordinates
        %   Site longitude coordinates in decimal degrees, specified as a numeric
        %   scalar. Coordinates are defined with reference to Earth ellipsoid model
        %   WGS-84. This property applies when CoordinateSystem is set to
        %   'geographic'.
        Longitude(1,1) {validateLongitude}
        %Antenna   Antenna element or array object
        %   Antenna element or array, specified as an object or 'isotropic'. The
        %   default value is 'isotropic', which corresponds to an antenna that
        %   radiates uniformly in all directions.
        Antenna(1,:)
        %AntennaPosition   Antenna position Cartesian coordinates (m)
        %   Position of antenna center, specified as a 3-by-1 vector representing
        %   [x;y;z] Cartesian coordinates (in meters). This property applies when
        %   CoordinateSystem is set to 'cartesian'. The default value is [0;0;0].
        AntennaPosition(:,1) {validateAntennaPosition}
        %AntennaAngle   Angle of antenna's local x-axis (degrees)
        %   Angle of antenna's local Cartesian coordinate system x-axis, specified
        %   as a numeric scalar representing azimuth angle or as a 2-by-1 vector
        %   representing azimuth and elevation angles. Azimuth angle is measured
        %   counter-clockwise to the antenna's x-axis, either from east (for
        %   geographic sites), or from the global x-axis around the global z-axis
        %   (for Cartesian sites). Elevation angle is measured from the horizontal
        %   (or X-Y) plane to the antenna's x-axis in the range -90 to 90. The
        %   default value is 0.
        AntennaAngle(:,1) {validateAntennaAngle}
        %AntennaHeight   Antenna height above surface (m)
        %   Height of antenna above ground or building surface, specified as a
        %   numeric non-negative scalar. If the site coincides with a building, the
        %   height is measured in meters from the top of the building to the center
        %   of the antenna. Otherwise, the height is measured in meters from ground
        %   elevation to the center of the antenna. This property applies when
        %   CoordinateSystem is set to 'geographic'. The default value is 10 for
        %   txsite and 1 for rxsite.
        AntennaHeight(1,1) {validateAntennaHeight}
        %SystemLoss   System loss (dB)
        %   System loss in dB, specified as non-negative numeric scalar. System
        %   loss may include transmission line loss and any other miscellaneous
        %   system losses. The default value is 0.
        SystemLoss(1,1) {validateSystemLoss}
    end
    
    properties(Dependent, Access = {?rfprop.AntennaSite, ?rfprop.PropagationModel, ?comm.Ray})
        UID % Unique identifier string
    end
    
    properties(Access = private)
        pCoordinateSystem = 'geographic'
        pLatitude = 42.3001
        pLongitude = -71.3504
        pAntennaPosition = [0;0;0]
        pAntennaAngle = 0
        pAntennaHeight = 10
        pSystemLoss = 0
    end
    
    % Persistent instance-specific properties (do not copy)
    properties(NonCopyable, Access = private)
        pName
        pAntenna
        pIsCopyWithDefaultName = false
    end
    
    properties(Hidden, Constant, NonCopyable)
        SiteGraphicsTemplate = struct("marker", struct("Location",[NaN, NaN, NaN], "ID", NaN), ...
            "pattern", [], ...
            "los", struct(), ...
            "link", struct(), ...
            "contour", [], ...
            "rays", struct(), ...
            "legend", [], ...
            "infoboxLegend", [])
    end
    
    properties(Constant, Hidden)
        CoordinateSystemChoices = {'geographic','cartesian'}
    end
    
    properties(Hidden, SetAccess = protected)
        IconSize = [36,36]
        IconAlignment = 'top'
        Icon = []
        ClusterMarkers = false
        ShowAntennaHeight = 'auto'
    end
    
    % Transient instance-specific properties (do not copy or load)
    properties(NonCopyable, Transient, Access = private)
        pUID = initLock
        pCopyNameCounter = 1
    end
    
    methods(Abstract, Access = protected)
        props = getSiteSpecificProperties(site)
    end
    
    methods
        function set.Name(site, name)
            if isstring(name)
                name = char(name);
            end
            site.pName = name;
            
            % Reset variables used to generate Name for copy
            site.pIsCopyWithDefaultName = false;
            site.pCopyNameCounter = 1;
        end
        
        function name = get.Name(site)
            name = site.pName;
            
            % For performance, generate char Name on first access
            if isnumeric(name) % Initialized to number
                name = message('shared_channel:rfprop:SiteName', name).getString;
                site.pName = name;
            end
        end
        
        function set.CoordinateSystem(site, coordSystem)
            site.pCoordinateSystem = validateCoordinateSystem(coordSystem);
        end
        
        function coordSystem = get.CoordinateSystem(site)
            coordSystem = site.pCoordinateSystem;
        end
        
        function set.Latitude(site, lat)
            site.pLatitude = lat;
        end
        
        function lat = get.Latitude(site)
            lat = site.pLatitude;
        end
        
        function set.Longitude(site, lon)
            site.pLongitude = lon;
        end
        
        function lon = get.Longitude(site)
            lon = site.pLongitude;
        end
        
        function set.Antenna(site, ant)
            if ischar(ant) || isstring(ant)
                ant = validatestring(ant,{'isotropic'},'','Antenna');
            else
                validateAntennaObject(ant,'scalar')
                site.validateAntennaTransmitterFrequency(ant)
            end
            site.pAntenna = ant;
        end
        
        function ant = get.Antenna(site)
            ant = site.pAntenna;
        end
        
        function set.AntennaPosition(site, pos)
            site.pAntennaPosition = pos;
        end
        
        function pos = get.AntennaPosition(site)
            pos = site.pAntennaPosition;
        end
        
        function set.AntennaAngle(site, ang)
            site.pAntennaAngle = ang;
        end
        
        function ang = get.AntennaAngle(site)
            ang = site.pAntennaAngle;
        end
        
        function set.AntennaHeight(site, ht)
            site.pAntennaHeight = ht;
        end
        
        function ht = get.AntennaHeight(site)
            ht = site.pAntennaHeight;
        end
        
        function set.SystemLoss(site, sloss)
            site.pSystemLoss = sloss;
        end
        
        function sloss = get.SystemLoss(site)
            sloss = site.pSystemLoss;
        end
        
        function uid = get.UID(site)
            uid = site.pUID;
            
            % For performance, generate char UID on first access
            if isnumeric(uid) % Initialized to number
                uid = sprintf('site%u', uid);
                site.pUID = uid;
            end
        end
        
        function sites = horzcat(varargin)
            %horzcat   Overload horizontal concatenation
            
            % All inputs must be same class
            siteClass = class(varargin{1});
            for sitesInd = 1:numel(varargin)
                validateattributes(varargin{sitesInd},{siteClass},{});
            end
            
            % All inputs must have same CoordinateSystem
            rfprop.internal.Validators.validateCoordinateSystem(varargin{:});
            sites = builtin('horzcat',varargin{:});
        end
        
        function sites = vertcat(varargin)
            %vertcat   Overload vertical concatenation
            
            % All inputs must be same class
            siteClass = class(varargin{1});
            for sitesInd = 1:numel(varargin)
                validateattributes(varargin{sitesInd},{siteClass},{});
            end
            
            % All inputs must have same CoordinateSystem
            rfprop.internal.Validators.validateCoordinateSystem(varargin{:});
            sites = builtin('vertcat',varargin{:});
        end
    end
    
    methods(Access = protected)
        function group = getPropertyGroups(site)
            %getPropertyGroups   Get property groups for object display
            
            showGeographicProps = isempty(site) || strcmp(site(1).CoordinateSystem,'geographic');
            if showGeographicProps
                propList = {
                    'Name','CoordinateSystem','Latitude','Longitude','Antenna', ...
                    'AntennaAngle','AntennaHeight','SystemLoss'};
            else
                propList = {
                    'Name','CoordinateSystem','Antenna','AntennaPosition', ...
                    'AntennaAngle','SystemLoss'};
            end
            propList = [propList, getSiteSpecificProperties(site)];
            group = matlab.mixin.util.PropertyGroup(propList);
        end
        
        function newSite = copyElement(site)
            %copyElement   Copy antenna site
            
            % Make a shallow copy of all copyable properties
            newSite = copyElement@matlab.mixin.Copyable(site);
            
            % Make deep copy of Antenna if using object
            ant = site.Antenna;
            if ischar(ant) || rfprop.AntennaSite.isCommAntenna(ant)
                newSite.Antenna = ant;
            elseif rfprop.AntennaSite.isPhasedAntenna(ant)
                newSite.Antenna = clone(ant);
            else
                newSite.Antenna = copy(ant);
            end
            
            % Set Name
            copyNameNum = site.pCopyNameCounter;
            site.pCopyNameCounter = copyNameNum + 1;
            if site.pIsCopyWithDefaultName
                % Generate Name using scheme: "NAME(1)", "NAME(2)",
                % "NAME(3)", etc. Since NAME already includes "Copy", this
                % helps keep generated Name succinct.
                newName = message('shared_channel:rfprop:SiteCopyOfCopyName', ...
                    site.Name, copyNameNum).getString;
            else
                % Generate Name using scheme: "NAME (Copy)", "NAME (Copy
                % 2)", "NAME (Copy 3)", etc.
                if copyNameNum < 2
                    newName = message('shared_channel:rfprop:SiteFirstCopyName', ...
                        site.Name).getString;
                else
                    newName = message('shared_channel:rfprop:SiteRepeatCopyName', ...
                        site.Name, copyNameNum).getString;
                end
            end
            newSite.Name = newName;
            
            % Set private properties            
            newSite.pIsCopyWithDefaultName = true;            
            setSiteIdentifiers(newSite);
        end
        
        function validateAntennaTransmitterFrequency(varargin)
            % Allow txsite to validate frequency when updating Antenna
        end
        
        function validateNumColumns(sites, v, prop)
            %validateNumColumns   Validate number of columns of property value
            
            try
                if ischar(v) % Treat character vector as scalar text
                    numCols = 1;
                else
                    numCols = size(v,2);
                end
                
                % Validate that v has either 1 column or same number of
                % columns as number of sites
                numSites = numel(sites);
                if ((numCols > 1) && (numCols ~= numSites))
                    error(message('shared_channel:rfprop:InvalidPropertyLength', prop, numSites));
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function inputs = validateSiteProperties(sites, inputs)
            %validateSiteProperties   Validate site properties specified in inputs
            
            try
                % Validate all property values. For performance, only validate
                % non-default values.
                usingDefaults = inputs.UsingDefaults;
                if ~ismember('Name', usingDefaults)
                    name = inputs.Name;
                    sites.validateNumColumns(name, 'Name')
                    % Name validation enforces "scalartext", so validate each
                    % element if cellstr or string array
                    if iscellstr(name)
                        validateNameArray(name)
                        for k = 1:numel(name)
                            validateName(name{k});
                        end
                    elseif isstring(name)
                        validateNameArray(name)
                        for k = 1:numel(name)
                            validateName(name(k));
                        end
                    else
                        % Must be character vector
                        validateName(name);
                    end
                end
                if ~ismember('CoordinateSystem', usingDefaults)
                    inputs.CoordinateSystem = validateCoordinateSystem(inputs.CoordinateSystem);
                end
                if ~ismember('Latitude', usingDefaults)
                    lat = inputs.Latitude;
                    sites.validateNumColumns(lat, 'Latitude')
                    validateLatitude(lat)
                end
                if ~ismember('Longitude', usingDefaults)
                    lon = inputs.Longitude;
                    sites.validateNumColumns(lon, 'Longitude')
                    validateLongitude(lon)
                end
                if ~ismember('Antenna', usingDefaults)
                    ant = inputs.Antenna;
                    sites.validateNumColumns(ant, 'Antenna')
                    % Antenna validation enforces object type, so validate each
                    % element if cell (size was verified in validateNumColumns)
                    if isstring(ant)
                        ant = cellstr(ant);
                    end
                    if iscell(ant)
                        validateAntennaCellArray(ant)
                        for k = 1:numel(ant)
                            antk = ant{k};
                            if ischar(antk) || isstring(antk)
                                ant{k} = validatestring(antk,{'isotropic'},'','Antenna');
                            else
                                validateAntennaObject(antk,'scalar');
                            end
                        end
                        inputs.Antenna = ant;
                    elseif ischar(ant)
                        inputs.Antenna = validatestring(ant,{'isotropic'},'','Antenna');
                    else
                        validateAntennaObject(ant,'row');
                    end
                end
                if ~ismember('AntennaPosition', usingDefaults)
                    pos = inputs.AntennaPosition;
                    sites.validateNumColumns(pos, 'AntennaPosition')
                    validateAntennaPosition(pos)
                end
                if ~ismember('AntennaAngle', usingDefaults)
                    ang = inputs.AntennaAngle;
                    sites.validateNumColumns(ang, 'AntennaAngle')
                    validateAntennaAngle(ang)
                end
                if ~ismember('AntennaHeight', usingDefaults)
                    ht = inputs.AntennaHeight;
                    sites.validateNumColumns(ht, 'AntennaHeight')
                    validateAntennaHeight(ht)
                end
                if ~ismember('SystemLoss', usingDefaults)
                    sloss = inputs.SystemLoss;
                    sites.validateNumColumns(sloss, 'SystemLoss')
                    validateSystemLoss(sloss)
                end
            catch e
                throwAsCaller(e)
            end
        end
        
        function setSiteProperties(sites, inputs)
            %setSiteProperties   Set site properties from inputs
            
            % Get inputs
            usingDefaults = inputs.UsingDefaults;
            name = inputs.Name;
            coordSystem = inputs.CoordinateSystem;
            lat = inputs.Latitude;
            lon = inputs.Longitude;
            ant = inputs.Antenna;
            pos = inputs.AntennaPosition;
            ang = inputs.AntennaAngle;
            ht = inputs.AntennaHeight;
            sloss = inputs.SystemLoss;
            
            % Set properties on each site. For performance, array values
            % are validated once (in validateSiteProperties) and set into
            % private fields, thereby preventing repeated validation for
            % each site. Most values need to be set since default values
            % are defined in inputParser.
            usingDefaultName = ismember('Name',usingDefaults);
            isScalarName = isscalar(name) || ischar(name);
            isScalarLat = isscalar(lat);
            isScalarLon = isscalar(lon);
            isScalarAnt = isscalar(ant) || ischar(ant);
            isScalarPos = size(pos,2) == 1; % Treat [x;y;z] vector as scalar
            isScalarAng = size(ang,2) == 1; % Treat [az;el] vector as scalar
            isScalarHt = isscalar(ht);
            isScalarSLoss = isscalar(sloss);
            for k = 1:numel(sites)
                % For Name, an optimization is to generate default on first
                % get, so only set value here if non-default. Convert all
                % values to char to guarantee char vector stored.
                if ~usingDefaultName
                    if iscell(name)
                        if isScalarName
                            sites(k).pName = char(name{1});
                        else
                            sites(k).pName = char(name{k});
                        end
                    else
                        if isScalarName
                            sites(k).pName = char(name);
                        else
                            sites(k).pName = char(name(k));
                        end
                    end
                end
                sites(k).pCoordinateSystem = coordSystem;
                if isScalarLat
                    sites(k).pLatitude = lat;
                else
                    sites(k).pLatitude = lat(k);
                end
                if isScalarLon
                    sites(k).pLongitude = lon;
                else
                    sites(k).pLongitude = lon(k);
                end
                if iscell(ant)
                    if isScalarAnt
                        sites(k).pAntenna = ant{1};
                    else
                        sites(k).pAntenna = ant{k};
                    end
                else
                    if isScalarAnt
                        sites(k).pAntenna = ant;
                    else
                        sites(k).pAntenna = ant(k);
                    end
                end
                if isScalarPos
                    sites(k).pAntennaPosition = pos;
                else
                    sites(k).pAntennaPosition = pos(:,k); % Column index to allow [x;y;z] spec
                end
                if isScalarAng
                    sites(k).pAntennaAngle = ang;
                else
                    sites(k).pAntennaAngle = ang(:,k); % Column index to allow [az;el] spec
                end
                if isScalarHt
                    sites(k).pAntennaHeight = ht;
                else
                    sites(k).pAntennaHeight = ht(k);
                end
                if isScalarSLoss
                    sites(k).pSystemLoss = sloss;
                else
                    sites(k).pSystemLoss = sloss(k);
                end
                
                % Initialize other NonCopyable properties. This is required
                % since array creation calls the constructor once and then
                % copies the object to each element of the array, leaving
                % these properties as empty [].
                sites(k).pIsCopyWithDefaultName = false;
                sites(k).pCopyNameCounter = 1;
            end
        end
        
        function setSiteIdentifiers(sites, inputs)
            %setSiteIdentifiers   Set site identifier properties
            
            usingDefaultName = (nargin > 1) && ismember('Name', inputs.UsingDefaults);
            for k = 1:numel(sites)
                % Set private Name and UID fields to site numbers. For
                % performance, defer creation of corresponding strings
                % until requested.
                if usingDefaultName
                    % Only assign name number if using default
                    sites(k).pName = rfprop.AntennaSite.siteNumber('defaultname','next');
                end
                sites(k).pUID = rfprop.AntennaSite.siteNumber('id','next');
            end
        end
    end
    
    methods
        show(sites, varargin)
        hide(sites, varargin)
        [az, el] = angle(sourceSites, targetSites, travelPath, varargin)
        d = distance(sourceSites, targetSites, travelPath, varargin)
        varargout = location(site, dist, az, varargin)
        varargout = los(observerSite, targetSites, varargin)
        Z = elevation(sites, varargin)
    end
    
    methods(Hidden)
        [d,az,el] = distanceangle(sourceSites, targetSites, varargin)
        [X,Y,Z] = position(sourceSites, targetSites, varargin)
        [d, heading] = greatCircleDistance(sourceSites, targetSites, varargin)
        [lats, lons, actualResolutions, d] = sampleGreatCircle(sourceSites, targetSites, resolution, varargin)
        G = gain(site, fq, az, el, varargin)
        [G,az,el] = gainPattern(site,fq)
        Z = terrain(sites, varargin)
    end
    
    methods(Static, Hidden)
        function sitenum = siteNumber(counter, newSitenum)
            %siteNumber   Site number for default name or ID
            
            % Use nameCounter to track number for default Name property and
            % idCounter to track number for private site ID. These numbers
            % are not necessarily the same, since the nameCounter only
            % increments when a default name is generated, whereas
            % idCounter increments for every site created.
            persistent nameCounter;
            persistent idCounter;
            if isempty(nameCounter)
                nameCounter = 0;
            end
            if isempty(idCounter)
                idCounter = 0;
            end
            
            % Set counter value to input or else increment if 'next'
            if nargin > 1
                if ischar(newSitenum) && strcmp(newSitenum,'next')
                    if strcmp(counter, 'defaultname')
                        nameCounter = nameCounter + 1;
                    else
                        idCounter = idCounter + 1;
                    end
                else
                    if strcmp(counter, 'defaultname')
                        nameCounter = newSitenum;
                    else
                        idCounter  = newSitenum;
                    end
                end
            end
            
            % Return counter value
            if strcmp(counter, 'defaultname')
                sitenum = nameCounter;
            else
                sitenum = idCounter;
            end
        end
    end
    
    methods(Static, Access=protected)        
        function inputs = parseInputs(p, varargin)
            %parseInputs   Parse inputs
            
            try
                % Get inputs from parser
                p.parse(varargin{:});
                inputs = p.Results;
                
                % Convert geographic site properties specified as column vectors
                name = inputs.Name;
                if ~isscalar(name) && ~ischar(name) && iscolumn(name)
                    inputs.Name = name';
                end
                lat = inputs.Latitude;
                if ~isscalar(lat) && iscolumn(lat)
                    inputs.Latitude = lat';
                end
                lon = inputs.Longitude;
                if ~isscalar(lon) && iscolumn(lon)
                    inputs.Longitude = lon';
                end
                
                % Store parameters that are using defaults
                inputs.UsingDefaults = p.UsingDefaults;
            catch e
                throwAsCaller(e)
            end
        end
        
        function numSites = numSites(inputs)
            %numSites   Number of sites specified in inputs
            
            % Name specified as char vector is one site
            name = inputs.Name;
            if ischar(name) 
                numNames = 1;
            else
                numNames = size(name,2);
            end
            
            % Antenna specified as char vector is one site
            ant = inputs.Antenna;
            if ischar(ant)
                numAntennas = 1;
            else
                numAntennas = size(inputs.Antenna,2);
            end
            
            % Number of sites is maximum number of columns of any input
            numSites = max([numNames, ...
                size(inputs.Latitude,2), size(inputs.Longitude,2), ...
                numAntennas, size(inputs.AntennaPosition,2), size(inputs.AntennaAngle,2), ...
                size(inputs.AntennaHeight,2), size(inputs.SystemLoss,2)]);
        end
        
        function isPhased = isPhasedAntenna(ant)
            isPhased = isa(ant,'phased.internal.AbstractAntennaElement') || ...
                isa(ant,'phased.internal.AbstractArray') || ...
                isa(ant,'phased.internal.AbstractSubarray');
        end
        
        function isEM = isElectromagneticAntenna(ant)
            isEM = isa(ant,'em.Antenna') || isa(ant,'em.Array') || ...
                isa(ant,'installedAntenna') || isa(ant,'customAntennaStl');
        end
        
        function isComm = isCommAntenna(ant)
            isComm = isa(ant, 'arrayConfig');
        end
    end
end

function v = initLock

% Protect persistent variables from clear. Call mlock here so that it is
% called just once per session.
mlock;
v = 0;
end

function validateNameArray(name)
validateattributes(name, {'cell','string'}, {'vector'}, '', 'Name');
end

function validateName(name)
% Use scalartext attribute, which does not allow string array
% or cell array validation
validateattributes(name, {'char','string'}, {'scalartext'}, '', 'Name');
end

function coordSystem = validateCoordinateSystem(coordSystem)

coordSystem = validatestring(coordSystem, ...
    rfprop.AntennaSite.CoordinateSystemChoices,'','CoordinateSystem');
end

function validateLatitude(lat)

% Note that >= and <= imply real, positive, finite
validateattributes(lat, {'numeric'}, ...
    {'real','nonsparse','vector','>=',-90,'<=',90}, '', 'Latitude');
end

function validateLongitude(lon)

% While Latitude is constrained to range [-90,90], Longitude is
% unconstrained and allowed to wrap
validateattributes(lon, {'numeric'}, ...
    {'real','finite','nonsparse','vector'}, '', 'Longitude');
end

function validateAntennaCellArray(name)
validateattributes(name, {'cell'}, {'row'}, '', 'Antenna');
end

function validateAntennaObject(ant,expSize)

% Validate object type before validateattributes in order to
% get good error message
if ~rfprop.AntennaSite.isElectromagneticAntenna(ant) && ...
        ~rfprop.AntennaSite.isPhasedAntenna(ant) && ...
        ~rfprop.AntennaSite.isCommAntenna(ant) 
    error(message('shared_channel:rfprop:InvalidAntennaObject', 'Antenna'));
end

if nargin < 2
    expSize = 'scalar';
end
validateattributes(ant, {'handle','arrayConfig'}, {expSize}, '', 'Antenna');
end

function validateAntennaAngle(ang)

% Value may be one value [az] or two [az; el]
isAzimuth = (size(ang,1) == 1);
if isAzimuth
    nrows = 1;
else
    nrows = 2;
end
validateattributes(ang, {'numeric'}, ...
    {'real','finite','nonsparse','nrows',nrows}, '', 'AntennaAngle');

% Enforce elevation is within vertical range [-90,90]. Azimuth
% does not need validation since it is allowed to wrap.
if ~isAzimuth
    el = ang(2,:);
    if any(el < -90) || any(el > 90)
        error(message('shared_channel:rfprop:InvalidAntennaElevationRange', 'AntennaAngle'));
    end
end
end

function validateAntennaPosition(pos)

validateattributes(pos, {'numeric'}, ...
    {'real','finite','nonsparse','nrows',3}, '', 'AntennaPosition');
end

function validateAntennaHeight(ht)

validateattributes(ht, {'numeric'}, ...
    {'nonnegative','real','finite','nonsparse','row'}, '', 'AntennaHeight');
end

function validateSystemLoss(sloss)

validateattributes(sloss, {'numeric'}, ...
    {'nonnegative','real','finite','nonsparse','row'}, '', 'SystemLoss');
end