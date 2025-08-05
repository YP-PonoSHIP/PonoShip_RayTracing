classdef propagationData < matlab.mixin.internal.Scalar
%
    
% Copyright 2019-2024 The MathWorks, Inc.
    
    properties
        Name {validateName} = 'Propagation Data'
    end
    
    properties (SetAccess = private)
        Data
    end
    
    properties(Dependent)
        DataVariableName
    end
    
    properties(Hidden)
        RemoveDuplicatesFcn(1,1) function_handle = @mean
        EnableAltitude(1,1) logical = false
        AltitudeReference(1,:) char {mustBeAltitudeReference(AltitudeReference)} = 'geoid'
    end
    
    properties(Access = private)
        pDataVariableNames
        pDataVariableName
        pLatitudes
        pLongitudes
        pAltitudes
        pTransmitterLatitude
        pTransmitterLongitude
        pTransmitterMaxRange
        pImageSize
        pMATLABVersion
    end
    
    % Transient instance-specific properties (do not load)
    properties(Transient, Access = private)
        pUID
        plotID
        contourID
    end

    properties(Constant, Hidden)
        LatitudeVariableName = 'latitude'
        LongitudeVariableName = 'longitude'
        AltitudeVariableName = 'altitude'
        AltitudeReferenceChoices = {'geoid','ellipsoid','ground','surface'}
        NumRowsDataDisp = 5
    end
    
    methods
        function pd = set.Name(pd, name)
            if isstring(name)
                name = char(name);
            end
            pd.Name = name;
        end
        
        function pd = set.DataVariableName(pd, variableName)
            pd.pDataVariableName = validateDataVariableName(pd, variableName, 'propagationData');
        end
        
        function value = get.DataVariableName(pd)
            value = pd.pDataVariableName;
        end
        
        function pd = propagationData(varargin)
            
            narginchk(1,inf)
            
            % Parse inputs and get data table
            p = inputParser;
            p.addParameter('Name', message('shared_channel:rfprop:PropagationDataTitle').getString);
            p.addParameter('RemoveDuplicatesFcn', @mean);
            p.addParameter('DataVariableName', '');
            p.addParameter('EnableAltitude', false);
            p.addParameter('AltitudeReference', 'geoid');
            if isnumeric(varargin{1})
                narginchk(4,inf)
                dataTable = parseCustomSyntax(p,varargin);
            elseif istable(varargin{1})
                dataTable = parseTableSyntax(p,varargin);
            else
                dataTable = parseFileSyntax(p,varargin);
            end
            
            % Get Latitude and Longitude variables
            allVariableNames = dataTable.Properties.VariableNames;
            [isLatitude, isLongitude] = validateLatitudeLongitudeVariableNames(allVariableNames);
            lats = dataTable{:,isLatitude};
            lons = dataTable{:,isLongitude};
            
            % Get optional Altitude variable
            isAltitude = cellfun(isVariableMatchFcn(pd.AltitudeVariableName),allVariableNames);
            hasCustomAltitude = p.Results.EnableAltitude && (nnz(isAltitude) == 1);
            if hasCustomAltitude
                alts = dataTable{:,isAltitude};
                altsref = p.Results.AltitudeReference;
            else
                % If no Altitude data found, set to correspond to surface
                alts = zeros(numel(lats),1);
                altsref = 'ground';
            end
            
            % Remove missing location rows
            if isnumeric(lats) && isnumeric(lons) && isnumeric(alts)
                locs = [lats lons alts];
                if any(ismissing(locs),'all')
                    [locs,missingInd] = rmmissing(locs);
                    lats = locs(:,1);
                    lons = locs(:,2);
                    alts = locs(:,3);
                    dataTable = dataTable(~missingInd,:);
                    warning(message('shared_channel:rfprop:PropagationDataMissingLocation'));
                end
                
                pd.pLatitudes = lats;
                pd.pLongitudes = wrapTo180(lons);
                pd.pAltitudes = alts;
            end
            validateLatitude(lats);
            validateLongitude(lons,numel(lats));
            
            % Get numeric data variables
            if hasCustomAltitude
                allVariableNames = allVariableNames(~isLatitude & ~isLongitude & ~isAltitude);
            else
                allVariableNames = allVariableNames(~isLatitude & ~isLongitude);
            end
            variableNames = {};
            for nameInd = 1:numel(allVariableNames)
                variableName = allVariableNames{nameInd};
                if isValidDataVariable(dataTable{:,variableName})
                    variableNames{end+1} = variableName; %#ok<AGROW>
                end 
            end
            pd.pDataVariableNames = variableNames;
            
            % Validate that at least one valid data variable exists
            if isempty(variableNames)
                error(message('shared_channel:rfprop:PropagationDataNoValidData'))
            end

            % Set property values
            inputs = p.Results;
            pd.Name = inputs.Name;
            pd.RemoveDuplicatesFcn = inputs.RemoveDuplicatesFcn;
            pd.Data = dataTable;
            if ismember('DataVariableName',p.UsingDefaults)
                dataVariableName = variableNames{1};
            else
                dataVariableName = inputs.DataVariableName;
            end
            pd.DataVariableName = dataVariableName;
            pd.EnableAltitude = inputs.EnableAltitude;
            pd.AltitudeReference = altsref;
            
            % Set prop data ID and MATLAB version
            mlock % Lock persistent variables
            pd.pUID = sprintf('propdata%u', createUIDNumber); 
            pd.plotID = ['plot' pd.pUID];
            pd.contourID = ['contour' pd.pUID];
            pd.pMATLABVersion = version('-release');
        end
        
        function disp(pd)
            %disp   Display propagation data object
            
            % Display default object display
            builtin('disp',pd);    
            
            % Display top rows of data table
            if isscalar(pd)
                dataTable = pd.Data;
                numRows = min(pd.NumRowsDataDisp,size(dataTable,1));
                dataTableTop = head(dataTable,numRows);
                dataTableHeaderID = 'shared_channel:rfprop:PropagationDataDiplayDataHeader';
                dataTableHeader = message(dataTableHeaderID,sprintf('%u',numRows)).getString;
                fprintf('  %s\n\n', dataTableHeader);
                disp(dataTableTop);
            end
        end
                
        function V = interp(pd, latq, lonq, varargin)
            %interp   Interpolate RF propagation data
            %   V = interp(PD,LAT,LON) interpolates the propagation data object PD,
            %   returning a value in V for each of the query points in vectors LAT and
            %   LON. The interpolation is performed using scattered data interpolation.
            %   Values corresponding to query points outside of the data region are
            %   assigned NaN.
            %
            %   V = interp(___,Name,Value) interpolates the propagation data with
            %   additional options specified by one or more Name-Value pairs.
            %
            %   interp Name-Value pairs:
            %
            %   DataVariableName - Data variable to interpolate, specified as a
            %      character vector or string scalar corresponding to a variable name
            %      in the Data table of PD. The default value is the DataVariableName
            %      property value of PD.
            %
            %   Method - Method used to interpolate the data, specified as one of the
            %      following options:
            %         'linear'  - Linear interpolation
            %         'nearest' - Nearest neighbor interpolation
            %         'natural' - Natural neighbor interpolation
            %      The default value is 'natural'.
            %
            %   % Example: Plot area in service of all transmitter sites
            %
            %   % Define names and locations of sites around Boston
            %   names = ["Fenway Park","Faneuil Hall","Bunker Hill Monument"];
            %   lats = [42.3467,42.3598,42.3763];
            %   lons = [-71.0972,-71.0545,-71.0611];
            %
            %   % Create transmitter site array
            %   txs = txsite("Name", names,...
            %      "Latitude",lats,...
            %      "Longitude",lons, ...
            %      "TransmitterFrequency",2.5e9);
            %
            %   % Compute received power data for each transmitter site
            %   maxr = 20000;
            %   pd1 = coverage(txs(1),"MaxRange",maxr);
            %   pd2 = coverage(txs(2),"MaxRange",maxr);
            %   pd3 = coverage(txs(3),"MaxRange",maxr);
            %
            %   % Compute rectangle containing locations of all data
            %   locs = [location(pd1); location(pd2); location(pd3)];
            %   [minlatlon, maxlatlon] = bounds(locs);
            %
            %   % Create grid of locations over rectangle
            %   gridlength = 300;
            %   latv = linspace(minlatlon(1),maxlatlon(1),gridlength);
            %   lonv = linspace(minlatlon(2),maxlatlon(2),gridlength);
            %   [lons,lats] = meshgrid(lonv,latv);
            %   lats = lats(:);
            %   lons = lons(:);
            %
            %   % Get data for each transmitter at grid locations using interpolation
            %   v1 = interp(pd1,lats,lons);
            %   v2 = interp(pd2,lats,lons);
            %   v3 = interp(pd3,lats,lons);
            %
            %   % Create propagation data containing minimum received power values
            %   minReceivedPower = min([v1 v2 v3],[],2,"includenan");
            %   pd = propagationData(lats,lons,"MinReceivedPower",minReceivedPower);
            %
            %   % Plot minimum received power, which shows weakest signal received
            %   % from any transmitter site. The area shown may correspond to the
            %   % service area of triangulation using the three transmitter sites.
            %   sensitivity = -110;
            %   contour(pd,"Levels",sensitivity:-5,"Type","power")
            %
            % See also plot, contour, location, getDataVariable
            
            p = inputParser;
            p.addParameter('DataVariableName',pd.DataVariableName);
            p.addParameter('Method','natural');
            p.parse(varargin{:});
            
            % Validate and get parameters
            dataVariableName = p.Results.DataVariableName;
            [data, lats, lons] = pd.getDataVariable(dataVariableName);
            interpolationMethod = validatestring(p.Results.Method, ...
                {'linear','nearest','natural'}, 'interp', 'Method');
            validateLatitude(latq);
            validateLongitude(lonq,numel(latq));
            
            % Create interpolant and query. Catch known reasons for
            % interpolation failure and throw as errors.
            oldTooFewWarnState = warning('error','MATLAB:scatteredInterpolant:TooFewPtsInterpWarnId');
            tooFewWarnCleanup = onCleanup(@() warning(oldTooFewWarnState));
            oldEmptyWarnState = warning('error','MATLAB:scatteredInterpolant:InterpEmptyTri2DWarnId');
            emptyWarnCleanup = onCleanup(@() warning(oldEmptyWarnState));
            F = scatteredInterpolant(double(lats), double(lons), double(data), ...
                interpolationMethod, 'none');
            try
                V = F(double(latq(:)), wrapTo180(double(lonq(:))));
            catch interpErr
                err = MException(message('shared_channel:rfprop:PropagationDataInterpolationError'));
                err = err.addCause(interpErr);
                throw(err)
            end
        end
        
        function varargout = location(pd)
            %location   Coordinates of RF propagation data
            %   L = location(PD) returns the location coordinates L of the data points
            %   in propagation data object PD as an M-by-2 vector [LAT LON] of degrees
            %   latitude LAT and longitude LON. The number of rows M corresponds to the
            %   number of rows in the Data table with valid latitude and longitude
            %   values. Duplicate locations are not removed. The output LON is wrapped
            %   so that values are in the range [-180 180].
            %
            %   [LAT,LON] = location(PD) returns the latitude LAT and longitude LON
            %   of propagation data object PD.
            %
            %   % Example: Plot area in service of all transmitter sites
            %
            %   % Define names and locations of sites around Boston
            %   names = ["Fenway Park","Faneuil Hall","Bunker Hill Monument"];
            %   lats = [42.3467,42.3598,42.3763];
            %   lons = [-71.0972,-71.0545,-71.0611];
            %
            %   % Create transmitter site array
            %   txs = txsite("Name", names,...
            %      "Latitude",lats,...
            %      "Longitude",lons, ...
            %      "TransmitterFrequency",2.5e9);
            %
            %   % Compute received power data for each transmitter site
            %   maxr = 20000;
            %   pd1 = coverage(txs(1),"MaxRange",maxr);
            %   pd2 = coverage(txs(2),"MaxRange",maxr);
            %   pd3 = coverage(txs(3),"MaxRange",maxr);
            %
            %   % Compute rectangle containing locations of all data
            %   locs = [location(pd1); location(pd2); location(pd3)];
            %   [minlatlon, maxlatlon] = bounds(locs);
            %
            %   % Create grid of locations over rectangle
            %   gridlength = 300;
            %   latv = linspace(minlatlon(1),maxlatlon(1),gridlength);
            %   lonv = linspace(minlatlon(2),maxlatlon(2),gridlength);
            %   [lons,lats] = meshgrid(lonv,latv);
            %   lats = lats(:);
            %   lons = lons(:);
            %
            %   % Get data for each transmitter at grid locations using interpolation
            %   v1 = interp(pd1,lats,lons);
            %   v2 = interp(pd2,lats,lons);
            %   v3 = interp(pd3,lats,lons);
            %
            %   % Create propagation data containing minimum received power values
            %   minReceivedPower = min([v1 v2 v3],[],2,"includenan");
            %   pd = propagationData(lats,lons,"MinReceivedPower",minReceivedPower);
            %
            %   % Plot minimum received power, which shows weakest signal received
            %   % from any transmitter site. The area shown may correspond to the
            %   % service area of triangulation using the three transmitter sites.
            %   sensitivity = -110;
            %   contour(pd,"Levels",sensitivity:-5,"Type","power")
            %
            % See also getDataVariable, interp
            
            lats = pd.pLatitudes;
            lons = pd.pLongitudes;
            alts = pd.pAltitudes;
            
            % Assign outputs
            varargout = cell(1,nargout);
            if nargout <= 1 % Can be 0 to assign into ans
                varargout{1} = [lats lons];
            else
                varargout{1} = lats;
                varargout{2} = lons;
                if pd.EnableAltitude
                    varargout{3} = alts;
                end
            end
        end
        
        function [v, lats, lons, alts] = getDataVariable(pd, variableName)
            %getDataVariable   Get data variable values
            %   V = getDataVariable(PD) returns the values of the data points in
            %   propagation data object PD corresponding to the DataVariableName
            %   property. The output is an M-by-1 double vector, where M is the number
            %   of data points used for plotting. The data is processed so that missing
            %   values are removed and duplicate location data are replaced with their
            %   mean value.
            %
            %   [V,LAT,LON] = getDataVariable(PD) additionally returns location
            %   coordinates LAT and LON of the data points in propagation data object
            %   PD corresponding to the DataVariableName property. The outputs are
            %   M-by-1 double vectors, where M is the number of data points used for
            %   plotting. The outputs LAT and LON are specified in degrees, and LON is
            %   wrapped so that values are in the range [-180 180].
            %
            %   [___] = getDataVariable(PD,VARNAME) returns the values of the data
            %   points corresponding to the VARNAME variable, which is specified as a
            %   character vector or string scalar corresponding to a variable name in
            %   the Data table. The variable name VARNAME must correspond to a variable
            %   with numeric data and cannot correspond to latitude or longitude.
            %
            % See also location, interp
            
            % Get variable name from input or property
            if nargin < 2
                variableName = pd.DataVariableName;
            else
                variableName = validateDataVariableName(pd, variableName, 'getDataVariable');
            end

            % Get variable values and locations
            v = pd.Data{:,variableName};
            lats = pd.pLatitudes;
            lons = pd.pLongitudes;
            alts = pd.pAltitudes;
            
            % Process data to remove missing and duplicate values
            [v,lats,lons,alts] = pd.removeMissingData(v,lats,lons,alts);
            [v,lats,lons,alts] = pd.removeDuplicateData(v,lats,lons,alts,pd.RemoveDuplicatesFcn);
        end
    end
    
    methods(Hidden)
        function pd = setContourProperties(pd,txlat,txlon,txmaxrange,imageSize,uid)
            pd.pTransmitterLatitude = txlat;
            pd.pTransmitterLongitude = txlon;
            pd.pTransmitterMaxRange = txmaxrange;
            pd.pImageSize = imageSize;
            
            % An ID is passed when association is required with a plot
            % created for a site object
            if nargin > 5
                pd.contourID = uid;
            end
        end
        
        function varname = validateDataVariableName(pd, variableName, fcnName)
            try
                varname = validatestring(variableName, pd.pDataVariableNames, ...
                    fcnName, 'DataVariableName');
            catch e
                if (ischar(variableName) || isStringScalar(variableName)) && ...
                        ismember(lower(variableName), lower(pd.Data.Properties.VariableNames))
                    error(message('shared_channel:rfprop:PropagationDataInvalidVariable',variableName))
                else
                    rethrow(e)
                end
            end
        end
    end
    
    methods(Static, Hidden)
        function pd = loadobj(s)
            %loadobj  Load propagationData object
            
            if isstruct(s)
                pd = propagationData(s.Data, ...
                    'Name', s.Name, ...                   
                    'DataVariableName', s.DataVariableName, ...
                    'RemoveDuplicatesFcn', s.RemoveDuplicatesFcn, ...
                    'EnableAltitude', s.EnableAltitude, ...
                    'AltitudeReference',s.AltitudeReference);
            else
                % Assign loaded object and generate new unique ID
                pd = s;
                pd.pUID = sprintf('propdata%u', createUIDNumber);
                pd.plotID = ['plot' pd.pUID];
                pd.contourID = ['contour' pd.pUID];
            end
        end

        function names = varnamechoices(pd)
            names = pd.pDataVariableNames;
        end
        
        function [v,lats,lons,alts] = removeMissingData(v,lats,lons,alts)
            
            if any(ismissing(v),'all')
                vlocs = rmmissing([v lats lons alts]);
                v = vlocs(:,1);
                lats = vlocs(:,2);
                lons = vlocs(:,3);
                alts = vlocs(:,4);
            end
        end
        
        function [uniqueData,lats,lons,alts] = removeDuplicateData(v,lats,lons,alts,removeDuplicatesFcn)
            
            % Get unique locations
            [locs,~,ic] = unique([lats lons alts],'rows','stable');            
            
            % Return early if no duplicate data
            if size(locs,1) == numel(lats)
                uniqueData = v;
                return
            end
            
            % Find all data associated with each location and replace with
            % single value using function
            numLoc = size(locs,1);
            uniqueData = zeros(numLoc,1);
            for locInd = 1:numLoc
                locData = v(ic == locInd);
                if numel(locData) > 1
                    locData = removeDuplicatesFcn(locData);
                end
                uniqueData(locInd) = locData;
            end
            lats = locs(:,1);
            lons = locs(:,2);
            alts = locs(:,3);
        end
        
        function dataOut = discretizeDataToLevels(dataIn, levels)
            %discretizeDataToLevels   Discretize data using levels
            
            datalevels = sort(levels);
            maxBin = max(max(datalevels(:)),max(dataIn(:))) + 1; % Need max bin edge to include all values
            bins = [datalevels(:); maxBin];
            dataOut = discretize(dataIn,bins,datalevels);
        end
    end
    
    methods
        contour(pd, varargin)
        plot(pd, varargin)
    end
end

function fcn = isVariableMatchFcn(varName)
fcn = @(x)(numel(x) >=3) && startsWith(varName,lower(x));
end

function propDataTable = parseCustomSyntax(p,inputs)

try
    % The 3rd input must be either a dataVariableName, name-value pair, or altitude
    validAltitudeOrVariableName(inputs{3});
    p.addRequired('Latitude');
    p.addRequired('Longitude');
    p.addOptional('Altitude',[],@(x)isnumeric(x) && isvector(x));
    p.PartialMatching = false; % Match full property names only (case-insensitive)
    p.KeepUnmatched = true;
    p.parse(inputs{:});
    
    % Get lat/lon/alt and verify sizes match
    lat = p.Results.Latitude;
    lon = p.Results.Longitude;
    alt = p.Results.Altitude;
    numLocations = numel(lat);
    validateattributes(lon, {'numeric'}, ...
        {'numel',numLocations}, 'propagationData', 'Longitude');
    hasAltitude = p.Results.EnableAltitude && ~isempty(alt) && isnumeric(alt);
    if hasAltitude
        validateAltitude(alt,numLocations);
    end
    
    % Initialize variable names/values for data table
    if hasAltitude
        tableVariableNames = {'Latitude' 'Longitude' 'Altitude'};
        tableVariableValues = {lat(:) lon(:) alt(:)};
        paramNames = p.Parameters;
    else
        tableVariableNames = {'Latitude' 'Longitude'};
        tableVariableValues = {lat(:) lon(:)};
        paramNames = setdiff(p.Parameters,'Altitude');
    end
    
    % Build variable names/values from inputs
    inputVariables = inputs(numel(tableVariableNames)+1:end);    
    for varInd = 1:2:numel(inputVariables)
        variableName = inputVariables{varInd};
        if ~matlab.internal.datatypes.isScalarText(variableName)
            error(message('shared_channel:rfprop:PropagationDataInvalidVariableName'));
        elseif ismember(variableName, tableVariableNames)
            error(message('shared_channel:rfprop:PropagationDataDuplicateVariableName'));
        elseif strcmpi(variableName,'Data')
            error(message('shared_channel:rfprop:PropagationDataReadOnlyDataProperty'));
        elseif ismember(variableName, paramNames)
            continue
        end
        
        variableValue = inputVariables{varInd+1};
        if ~isValidDataVariable(variableValue)
            error(message('shared_channel:rfprop:PropagationDataInvalidVariable',variableName))
        elseif ~isvector(variableValue) || ~isequal(numel(variableValue),numLocations)
            error(message('shared_channel:rfprop:PropagationDataInvalidVariableSize',variableName))
        else
            tableVariableNames = [tableVariableNames variableName]; %#ok<AGROW>
            tableVariableValues = [tableVariableValues variableValue(:)]; %#ok<AGROW>
        end
    end
    
    % Create table from data variables
    propDataTable = table(tableVariableValues{:},'VariableNames',tableVariableNames);
catch e
    throwAsCaller(e);
end
end

function propDataTable = parseTableSyntax(p,inputs)

try
    p.addRequired('Table');
    p.parse(inputs{:});
    propDataTable = p.Results.Table;
catch e
    throwAsCaller(e);
end
end

function propDataTable = parseFileSyntax(p,inputs)

% Parse inputs and get file name
try
    p.addRequired('FileName');
    p.parse(inputs{:});
    fileName = p.Results.FileName;
catch e
    throwAsCaller(e);
end

% Get readtable import options
try
    opts = detectImportOptions(fileName,'PreserveVariableNames',true);
catch readErr
    % Rethrow common "file not found" error but otherwise throw custom
    % error suggesting use of readtable to customize import.
    if strcmpi(readErr.identifier, 'MATLAB:textio:textio:FileNotFound')
        throwAsCaller(readErr)
    else
        err = MException(message('shared_channel:rfprop:PropagationDataReadTableError',fileName));
        err = err.addCause(readErr);
        throwAsCaller(err)
    end
end

% Validate variable names for latitude and longitude and read table
try
    [~,~,pass] = validateLatitudeLongitudeVariableNames(opts.VariableNames,false);
    if pass
        propDataTable = readtable(fileName,'PreserveVariableNames',true);
    else
        % Try again assuming first row contains variable names
        opts = detectImportOptions(fileName,'PreserveVariableNames',true,'NumHeaderLines',0);
        validateLatitudeLongitudeVariableNames(opts.VariableNames);
        propDataTable = readtable(fileName,'PreserveVariableNames',true,"NumHeaderLines",0);
    end
catch e
    throwAsCaller(e);
end
end

function [isLat, isLon, pass] = validateLatitudeLongitudeVariableNames(allVariableNames, throwError)

if nargin < 2
    throwError = true;
end

latVarName = propagationData.LatitudeVariableName;
isLat = cellfun(isVariableMatchFcn(latVarName),allVariableNames);
lonVarName = propagationData.LongitudeVariableName;
isLon = cellfun(isVariableMatchFcn(lonVarName),allVariableNames);
hasOneLat = (nnz(isLat) == 1);
hasOneLon = (nnz(isLon) == 1);
pass = hasOneLat && hasOneLon;

if throwError && ~pass
    if ~hasOneLat
        error(message('shared_channel:rfprop:PropagationDataInvalidLatitudeVariable'));
    end
    if ~hasOneLon
        error(message('shared_channel:rfprop:PropagationDataInvalidLongitudeVariable'));
    end
end
end

function validateName(name)
% Use scalartext attribute, which does not allow string array
% or cell array validation
validateattributes(name, {'char','string'}, {'scalartext'}, '', 'Name');
end

function validateLatitude(lat)

% Note that >= and <= imply real, positive, finite
validateattributes(lat, {'numeric'}, ...
    {'real','finite','nonsparse','vector','>=',-90,'<=',90}, 'propagationData', 'Latitude');
end

function validateLongitude(lon, numloc)

% While Latitude is constrained to range [-90,90], Longitude is
% unconstrained and allowed to wrap
validateattributes(lon, {'numeric'}, ...
    {'real','finite','nonsparse','vector','numel',numloc}, 'propagationData', 'Longitude');
end

function validateAltitude(alt,numloc)

validateattributes(alt, {'numeric'}, ...
    {'real','finite','nonsparse','vector','numel',numloc}, 'propagationData', 'Altitude');
end

function validAltitudeOrVariableName(param)
% Throws an error if the input is not a valid altitude or valid
% dataVariableName/Name-value pair
validAltitude = @(x)isnumeric(x) && isvector(x);
if ~validAltitude(param) && ~isvarname(param)
    if matlab.internal.datatypes.isScalarText(param)
        error(message('MATLAB:InputParser:UnmatchedNotAValidFieldName', param));
    else
        error(message('MATLAB:InputParser:ParamMustBeChar'));
    end
end
end

function mustBeAltitudeReference(param)
    mustBeMember(param, propagationData.AltitudeReferenceChoices);
end

function isvar = isValidDataVariable(data)

isvar = isnumeric(data) && ~issparse(data) && isreal(data) && any(~isnan(data(:)));
end

function UID = createUIDNumber
%createUIDNumber   Create unique identifier number

persistent uidCounter;
if isempty(uidCounter)
    uidCounter = 0;
end

% Set counter value to input or else increment if 'next'
uidCounter = uidCounter + 1;
UID = uidCounter;
end