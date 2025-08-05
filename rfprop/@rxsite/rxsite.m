classdef (Sealed) rxsite < rfprop.AntennaSite
%

% Copyright 2017-2024 The MathWorks, Inc.

    properties(Dependent)
        ReceiverSensitivity(1,1) {validateSensitivity} % dBm
    end
    
    properties(Access = private)
        pReceiverSensitivity = -100
    end
    
    properties(Hidden, Constant)
        DefaultMarkerIcon = '/toolbox/shared/globe/globeviewer/release/globeviewer/images/receiver.svg'
    end
    
    methods
        function set.ReceiverSensitivity(rx, sens)
            rx.pReceiverSensitivity = sens;
        end
        
        function sens = get.ReceiverSensitivity(rx)
            sens = rx.pReceiverSensitivity;
        end
    end

    methods
        function rxs = rxsite(varargin)
            
            % Process Name-Value parameters. Note that some superclass
            % properties may have rx-specific default values.
            p = inputParser;
            if mod(numel(varargin),2)
                % Validator function is necessary for inputParser to allow string
                % option instead of treating it like parameter name
                p.addOptional('CoordinateSystem', 'geographic', ...
                    @(x)matlab.internal.datatypes.isScalarText(x));
            else
                p.addParameter('CoordinateSystem', 'geographic');
            end
            p.addParameter('Name', '');
            p.addParameter('Latitude', 42.3021); % Lakeside
            p.addParameter('Longitude', -71.3764);
            p.addParameter('Antenna', 'isotropic');
            p.addParameter('AntennaPosition', [0;0;0]);
            p.addParameter('AntennaAngle', 0);
            p.addParameter('AntennaHeight', 1);
            p.addParameter('SystemLoss', 0);
            p.addParameter('ReceiverSensitivity', -100);
            inputs = rfprop.AntennaSite.parseInputs(p, varargin{:});
            
            % Get inputs
            usingDefaults = inputs.UsingDefaults;
            sens = inputs.ReceiverSensitivity;
            
            % Allocate sites array by assigning into last rx. Compute size
            % based on longest length property value. All other values must
            % either match length or be single site value.
            numSites = max([rfprop.AntennaSite.numSites(inputs), ...
                size(sens,2)]);
            siteNum = rxs.siteNumber('defaultname');
            rxs(numSites).Name = '';
            
            % Reset counter since MATLAB creates a default object to create
            % the site array, which increments counter as side effect
            rxs.siteNumber('defaultname', siteNum);
            
            % Validate all property values. For performance, only validate
            % non-default values.
            inputs = rxs.validateSiteProperties(inputs);
            if ~ismember('ReceiverSensitivity', usingDefaults)
                rxs.validateNumColumns(sens, 'ReceiverSensitivity')
                validateSensitivity(sens)
            end
            
            % Set properties on each site. For performance, array values
            % are validated once (above) and set into private fields,
            % thereby preventing repeated validation for each site. All
            % values need to be set since default values are defined in
            % inputParser.
            rxs.setSiteProperties(inputs);
            isScalarSens = isscalar(sens);
            for k = 1:numel(rxs)
                if isScalarSens
                    rxs(k).pReceiverSensitivity = sens;
                else
                    rxs(k).pReceiverSensitivity = sens(k);
                end
            end
            
            % Set names/IDs last so that counters are not incremented if
            % any error occurs above
            rxs.setSiteIdentifiers(inputs);
        end
    end
    
    methods
        ss = sigstrength(rxs, txs, varargin)
        varargout = link(rxs, txs, varargin)
        r = sinr(rxs, txs, varargin)
    end
    
    methods(Access = protected)
        function props = getSiteSpecificProperties(~)
            props = {'ReceiverSensitivity'};
        end
    end
    
    methods(Static, Hidden)
        function rx = loadobj(s)
            %loadobj  Load rxsite object
            
            % Automatically convert empty to dipole antenna. An empty
            % Antenna may appear if loading an object from a previous
            % release where dipole was default but only constructed on
            % first get.
            if isequal(s.Antenna,[])
                s.Antenna = design(dipole,1.9e9);
            end
            
            if isstruct(s)
                % Handle fields added in R2020b
                if isfield(s,'CoordinateSystem')
                    coordSys = s.CoordinateSystem;
                else
                    coordSys = 'geographic';
                end
                if isfield(s,'AntennaPosition')
                    antPos = s.AntennaPosition;
                else
                    antPos = [0;0;0];
                end
                
                rx = rxsite('Name', s.Name, ...
                    'CoordinateSystem', coordSys, ...
                    'Latitude', s.Latitude, ...
                    'Longitude', s.Longitude, ...
                    'Antenna', s.Antenna, ...
                    'AntennaPosition', antPos, ...
                    'AntennaAngle', s.AntennaAngle, ...
                    'AntennaHeight', s.AntennaHeight, ...
                    'SystemLoss', s.SystemLoss, ...
                    'ReceiverSensitivity', s.ReceiverSensitivity);
            else
                % Assign loaded object and generate new unique ID
                rx = s;
                setSiteIdentifiers(rx);
            end
        end
    end
end

function validateSensitivity(sens)

try
    validateattributes(sens, {'numeric'}, ...
        {'real','finite','nonsparse','row'}, 'rxsite', 'ReceiverSensitivity');
catch e
    throwAsCaller(e);
end
end