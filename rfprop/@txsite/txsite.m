classdef (Sealed) txsite < rfprop.AntennaSite
%
    
% Copyright 2017-2024 The MathWorks, Inc.
    
    properties(Dependent)
        TransmitterFrequency(1,1) % Hz
        TransmitterPower(1,1) {validatePower} % W
    end
    
    properties(Access = private)
       pTransmitterFrequency = 1.9e9
       pTransmitterPower = 10
    end
    
    properties(Hidden, Constant)
        DefaultMarkerIcon = '/toolbox/shared/globe/globeviewer/release/globeviewer/images/transmitter.svg'
    end
    
    methods
        function set.TransmitterFrequency(tx, fq)
            tx.validateTransmitterFrequency(fq)
            tx.pTransmitterFrequency = fq;
        end
        
        function fq = get.TransmitterFrequency(tx)
            fq = tx.pTransmitterFrequency;
        end
        
        function set.TransmitterPower(tx, pow)
            tx.pTransmitterPower = pow;
        end
        
        function pow = get.TransmitterPower(tx)
            pow = tx.pTransmitterPower;
        end
    end
    
    methods
        varargout = coverage(txs, varargin)
    end
    
    methods(Hidden)
        maxrange = range(tx, ss, Grx_db, type, pm)
        contourmap(txs, varargin)
        [datalats, datalons, rxs, ss] = radialReceiverLocationsLayoutData(txs, varargin)
    end
    
    methods
        function txs = txsite(varargin)
            
            % Process Name-Value parameters. Note that some superclass
            % properties may have tx-specific default values.
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
            p.addParameter('Latitude', 42.3001); % Apple Hill
            p.addParameter('Longitude', -71.3504);
            p.addParameter('Antenna', 'isotropic');
            p.addParameter('AntennaPosition', [0;0;0]);
            p.addParameter('AntennaAngle', 0);
            p.addParameter('AntennaHeight', 10);
            p.addParameter('SystemLoss', 0);
            p.addParameter('TransmitterFrequency', 1.9e9);
            p.addParameter('TransmitterPower', 10);            
            inputs = rfprop.AntennaSite.parseInputs(p, varargin{:});
            
            % Get inputs
            usingDefaults = inputs.UsingDefaults;
            fq = inputs.TransmitterFrequency;
            pow = inputs.TransmitterPower;
                        
            % Allocate sites array by assigning into last tx. Compute size
            % based on longest length property value. All other values must
            % either match length or be single site value.
            numSites = max([rfprop.AntennaSite.numSites(inputs), ...
                size(fq,2), size(pow,2)]);
            siteNum = txs.siteNumber('defaultname');
            txs(numSites).Name = '';
            
            % Reset counter since MATLAB creates a default object to create
            % the site array, which increments counter as side effect
            txs.siteNumber('defaultname', siteNum);
            
            % Validate all property values. For performance, only validate
            % non-default values.
            inputs = txs.validateSiteProperties(inputs);
            usingDefaultFrequency = ismember('TransmitterFrequency', usingDefaults);
            if ~usingDefaultFrequency
                txs.validateNumColumns(fq, 'TransmitterFrequency')
                validateattributes(fq, {'numeric'}, ...
                    {'positive','real','finite','nonsparse','row'}, ...
                    'txsite', 'TransmitterFrequency');
                % Defer range validation until antenna is assigned
            end
            if ~ismember('TransmitterPower', usingDefaults)
                txs.validateNumColumns(pow, 'TransmitterPower')
                validatePower(pow)
            end
            
            % Set properties on each site. For performance, array values
            % are validated once (above) and set into private fields,
            % thereby preventing repeated validation for each site. All
            % values need to be set since default values are defined in
            % inputParser.
            txs.setSiteProperties(inputs);
            isScalarFq = isscalar(fq);
            isScalarPow = isscalar(pow);
            for k = 1:numel(txs)
                if isScalarFq
                    txs(k).pTransmitterFrequency = fq;
                else
                    txs(k).pTransmitterFrequency = fq(k);
                end
                if isScalarPow
                    txs(k).pTransmitterPower = pow;
                else
                    txs(k).pTransmitterPower = pow(k);
                end
            end
            
            % Validate frequency value for antenna
            if ~usingDefaultFrequency
                ant = inputs.Antenna;
                if isempty(ant) || isscalar(ant)
                    if iscell(ant)
                        ant = ant{1};
                    end
                    txs.validateTransmitterFrequency(fq, ant);
                else
                    for k = 1:numel(txs)
                        if isScalarFq
                            txs.validateTransmitterFrequency(fq, txs(k).Antenna);
                        else
                            txs.validateTransmitterFrequency(fq(k), txs(k).Antenna);
                        end
                    end
                end
            end
            
            % Set names/IDs last so that counters are not incremented if
            % any error occurs above
            txs.setSiteIdentifiers(inputs);
        end
    end
    
    methods(Access = protected)
        function props = getSiteSpecificProperties(~)
            props = {'TransmitterFrequency','TransmitterPower'};
        end
        
        function validateAntennaTransmitterFrequency(tx, ant)
            try
                validateTransmitterFrequency(tx, tx.TransmitterFrequency, ant)
            catch e
                throwAsCaller(e);
            end
        end
        
        function validateTransmitterFrequency(tx, fq, ant)
            try
                if nargin < 3
                    ant = tx.Antenna;
                end
                
                % Enforce frequency range for EM antenna objects only
                if rfprop.AntennaSite.isElectromagneticAntenna(ant)
                    attributes = {'real','finite','nonsparse','>=',1e3,'<=',200e9};
                else
                    attributes = {'positive','real','finite','nonsparse'};
                end
                validateattributes(fq,{'numeric'},attributes,'txsite','TransmitterFrequency');
            catch e
                throwAsCaller(e);
            end
        end
    end
    
    methods(Static, Hidden)
        function tx = loadobj(s)
            %loadobj  Load txsite object
            
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
                
                tx = txsite('Name', s.Name, ...
                    'CoordinateSystem', coordSys, ...
                    'Latitude', s.Latitude, ...
                    'Longitude', s.Longitude, ...
                    'Antenna', s.Antenna, ...
                    'AntennaPosition', antPos, ...
                    'AntennaAngle', s.AntennaAngle, ...
                    'AntennaHeight', s.AntennaHeight, ...
                    'SystemLoss', s.SystemLoss, ...
                    'TransmitterFrequency', s.TransmitterFrequency, ...
                    'TransmitterPower', s.TransmitterPower);
            else
                % Assign loaded object and generate new unique ID
                tx = s;
                setSiteIdentifiers(tx);
            end
        end
    end
end

function validatePower(pow)

try
    validateattributes(pow, {'numeric'}, ...
        {'positive','real','finite','nonsparse','row'}, ...
        'txsite', 'TransmitterPower');
catch e
    throwAsCaller(e);
end
end