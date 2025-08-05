function r = sinr(rxs, txs, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate site inputs
validateattributes(rxs,{'rxsite'},{'nonempty'},'sinr','',1);
validateattributes(txs,{'txsite'},{'nonempty'},'sinr','',2);

% Process optional inputs
p = inputParser;
if nargin > 2 && mod(numel(varargin),2)
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('SignalSource', 'strongest');
p.addParameter('ReceiverNoisePower', -107);
p.addParameter('ReceiverGain', []);
p.addParameter('ReceiverAntennaHeight', []);
p.addParameter('Map', []);
p.addParameter('TransmitterAntennaSiteCoordinates', []);

% Add plot-only parameters which are silently ignored
p.addParameter('MaxRange', []);
p.addParameter('Values', []);
p.addParameter('Resolution', []);
p.addParameter('Colormap', []);
p.addParameter('ColorLimits', []);
p.addParameter('Transparency', []);
p.addParameter('ShowLegend', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.TransmitterAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(rxs, txs);
else
    usingCartesian = strcmp(p.Results.TransmitterAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

% Validate and get parameters
numRxs = numel(rxs);
if usingCartesian
    [map, ~, ~, mapStruct] = rfprop.internal.Validators.validateCartesianMap(p);
    pm = rfprop.internal.Validators.validateCartesianPropagationModel(p, map, 'sinr');
else
    [map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'sinr');
    pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'sinr');
end
rfprop.internal.Validators.validateMaxNumReflections(pm, 'sinr');
rfprop.internal.Validators.validateMaxNumDiffractions(pm, 'sinr'); % Disable 2 order diffractions
sigSource = validateSignalSource(p, numRxs);
noisePower = validateReceiverNoisePower(p);
[rxGain, usingDefaultGain] = validateReceiverGain(p);
[rxHeight, usingDefaultHeight] = validateReceiverAntennaHeight(p);

% Create vector array of all txs
txs = txs(:);
usingSiteSigSource = isa(sigSource,'txsite');
if usingSiteSigSource
    % Include signal sources in txsite list
    txs = union(txs,sigSource,'stable');
end

txsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TransmitterAntennaSiteCoordinates, txs, map, 'sinr');

% Compute SignalStrength values if not passed in
if isempty(txsCoords.CustomData) || ~isfield(txsCoords.CustomData,'SignalStrength')
    % If specified, temporarily override AntennaHeight for rxs. The height is
    % used in distance/angle calculations as well as by the propagation model.
    if ~usingDefaultHeight
        originalRxHeights = [rxs.AntennaHeight];
        setAntennaHeights(rxs,rxHeight);
        restoreRxAntennaHeights = onCleanup(@()setAntennaHeights(rxs,originalRxHeights));
    end
    
    % Calculate signal strengths (in dBm) at receivers due to transmitters
    args = {'Type', 'power', ...
        'Map', mapStruct, ...
        'TransmitterAntennaSiteCoordinates', txsCoords};
    if ~usingDefaultGain
        args = [args, 'ReceiverGain', rxGain];
    end
    ss = sigstrength(rxs, txs, pm, args{:});
else
    % Get signal strength data from input and derive number of receiver
    % sites corresponding to data
    ss = txsCoords.CustomData.SignalStrength;
    numRxs = size(ss,2);
end

% Assign txs signal source index for each receiver
if usingSiteSigSource    
    sigSourceIndices = zeros(numRxs,1);
    if isscalar(sigSource)
        % Signal source is the same for all receivers
        sigSourceIndices(:) = find(txs == sigSource);
    else
        % Signal source is defined for each receiver
        for rxInd = 1:numRxs
            sigSourceIndices(rxInd) = find(txs == sigSource(rxInd));
        end
    end
end
txFqs = [txs.TransmitterFrequency];

% Convert values from dBm to W to allow simple arithmetic
ss = dbm2pow(ss);
noisePower = dbm2pow(noisePower);

% Calculate SINR for each receiver
r = zeros(1,numRxs);
for rxInd = 1:numRxs
    
    % Get signal strengths from transmitters to this receiver
    sigStrengths = ss(:,rxInd)';    

    % Compute signal
    if strcmp(sigSource, 'strongest')
        % Use signal source that has greatest strength
        [sigSourcePower, sigSourceInd] = max(sigStrengths);
    else
        sigSourceInd = sigSourceIndices(rxInd);
        sigSourcePower = sigStrengths(sigSourceInd);
    end
        
    % Get interference frequencies and powers (remove signal source from txs)
    intSourceFqs = txFqs;
    intSourcePowers = sigStrengths;
    intSourceFqs(sigSourceInd) = [];
    intSourcePowers(sigSourceInd) = [];
    
    % Compute interference power from sources with matching frequency
    sigSourceFq = txFqs(sigSourceInd);
    intSourceInd = (intSourceFqs == sigSourceFq);
    intSourcePowers = intSourcePowers(intSourceInd);
    interferencePower = sum(intSourcePowers);
    
    % Compute SINR in dB
    r0 = sigSourcePower / (interferencePower + noisePower);
    rdb = 10*log10(r0); % Convert to dB
    
    % Assign output
    r(rxInd) = rdb;
end

% Restore original antenna heights
if ~usingDefaultHeight
    delete(restoreRxAntennaHeights);
end
end

function setAntennaHeights(sites, hts)

applyScalarValue = isscalar(hts);
for rxInd = 1:numel(sites)
    if applyScalarValue
        sites(rxInd).AntennaHeight = hts;
    else
        sites(rxInd).AntennaHeight = hts(rxInd);
    end
end
end

function sigsource = validateSignalSource(p, numRxs)

try
    sigsource = p.Results.SignalSource;
    if ischar(sigsource) || isstring(sigsource)
        sigsource = validatestring(sigsource, {'strongest'}, ...
            'sinr','SignalSource');
    elseif isscalar(sigsource)
        validateattributes(sigsource,{'txsite'}, {'scalar'}, ...
            'sinr','SignalSource');
    else
        validateattributes(sigsource,{'txsite'}, {'numel',numRxs}, ...
            'sinr','SignalSource');
        sigsource = sigsource(:); % Guarantee column vector
    end
catch e
    throwAsCaller(e);
end
end

function noisePower =  validateReceiverNoisePower(p)

try
    noisePower = p.Results.ReceiverNoisePower;
    validateattributes(noisePower, {'numeric'}, {'real','finite','nonnan','nonsparse','scalar'}, ...
        'sinr', 'ReceiverNoisePower');
catch e
    throwAsCaller(e);
end
end

function [rxGain, usingDefaultGain] = validateReceiverGain(p)

try
    rxGain = p.Results.ReceiverGain;
    usingDefaultGain = ismember('ReceiverGain',p.UsingDefaults);
    if ~usingDefaultGain
        validateattributes(rxGain,{'numeric'}, {'real','finite','nonnan','nonsparse','scalar'}, ...
            'sinr', 'ReceiverGain');
    end
catch e
    throwAsCaller(e);
end
end

function [rxHeight, usingDefaultHeight] = validateReceiverAntennaHeight(p)

try
    rxHeight = p.Results.ReceiverAntennaHeight;
    usingDefaultHeight = ismember('ReceiverAntennaHeight',p.UsingDefaults);
    if ~usingDefaultHeight
        validateattributes(rxHeight,{'numeric'}, {'real','finite','nonnan','nonsparse','scalar','nonnegative', ...
            '<=',rfprop.Constants.MaxPropagationDistance}, 'sinr', 'ReceiverAntennaHeight');
    end
catch e
    throwAsCaller(e);
end
end

function w = dbm2pow(dbm)
w = 10.^((dbm-30)./10);
end