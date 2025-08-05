function [pl, phase] = raypl(ray, varargin)
%

% Copyright 2019-2024 The MathWorks, Inc.

%#codegen

% Parse inputs
if isempty(coder.target)
    p = inputParser;
    % Use hidden ReflectionMaterials property if the N-V pair not specified
    p.addParameter("ReflectionMaterials", {'concrete'});
    p.addParameter('TransmitterPolarization', "none");
    p.addParameter('ReceiverPolarization', "none");
    p.addParameter('TransmitterAxes', eye(3));
    p.addParameter('ReceiverAxes', eye(3));
    p.addParameter('ValidateInputs', true);
    p.parse(varargin{:});
    
    mtlsChar = convertStringsToChars(p.Results.ReflectionMaterials);
    txPol = convertCharsToStrings(p.Results.TransmitterPolarization);
    rxPol = convertCharsToStrings(p.Results.ReceiverPolarization);
    txAxes = p.Results.TransmitterAxes;
    rxAxes = p.Results.ReceiverAxes;
    validateInputs = p.Results.ValidateInputs;
else
    pvPairs = struct(...
        'ReflectionMaterials', uint32(0), ...
        'TransmitterPolarization', uint32(0), ...
        'ReceiverPolarization', uint32(0), ...
        'TransmitterAxes', uint32(0), ...
        'ReceiverAxes', uint32(0), ...
        'ValidateInputs', uint32(0));
    popts = struct( ...
        'CaseSensitivity', false, ...
        'StructExpand', true, ...
        'PartialMatching', true);
    pStruct = coder.internal.parseParameterInputs( ...
        pvPairs, popts, varargin{:});
    mtlsChar = convertStringsToChars( ...
        coder.internal.getParameterValue( ...
        pStruct.ReflectionMaterials, {'concrete'}, varargin{:}));
    txPol = convertCharsToStrings( ...
        coder.internal.getParameterValue( ...
        pStruct.TransmitterPolarization, "none", varargin{:}));
    rxPol = convertCharsToStrings( ...
        coder.internal.getParameterValue( ...
        pStruct.ReceiverPolarization, "none", varargin{:}));
    txAxes = coder.internal.getParameterValue( ...
        pStruct.TransmitterAxes, eye(3), varargin{:});
    rxAxes = coder.internal.getParameterValue( ...
        pStruct.ReceiverAxes, eye(3), varargin{:});
    validateInputs = coder.internal.getParameterValue( ...
        pStruct.ValidateInputs, true, varargin{:});
end

% Validate ray object
if validateInputs
    validateattributes(ray, {'comm.Ray'}, {'scalar'}, ...
        'raypl', 'ray input');
    validateConfig(ray);
    coder.internal.errorIf( ...
        strcmp(ray.PathSpecification, 'Delay and angles'), ...
        'shared_channel:raypl:NotInLocSpec');
    coder.internal.errorIf( ...
        sum('R' == [ray.Interactions.Type], 2) ~= size(ray.Interactions, 2), ...
        'shared_channel:raypl:NotAllReflection');    
end

% Get carrier frequency
fc = ray.Frequency;

% Validate material input for NLOS
if ischar(mtlsChar)
    mtls = {mtlsChar};
else
    mtls = mtlsChar;
end

if validateInputs && ~ray.LineOfSight
    if iscell(mtls)
        % Force mtls to be a homogeneous cell array for code generation
        coder.varsize('mtls', [1 length(mtls)], [0 0]);

        validateattributes(mtls, {'cell'}, {'row'}, ...
            'raypl', 'material input');
        supportedMtls = {'plasterboard','ceiling-board','chipboard', ...
            'floorboard','concrete','brick','wood','glass', 'marble', 'plywood', ...
            'metal','water','vegetation','loam','perfect-reflector'};
        for i = 1:length(mtls)
            coder.internal.errorIf(~any(strcmpi(mtls{i}, supportedMtls)), ...
                'shared_channel:raypl:MaterialNameNotSupported', mtls{i});
        end
    else
        validateattributes(mtls, {'numeric'}, ...
            {'real','nrows',2}, ...
            'raypl', 'material input');
        mtls = double(mtls);
    end
    
    % Check number of materials matches number of reflections
    numMtls = size(mtls, 2);
    if numMtls == ray.NumInteractions
        materials = mtls;
    else
        if numMtls == 1
            materials = repmat(mtls, 1, ray.NumInteractions);
        else
            coder.internal.error('shared_channel:raypl:NumRefNumMtlMismatch');
        end
    end
else
    materials = mtls;
end

% Validate Tx orientation axes
if validateInputs 
    validateAxes(txAxes, 'transmitter');
end

% Validate Rx orientation axes
if validateInputs 
    validateAxes(rxAxes, 'receiver');
end

% Validate Tx polarization input and derive Jones vector from it
[txJV,isTxPol] = getJonesVector(txPol, 'transmitter', validateInputs);

% Validate Rx polarization input and derive Jones vector from it
[rxJV,isRxPol] = getJonesVector(rxPol, 'receiver', validateInputs);

% For unpolarized antennas (e.g. "isotropic"), ignore txAxes' and rxAxes'
% rotation effect on local polarization
if ~isTxPol
    txAxes = eye(3);
end
if ~isRxPol
    rxAxes = eye(3);
end

% Get polarization matrix
polMtx = getPolMatrix(ray, materials, txAxes, rxAxes);

% Calculate total loss due to interactions, depending on polarization
if ~isTxPol || ~isRxPol
    % To handle non-/un-polarized antenna(s), which by definition have
    % random/unknown phase coherency between polarizations, treat both
    % antennas as unpolarized with no polarization mismatch or phase
    % coherency.
    intactLoss = sum(abs(polMtx),'all')/2; % = sqrt(2)/2*[1;1]' * abs(polMtx) * sqrt(2)/2*[1;1]
else
    intactLoss = rxJV' * polMtx * txJV;
end

% Calculate path loss in dB from interactions and free-space
lambda = 299792458/ray.Frequency; % wavelength in free-space
pl = -20*log10(abs((intactLoss))) + ...     % Reflection loss
      fspl(ray.PropagationDistance,lambda); % Free space loss = fspl
pl = max(0, pl);

% Calculate phase change in radians from reflections and free space
numCycles = ray.PropagationDelay*fc;
phase = mod(angle(intactLoss) ...   % Reflection phase change, exp(-iwt) convention
       + 2*pi*numCycles, 2*pi);     % Free space phase change, exp(-iwt) convention
end

function validateAxes(antennaAxes, antennaStr)
% Validate tx/rx orientation axes input

validateattributes(antennaAxes, {'double'}, ...
    {'real','finite','size',[3 3]}, ...
    'raypl', [antennaStr, ' orientation axes input']);
coder.internal.errorIf(...
    max(max(abs(antennaAxes'*antennaAxes-eye(3)))) > sqrt(eps), ...
    'shared_channel:raypl:AxesNotUnitary', antennaStr);

end

function [JV,isPol] = getJonesVector(pol, antennaStr, validateInputs)
% Validate tx/rx polarization input and derive Jones vector in the form of
% [H; V]

if validateInputs
    validateattributes(pol, {'string','numeric'}, {}, ...
        'raypl', [antennaStr, ' polarization input']);
end

if isstring(pol)
    if validateInputs
        validateattributes(pol, {'string'}, {'scalar'}, ...
            'raypl', [antennaStr, ' polarization input']);
        coder.internal.errorIf( ...
            ~any(strcmpi(pol, {'H','V','LHCP','RHCP','none'})), ...
            'shared_channel:raypl:PolTypeNotSupported', antennaStr);
    end
    switch pol  % This matches the top-of-file & online documentation
        case 'H'
            JV = [1; 0];
            isPol = true;
        case 'V'
            JV = [0; 1];
            isPol = true;
        case 'LHCP'
            JV = 1/sqrt(2)*[1; 1i]; % exp(-iwt) transmit convention
            isPol = true;
        case 'RHCP'
            JV = 1/sqrt(2)*[1; -1i]; % exp(-iwt) transmit convention
            isPol = true;
        otherwise % 'none' = randomly polarized
            JV = 1/sqrt(2)*[1; 1];
            isPol = false;
    end
else % Numeric
    if validateInputs
        validateattributes(pol, {'double'}, {'finite','size',[2,1]}, ...
            'raypl', [antennaStr, ' polarization input']);
        coder.internal.errorIf(abs(norm(pol) - 1) > sqrt(eps), ...
            'shared_channel:raypl:JonesVectorNotNormalized', antennaStr);
    end
    JV = pol;
    isPol = true;
end

end

% [EOF]