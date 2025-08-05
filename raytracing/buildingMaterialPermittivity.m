function [epsilon, sigma, complexEpsilon] = buildingMaterialPermittivity(mtls, fc)
%

% Copyright 2019-2024 The MathWorks, Inc.

%#codegen

coder.extrinsic('rfprop.PropagationModel.initializePropModels');

if coder.target('MATLAB')
    rfprop.PropagationModel.initializePropModels;        
else
    coder.const(@feval, 'rfprop.PropagationModel.initializePropModels');
    rfprop.PropagationModel.initializePropModels(); 
end

% Validate material name input
coder.internal.errorIf(~isrow(mtls) || ...
    (~ischar(mtls) && ~iscellstr(mtls) && ~isstring(mtls)), ...
    'shared_channel:buildingMaterialPermittivity:InvalidMaterialDataType');

% Validate carrier frequency input
validateattributes(fc, {'numeric'}, ...
    {'real','positive','finite','scalar'}, ...
    'buildingMaterialPermittivity', 'carrier frequency');

% Cast to double if any
fc = double(fc);

% Table 3 in [1]
mtlLib = {...
    % Name                  a       b       c       d       Fmin    Fmax    Extrap      
    'vacuum'                1.0     0       0.0     0.0     1e6     1e11    1;
    'concrete'              5.24    0       0.0462  0.7822  1e9     1e11    1;
    'brick'                 3.91    0       0.0238  0.16    1e9     40e9    1;
    'plasterboard'          2.73    0       0.0085  0.9395  1e9     1e11    1;
    'wood'                  1.99    0       0.0047  1.0718  1e6     1e11    1;
    'glass',                6.31    0       0.0036  1.3394  1e8     100e9   1;
    'glass',                5.79    0       0.0004  1.658   220e9   450e9   1;
    'ceiling-board'         1.48    0       0.0011  1.0750  1e9     100e9   1;
    'ceiling-board'         1.52    0       0.0029  1.029   220e9   450e9   1;
    'chipboard'             2.58    0       0.0217  0.7800  1e9     1e11    1;
    'plywood'               2.71    0       0.33    0       1e9     40e9    1;
    'marble'                7.074   0       0.0055  0.9262  1e9     60e9    1;
    'floorboard'            3.66    0       0.0044  1.3515  5e10    1e11    1;
    'metal',                1.0     0       1e7     0       1e9     1e11    1;
    'very-dry-ground'       3.0     0       0.00015 2.52    1e9     1e10    0;
    'medium-dry-ground'     15.0    -0.1    0.035   1.63    1e9     1e10    0;
    'wet-ground'            30.0    -0.4    0.15    1.30    1e9     1e10    0};

% Extract out the parameters from the cell array because coder doesn't
% support indexing into a hetergeneous cell array. Also the following line
% is an alternative of mtlLib(:,2:8) which coder doesn't support. 
mtlParams = cell2mat(reshape({mtlLib{:,2:8}}, [], 7));

% The following line is an alternative of mtlLib(:,1)' which coder doesn't
% support. 
mtlNames = {mtlLib{:,1}};

% Convert material names
if coder.target('MATLAB')
    mtlCell = cellstr(string(mtls));
else % mtls cannot be a string array in codegen
    if ~iscellstr(mtls) %#ok<ISCLSTR>
        mtlCell = {mtls};
    else
        mtlCell = mtls;
    end
end

% Find matching rows in the library table
numMtl = length(mtlCell);
matchingIdx = zeros(1, numMtl);
thisFreqLimits = zeros(1, 2);
for mtlIdx = 1:numMtl
    % Find matching material in table
    thisMatchingIdx = find(strcmpi(mtlCell{mtlIdx}, mtlNames));
    if isempty(thisMatchingIdx)
        % Error for an invalid query material
        coder.internal.error('shared_channel:buildingMaterialPermittivity:InvalidMaterial', mtlCell{mtlIdx});
    elseif ~isscalar(thisMatchingIdx)
        % Build an array of frequencies then find closest frequency
        frequencyArray = mtlParams(thisMatchingIdx,5:6)';
        [~,I] = min(abs(frequencyArray - fc), [], "all");
        I = I(end);
        thisMatchingIdx = thisMatchingIdx(ceil(I/2));
    end

    % Error if the query material is one of the 3 "ground" materials and
    % the frequency is outside their frequency range. See page 24 in [1].
    thisMatchingIdx = thisMatchingIdx(1); % codegen requires this to be forced scalar
    thisFreqLimits = mtlParams(thisMatchingIdx,5:6);
    if mtlParams(thisMatchingIdx,7) == 0
        coder.internal.errorIf(fc < thisFreqLimits(1) || fc > thisFreqLimits(2), ...
            'shared_channel:buildingMaterialPermittivity:FcOutOfRangeForGround', ...
            thisFreqLimits(1)/1e9, thisFreqLimits(2)/1e9, mtlCell{mtlIdx});
    end
    matchingIdx(mtlIdx) = thisMatchingIdx;
end

% Calculate permittivity [1, eq. (57)-(59)]
fcGHz = fc/1e9;
epsilon = ([mtlParams(matchingIdx,1)] .* (fcGHz.^[mtlParams(matchingIdx,2)]))';
sigma = ([mtlParams(matchingIdx,3)] .* (fcGHz.^[mtlParams(matchingIdx,4)]))';

% complexEpsilon = epsilon - 1i*sigma/(2*pi*fc*epsilon0). See eq. (9b) in [1]
epsilon0 = 8.854187817e-12; % absolute permittivity of free-space [ITU-R P.527-4]
complexEpsilon = epsilon - 1i*sigma/(2*pi*fc*epsilon0);

end

% [EOF]