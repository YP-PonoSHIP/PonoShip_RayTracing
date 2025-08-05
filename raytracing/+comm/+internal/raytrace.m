function rays = raytrace(env, mtl, txs, rxs, cfg, varargin)
%RAYTRACE Calculate rays between transmitter and receiver
%   RAYS = COMM.INTERNAL.RAYTRACE(ENV, MTL, TXS, RXS, CFG) calculates
%   reflection rays given an environmental description, ENV, surface
%   materials, MTL, transmitter sites, TXS, receiver sites, RXS, and a
%   configuration structure, CFG. The reflection rays are calculated using
%   either the image or SBR method.
%   
%   ENV is either a triangulation object or a structure with three fields: 
%     Triangulation:  Triangulation object for the environmental geometry
%     RayTracer:      Ray tracer object built from the triangulation object
%     SharpEdgeFlags: Logical matrix of the same size as
%                     Triangulation.ConnnectivityList indicating each edge
%                     of every triangle is sharp or not
%   
%   MTL is a row vector of structure that has one field: 
%     Material: Specified as a string scalar for the material name or a
%               double precision, 2-by-1 real vector for the relative
%               permittivity and conductivity. 
%
%   TXS/RXS is a vector of comm.internal.Site objects. CFG is a structure.
%
%   RAYS are a NUMTX-by-NUMRX cell array. Each element is a row vector of
%   comm.Ray objects containing all the rays from the corresponding
%   transmitter position to the corresponding receiver position.

% Copyright 2019-2024 The MathWorks, Inc.

%#ok<*AGROW>

narginchk(6, 7);

% Validate geometry input
validateattributes(env, {'struct','triangulation'}, {'numel',1}, ...
    'raytrace', 'environment');
if isstruct(env)
    TRI = env.Triangulation;
    validateattributes(TRI, {'triangulation'}, {'numel', 1}, ...
        'raytrace', 'triangulation geometry');
    RT = env.RayTracer;
    validateattributes(RT, {'matlabshared.internal.StaticSceneRayTracer'}, {'scalar'}, ...
        'raytrace', 'ray tracer');
    validateattributes(env.SharpEdgeFlags, {'logical'}, {'size',[size(TRI,1), 3]}, ...
        'raytrace', 'sharp edge flags');
else
    TRI = env;
    RT = matlabshared.internal.StaticSceneRayTracer(TRI);
    env = struct( ...
        "Triangulation", TRI, ...
        "RayTracer", RT, ...
        "SharpEdgeFlags", comm.internal.geometry.getSharpEdges(TRI));    
end

% Validate material input
validateattributes(mtl, {'struct'}, {'row'}, ...
    'raytrace', 'material');
for i = 1:length(mtl)
    thisMtl = mtl(i);
    if isstring(thisMtl.Material)
        validateattributes(thisMtl.Material, {'string'}, ...
            {'scalar'}, 'raytrace', 'material name');
    else
        validateattributes(thisMtl.Material, {'double'}, ...
            {'real','finite','size',[2 1]}, ...
            'raytrace', 'material properties');
    end     
end

% Validate Tx and Rx sites
validateattributes(txs, {'comm.internal.Site'}, {'vector'}, ...
    'raytrace', 'transmitter site input');
coder.internal.errorIf(numel([txs.Frequency]) < numel(txs), ...
    'shared_channel:raytrace:EmptyTxFreq');
validateattributes(rxs, {'comm.internal.Site'}, {'vector'}, ...
    'raytrace', 'receiver site input');

% Get number of Tx and Rx
numTx = numel(txs);
numRx = numel(rxs);

% Validate configuration structure
validateattributes(cfg, {'struct'}, {'scalar'}, ...
    'raytrace', 'configuration input');
isSBR = strcmp(cfg.Method, 'sbr');
if ~isSBR % Image method
    refOrder = cfg.NumReflections;
    % Validate reflection order input
    validateattributes(refOrder, {'double'}, ...
        {'real','nonnegative','integer','row'}, ...
        'raytrace', 'reflection orders');    
    refOrder = unique(refOrder);
else % SBR method
    maxRefOrder = cfg.MaxNumReflections;
    % Validate reflection order input
    validateattributes(maxRefOrder, {'double'}, ...
        {'real','nonnegative','integer','scalar'}, ...
        'raytrace', 'max reflection order');    
    
    maxDiffOrder = cfg.MaxNumDiffractions;
    % Validate diffraction order input
    validateattributes(maxDiffOrder, {'double'}, ...
        {'real','nonnegative','integer','scalar'}, ...
        'raytrace', 'max diffraction order');
end

% Validate MaxAbsolutePathLoss and MaxRelativePathLoss
if isfield(cfg, 'MaxAbsolutePathLoss')
    validateattributes(cfg.MaxAbsolutePathLoss, {'double'}, ...
        {'nonsparse','real','positive'});
else
    cfg.MaxAbsolutePathLoss = Inf;
end
if isfield(cfg, 'MaxRelativePathLoss')
    validateattributes(cfg.MaxRelativePathLoss, {'double'}, ...
        {'nonsparse','real','nonnegative'});
else
    cfg.MaxRelativePathLoss = 40;
end

sceneMesh = varargin{1};

if nargin == 7 && ~isempty(varargin{2})
    progressDlg = varargin{2};
    validateattributes(progressDlg, {'matlab.ui.dialog.ProgressDialog','huiprogressdlg'}, ...
        {'scalar'}, 'raytrace', 'ProgressDialog');
    updateProgressDlg = true;
else
    updateProgressDlg = false;
    progressDlg = [];
end

% Initialize return
rays = repmat({comm.Ray.empty}, numTx, numRx);

% Create a ray object "template" to avoid calling constructor for each ray
if ~isSBR
    tempRay = comm.Ray( ...
        'CoordinateSystem', 'Cartesian', ...    
        'PathLossSource', 'Custom');
end

% Process materials which depend on the transmitter frequencies
txFrequency = [txs(:).Frequency];
[uniqueFreq, ~, uniqueIdxFreq] = unique(txFrequency, 'stable');
numUniqueFreq = numel(uniqueFreq);
numMtl = numel(mtl);
mtlTxPermCond = nan(2, numUniqueFreq, numMtl);
eps0 = 8.854187817e-12; % vacuum permittivity (F/m)
cond2PermFactor = (1i/(2*pi*eps0))./uniqueFreq;
mtlPerm = nan(numMtl, numUniqueFreq);
for mtlIdx = 1:numMtl
    thisMaterial = mtl(mtlIdx).Material;

    % Calculate permittivity and update material name as necessary
    if isstring(thisMaterial)
        [epsilon, sigma] = arrayfun(@(mat,f) ...
            comm.internal.factoryMaterialPermittivity(mat, f),  ...
                repmat(thisMaterial, 1, numUniqueFreq), uniqueFreq);
        mtlTxPermCond(:, :, mtlIdx) = [epsilon; sigma];
        mtlPerm(mtlIdx, :) = epsilon - sigma.* ...
            cond2PermFactor; % complex relative permittivity, exp(jwt)
            
        mtl(mtlIdx).MaterialName = mtl(mtlIdx).Material;
    else % No need to update material name
        mtlTxPermCond(:,:,mtlIdx) = repmat(thisMaterial, 1, numUniqueFreq);
        mtlPerm(mtlIdx, :) = thisMaterial(1) - thisMaterial(2).* ...
            cond2PermFactor; % complex relative permittivity, exp(jwt)
    end
end
mtlPerm(isnan(mtlPerm)) = 1 - 1i*Inf;

% Pre-construct SBR variables
if isSBR
    % Angular separation -> tesselation frequency
    angsep = cfg.AngularSeparation;
    if isnumeric(angsep)
        % Compute tessellation frequency which is next power of 2 from target
        targetTessFreq = 69/angsep;
        tessFreq = pow2(nextpow2(targetTessFreq));
    else
        switch angsep
            case "low"
                tessFreq = 256;
            case "medium"
                tessFreq = 128;
            otherwise % 'high'
                tessFreq = 64;
        end
    end

    % Material triangular indices
    if numMtl == 1 % StaticScenePathTracer has simpler input for 1 net material
        mtlTriIdx = 1;
    else
        [~, tri] = Geometry.extractMesh(sceneMesh);
        mtlTriIdx = tri(:,4);
    end

    % Transmit and receive antenna objects for path tracer
    eye3 = eye(3);
    txAntenna = matlabshared.internal.Antenna;
    txAntenna.Transform = [eye3 [1; 1; 1]; 0 0 0 1];
    for rxIdx = numRx:-1:1
        rxAntennas(rxIdx) = matlabshared.internal.Antenna;
        rxAntennas(rxIdx).Transform = [eye3 rxs(rxIdx).Position; 0 0 0 1];
    end

    % Pre-allocated ray variables
    tempRayObject = comm.Ray( CoordinateSystem='Cartesian', ...
        LineOfSight=false, PathLossSource='Custom');
    tempRayInt = repmat(struct('Type', 'Reflection', ...
        'Location', [10; 10; 0], 'MaterialName', mtl(1).MaterialName), ...
        1, cfg.MaxNumReflections + cfg.MaxNumDiffractions);
    emptyRays = repmat({comm.Ray.empty}, 1, numRx);

    % Initialize loop variables
    thisTxFrequency = nan;
    thisTxPosition = nan;
end

% UseGPU
if isfield(cfg, "UseGPU")
    useGPU = cfg.UseGPU;
else
    useGPU = "off";
end

% Iterate through each TX
progBarInterval = 0:(1/numTx):1;
for txIdx = 1:numTx  
    thisTx = txs(txIdx);
    thisUniqueIdxFreq = uniqueIdxFreq(txIdx);

    % Update materials for this TX frequency
    thisMtl = mtl;
    for i = 1:numMtl
        thisMtl(i).Material = mtlTxPermCond(:, thisUniqueIdxFreq, i);
    end

    % Perform ray-tracing
    if isSBR % SBR method
        % Update transmit frequency
        prevTxFrequency = thisTxFrequency;
        thisTxFrequency = txFrequency(txIdx);
        if ~isequal(thisTxFrequency, prevTxFrequency)
            % Create path tracer object
            pt = matlabshared.internal.StaticScenePathTracer(env.Triangulation, ...
                mtlTriIdx, struct('permittivity', num2cell(mtlPerm(:,thisUniqueIdxFreq))), tessFreq);
            pt.Algorithm = matlabshared.internal.StaticScenePathTracerAlgorithm.Hybrid;
            pt.MaxNumReflections = cfg.MaxNumReflections;
            pt.MaxNumDiffractions = cfg.MaxNumDiffractions;
            pt.MaxAbsolutePathLoss = cfg.MaxAbsolutePathLoss;
            pt.MaxRelativePathLoss = cfg.MaxRelativePathLoss;
            pt.UseGPU = useGPU;

            % Update txAntenna and tempRayObject
            txAntenna.Frequency = thisTxFrequency;
            tempRayObject.Frequency = thisTxFrequency;
        end

        % Update transmit position
        prevTxPosition = thisTxPosition;
        thisTxPosition = thisTx.Position;
        if ~isequal(thisTxPosition, prevTxPosition)
            txAntenna.Transform(1:3,4) = thisTxPosition;
            tempRayObject.TransmitterLocation = thisTxPosition;
        end        

        % Perform SBR
        rays(txIdx, :) = findRaysSBR(thisMtl, thisTx, rxs, ...
            updateProgressDlg, progressDlg, progBarInterval(txIdx:txIdx+1), sceneMesh, ...
            pt, txAntenna, rxAntennas, tempRayObject, tempRayInt, emptyRays);
    else % Image method
        % Get viewshed for this Tx
        if any(refOrder)
            temp = comm.internal.geometry.viewshed(TRI, RT, thisTx.Position');
            visFacetTx = temp{1};
        else
            visFacetTx = [];
        end
        
        % Iterate through each Rx
        for rxIdx = 1:numRx
            % Skip this pair if Tx and Rx are at the same position
            if (norm(thisTx.Position - rxs(rxIdx).Position) < sqrt(eps))
                continue;
            end
            
            % Calculate rays for this pair of Tx and Rx
            rays{txIdx,rxIdx} = findRaysImage( ...
                env, thisMtl, thisTx, visFacetTx, ...
                rxs(rxIdx), refOrder, tempRay, sceneMesh);
            
            % Update progress bar
            if updateProgressDlg
                totIdx = ((txIdx-1)*numRx+rxIdx);
                if mod(totIdx,rfprop.Constants.ProgressDialogUpdateSkipFactor) == 0
                    % Update progress dialog using skip factor
                    updateProgressBar(progressDlg, totIdx/(numTx*numRx));
                else
                    updateProgressBar(progressDlg);
                end
            end
        end
    end
end

end

function rays = findRaysImage(env, mtl, tx, visFacetTx, rx, refOrder, tempRay, sceneMesh)

% Parse site inputs
txPos = tx.Position';
rxPos = rx.Position';

% Parse material input
mtlProp = [mtl.Material];

% Update template ray object
tempRay.TransmitterLocation = txPos';
tempRay.ReceiverLocation = rxPos';
tempRay.Frequency = tx.Frequency;

% Check LOS ray
if any(0 == refOrder) 
    [~, ~, isNLOS] = ...
        comm.internal.geometry.firstIntersect(env, txPos, rxPos, 'segment');
    if ~isNLOS
        rays = tempRay;
    else
        rays = comm.Ray.empty;
    end
else
    rays = comm.Ray.empty;
end

% All the rays to be found are NLOS
tempRay.LineOfSight = false;
tempRay.Interactions = comm.internal.RayInteraction;

% If only one material is being used, we can avoid using expensive
% geometry library indexing.
usingSingleMaterial = isscalar(mtl);

% Desired reflection orders to search other than LOS
if ~isempty(visFacetTx)
    % Iterate through the desired reflection orders
    for i = (refOrder(1) == 0)+1:length(refOrder)
        % Find the reflection rays recursively. The output refPts is of
        % size [numRays, 3, numRef]; The output refFacetIdx is of size
        % [numRays, numRef].
        [refPts, refFacetIdx] = ...
            comm.internal.geometry.calcImageReflectionPath( ... 
            env, txPos, visFacetTx, rxPos, refOrder(i), refOrder(i));
        
        % Formulate ray objects and materials for this order of reflections
        numRays = size(refPts, 1);       
        thisOrderRays = repmat(tempRay, 1, numRays);
        
        % Create interaction objects to capture all the reflections
        % including material information
        for k = 1:numRays
            numRefFacet = numel(refFacetIdx(k,:));
            refMtlIdx = zeros(1, numRefFacet);
            interactionMaterialNames = cell(1,numRefFacet);
            % There can be more than one interaction triangle so we
            % must iterate over the triangle IDs to aggregate materials
            for k2 = 1:numRefFacet
                if usingSingleMaterial
                    materialIndex = 1;
                else
                    interactionTriangle = sceneMesh.GetTriangle(refFacetIdx(k,k2)-1);
                    materialIndex = interactionTriangle.GetTriangleMaterialIndex;
                end
                interactionMaterialNames{k2} = mtl(materialIndex).MaterialName;
                refMtlIdx(k2) = materialIndex;
            end
            thisOrderRays(k).Interactions = ...
                comm.internal.RayInteraction.createRefPath( ...
                shiftdim(refPts(k,:,:)), mtlProp(:, refMtlIdx), ...
                txPos', rxPos');
            for k2 = 1:numel(thisOrderRays(k).Interactions)
                thisOrderRays(k).Interactions(k2).MaterialName = interactionMaterialNames{k2};
            end
        end
        
        % Sort all rays of this order according to distance
        if numRays > 1
            [~, idx] = sort([thisOrderRays.PropagationDistance]);    
            thisOrderRays = thisOrderRays(idx);
        end
        
        % Concatenate all the rays that have been found
        rays = [rays, thisOrderRays];  
    end
end

% Calculate path loss and phase shift for each ray
if ~isempty(rays)
    mtls = interactionMaterials(rays);
    rays = calcPathLossSingleSite(tx, rx, rays, mtls);
end

end

function rays = findRaysSBR(mtl, tx, rxs, updateProgressDlg, ...
    progressDlg, progBarInterval, sceneMesh, pt, txAntenna, rxAntennas, ...
    tempRayObject, tempRayInt, rays)

% Perform ray tracing in C++
try
    ptOut = trace(pt, txAntenna, rxAntennas);
catch
    error(message('shared_channel:raytrace:TraceError'))
end

% Find number of rays per receiver
numRx = numel(rxs);
numRaysPerRx = nan(numRx,1);
for rxIdx = 1:numRx
    numRaysPerRx(rxIdx) = length(ptOut{rxIdx});
end
tempRxRays = repmat(tempRayObject, 1, max(numRaysPerRx));
numRaysTotal = sum(numRaysPerRx);

% Assume finding all the rays (geometric part) takes 50% of the total time
if updateProgressDlg
    progBarDiff = diff(progBarInterval);
    progBarHalfPoint = progBarInterval(1) + 0.5*progBarDiff;
    updateProgressBar(progressDlg, progBarHalfPoint);

    % Variables for incrementing dialog every 1%
    dialogUpdateMinPercent = 0.01;
    numRaysProcessed = 0;
    numRaysDialogFactor = 0.5*progBarDiff;
    if numRaysDialogFactor <= dialogUpdateMinPercent
        % Not worth updating during loop
        numRaysUpdateDialog = Inf;
    else
        numRaysDialogFactor = numRaysDialogFactor/numRaysTotal;
        numRaysIncDialog = max(1, floor(dialogUpdateMinPercent ...
            / numRaysDialogFactor));
        numRaysUpdateDialog = 0;
    end
end

% Parse the results into comm.Ray objects for each Rx
sqrtEps = sqrt(eps);
hasMultipleMaterialsScene = ~isscalar(mtl);
for rxIdx = 1:numRx
    % Skip if no rays
    thisNumRays = numRaysPerRx(rxIdx);
    if thisNumRays == 0
        continue
    end

    % Skip this pair if Tx and Rx are at the same position
    if norm(tx.Position - rxs(rxIdx).Position) < sqrtEps
        continue;
    end

    % Parse rays
    thixRxRays = tempRxRays(1:thisNumRays);
    raysPol = nan(4,3,thisNumRays);
    for rayIdx = 1:thisNumRays
        % Get all the interactions including Tx & Rx for this ray
        thisPtOut = ptOut{rxIdx}(rayIdx);
        thisPtOutInt = thisPtOut.Interactions;
        
        % Log Rx location 
        thisPtOutIntLocationArray = vertcat(thisPtOutInt.Point)';
        thixRxRays(rayIdx).ReceiverLocation = thisPtOutIntLocationArray(:,end);

        % Log 3x3 polarization matrix, path loss, and phase
        raysPol(1:3,1:3,rayIdx) = thisPtOut.Polarization;
        raysPol(4,1:2,rayIdx) = [thisPtOut.PathLoss thisPtOut.PhaseShift];
        
        % Number of interactions excluding Tx and Rx 
        numInt = length(thisPtOutInt) - 2;
        if numInt == 0 % LOS
            thixRxRays(rayIdx).LineOfSight = true;
            continue
        end

        % Log interactions for NLOS ray
        thisRayInt = tempRayInt(1:numInt);
        for intIdx = 1:numInt
            % Log position
            thisRayInt(intIdx).Location = thisPtOutIntLocationArray(:,intIdx+1);

            % Log interaction type
            thisInt = thisPtOutInt(intIdx+1);
            if thisInt.Type == 4 % Diffraction
                thisRayInt(intIdx).Type = 'Diffraction';
            end

            % Log materials; diffractions can have two unique materials
            if hasMultipleMaterialsScene % Multiple materials in scene
                % Get material indices
                primID = thisInt.Interaction.PrimitiveID;
                numPrimID = numel(primID);
                materialIdx = ones(numPrimID, 1);
                for idx = 1:numPrimID
                    interactionTriangle = sceneMesh.GetTriangle(primID(idx) - 1); % Subtract 1 because primID is MATLAB style 1-based indexing, while geometrylibrary is 0 based indexing
                    materialIdx(idx) = interactionTriangle.GetTriangleMaterialIndex;
                end
                materialIdx(materialIdx==0) = 1;
                
                % Log material(s)
                materialIdx = unique(materialIdx, 'stable');
                if isscalar(materialIdx) % Single material
                    thisRayInt(intIdx).MaterialName = mtl(materialIdx).MaterialName;
                else
                    thisRayInt(intIdx).MaterialName = mtl(materialIdx(1)).MaterialName + ...
                        " / " + mtl(materialIdx(2)).MaterialName;
                end
            end
        end
        thixRxRays(rayIdx).Interactions = thisRayInt;
    end

    % Calculate path loss and phase shift
    rays{rxIdx} = calcPathLossSingleSite(tx, rxs(rxIdx), thixRxRays, raysPol);

    % Update progress bar or, for point-to-point analysis, sort rays by
    % interaction type and then by propagation distance
    if updateProgressDlg
        numRaysProcessed = numRaysProcessed + thisNumRays;
        if numRaysProcessed >= numRaysUpdateDialog
            numRaysUpdateDialog = numRaysProcessed + numRaysIncDialog;
            updateProgressBar(progressDlg, progBarHalfPoint + ...
                numRaysProcessed*numRaysDialogFactor);
        else
            updateProgressBar(progressDlg);
        end
    else
        rays{rxIdx} = sortRays(rays{rxIdx});   
    end
end

% Update progress bar to interval end
if updateProgressDlg
    updateProgressBar(progressDlg, progBarInterval(2));
end

end

function rays = calcPathLossSingleSite(tx, rx, rays, mtl_raysPol)
% Calculate path losses for all rays between one tx and one rx. mtl_raysPol
% is an array of materials for image method, or a 3x3 polarization matrix
% for SBR method.

if isempty(rays)
    return
end

fc = tx.Frequency;
txAxes = tx.OrientationAxes;
rxAxes = rx.OrientationAxes;

% Calculate Jones vector
bothAntPol = all([isAntennaPolarized(tx.Antenna);
    isAntennaPolarized(rx.Antenna)]);
if bothAntPol
    % Tx
    globalSph = [[rays.AngleOfDeparture]; ones(size(rays))];
    localSph = global2localcoord( ...
        globalSph, 'ss', [0; 0; 0], txAxes);
    txJV = comm.internal.calcJonesVector( ...
        tx.Antenna, fc, localSph(1:2,:), true);

    % Rx
    globalSph = [[rays.AngleOfArrival]; ones(size(rays))];
    localSph = global2localcoord( ...
        globalSph, 'ss', [0; 0; 0], rxAxes);
    rxJV = comm.internal.calcJonesVector( ...
        rx.Antenna, fc, localSph(1:2,:), true);
end

% Calculate path loss and phase shift for each ray via raypl for image
% method or 3x3 polarization matrix for SBR method
if iscell(mtl_raysPol) % Image method; use raypl & materials
    if bothAntPol
        for k = 1:length(rays)
            [rays(k).PathLoss, rays(k).PhaseShift] = raypl(rays(k), ...
                "TransmitterPolarization", txJV(:,k), ...
                "ReceiverPolarization",    rxJV(:,k), ...
                "TransmitterAxes",         txAxes, ...
                "ReceiverAxes",            rxAxes, ...
                "ReflectionMaterials",     mtl_raysPol{k}, ...
                "ValidateInputs",          false);
        end
    else % Tx or Rx is unpolarized. Treat both to be unpolarized.
        for k = 1:length(rays)
            [rays(k).PathLoss, rays(k).PhaseShift] = raypl(rays(k), ...
                "TransmitterPolarization", 'none', ...
                "ReceiverPolarization",    'none', ...
                "TransmitterAxes",         txAxes, ...
                "ReceiverAxes",            rxAxes, ...
                "ReflectionMaterials",     mtl_raysPol{k}, ...
                "ValidateInputs",          false);
        end
    end
else % SBR method; use 3x3 polarization matrix
    two_pi = 2*pi;
    if bothAntPol
        lambda = 299792458/fc; % wavelength in free-space
        for k = 1:length(rays)
            thisRay = rays(k);

            % Get the ray direction (out-going, in-coming) relative to
            % transmitter and receiver
            if thisRay.LineOfSight
                % For this case, out-going direction = in-coming direction
                rxDir = thisRay.ReceiverLocation - thisRay.TransmitterLocation;
                rxDir = rxDir/norm(rxDir);
                txDir = rxDir;
            else
                rxDir = thisRay.ReceiverLocation - thisRay.Interactions(end).Location;
                rxDir = rxDir/norm(rxDir);
                txDir = thisRay.Interactions(1).Location - thisRay.TransmitterLocation;
                txDir = txDir/norm(txDir);
            end

            % Path loss
            intactLoss = jonesToGlobal(rxJV(:,k), rxDir')' ...
                * mtl_raysPol(1:3,1:3,k) * jonesToGlobal(txJV(:,k), txDir');
            rays(k).PathLoss = -20*log10(abs((intactLoss))) + ...     % Reflection loss
                fspl(rays(k).PropagationDistance,lambda); % Free space loss = fspl

            % Phase shift
            numCycles = rays(k).PropagationDelay*fc;
            rays(k).PhaseShift = mod(-angle(intactLoss) ... % Reflection phase change, exp(jwt) convention
                + two_pi*numCycles, two_pi); % Free space phase change, exp(-iwt) convention
        end
    else
        for k = 1:length(rays)
            rays(k).PathLoss = mtl_raysPol(4,1,k);
            rays(k).PhaseShift = mod(-mtl_raysPol(4,2,k), two_pi); % Free space phase change, exp(-iwt) convention
        end
    end
end

end

function E = jonesToGlobal(J, dir)
% Convert Jones Vector & ray direction into a global 3x3 matrix
% J = [H; V]
% dir = [x y z], normalized, in global coordinates

% Create vertical & horizontal polarization basis vectors for this
% direction vector
[az, el] = cart2sph(dir(1), dir(2), dir(3));
[x, y, z] = sph2cart(az, el+pi/2, 1);
vPolDir = [x y z];
hPolDir = cross(dir, vPolDir);

E = [hPolDir; vPolDir; dir]' * [J;0];
end

function isPol = isAntennaPolarized(antenna)
% Determine if antenna is polarized
if isa(antenna,'phased.internal.AbstractAntennaElement') || ...
        isa(antenna,'phased.internal.AbstractArray') || ...
        isa(antenna,'phased.internal.AbstractSubarray')
    % Elements or arrays in PST
    isPol = isPolarizationCapable(antenna);
elseif isa(antenna, 'em.Antenna') || isa(antenna, 'em.Array') || ...
       isa(antenna,'installedAntenna') || isa(antenna,'customAntennaStl')
    % Elements or arrays in Antenna Toolbox
    isPol = true;
else % Default 'isotropic' or arrayConfig object
    isPol = false;
end
end

function updateProgressBar(progressDlg, varargin)

% Halt if cancel requested
if progressDlg.CancelRequested
    error(message('shared_channel:rfprop:ProgressDialogCancelled'));
end

% Update progress dialog
if ~isempty(varargin)   
    progressDlg.Value = varargin{1};
end

end

function sortedRays = sortRays(rays)
% Sort the rays with the same sequence of interaction types by their
% propagation distance 

numRays = length(rays);

if numRays < 2
    sortedRays = rays;
    return;
end

% Find and sort unique sequence of interaction types
intTypeSeq = cell(1, numRays);
for i = 1:numRays
    if rays(i).LineOfSight
        intTypeSeq{i} = '';
    else
        intTypeSeq{i} = [rays(i).Interactions.Type];
    end
end
uniIntTypeSeq = sort(unique(intTypeSeq));

% Sort the rays with the same sequence
sortedRays = rays;
curIdx = 0;
for i = 1:length(uniIntTypeSeq)
    rayIdx = strcmp(uniIntTypeSeq{i}, intTypeSeq);
    thisRays = rays(rayIdx);
    [~, idx] = sort([thisRays.PropagationDistance]);    
    sortedRays(curIdx+(1:length(thisRays))) = thisRays(idx);
    curIdx = curIdx + length(thisRays);
end

end

function mtls = interactionMaterials(rays)

numRays = numel(rays);
mtls = cell(1,numRays);
for rayInd = 1:numRays
    numInt = rays(rayInd).NumInteractions;
    mtls{rayInd} = cell(1,numInt);
    for intInd = 1:numInt
        mtls{rayInd}{intInd} = rays(rayInd).Interactions(intInd).Material;
    end
end
end