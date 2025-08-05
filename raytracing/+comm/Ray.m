classdef Ray < comm.internal.ConfigBase 
%

% Copyright 2019-2024 The MathWorks, Inc.
%#codegen

properties    
    PathSpecification(1,:) char = 'Locations'
    CoordinateSystem(1,:) char = 'Cartesian'
    SystemScale(1,1) {mustBeNonsparse, mustBeNumeric, mustBeFinite, mustBePositive} = 1
    TransmitterLocation {mustBeNonsparse, mustBeVector, mustBeNumeric, ...
        mustBeFinite, mustBeReal, mustBeVector} = [0; 0; 0] % [lat; lon; h] or [x; y; z]
    ReceiverLocation {mustBeNonsparse, mustBeVector, mustBeNumeric, ...
        mustBeFinite, mustBeReal} = [10; 10; 10] % [lat; lon; h] or [x; y; z]
    LineOfSight(1,1) logical = true
    Interactions = struct('Type', 'Reflection', ...
        'Location', [10; 10; 0], 'MaterialName', "")
end

properties (Dependent)
    PropagationDelay(1,1) double {mustBeFinite, mustBeNonnegative} % s
end

properties (SetAccess = private)
    PropagationDistance(1,1) double {mustBeFinite, mustBeNonnegative} % m
end

properties (Dependent)
    AngleOfDeparture(2,1) double {mustBeFinite, validateAngleRange(AngleOfDeparture)} % [az; el] °
    AngleOfArrival(2,1) double {mustBeFinite, validateAngleRange(AngleOfArrival)} % [az; el] °
end

properties
    Frequency(1,1) {mustBeNonsparse, mustBeNumeric, mustBeFinite, mustBePositive} = 1.9e9 % Hz
    PathLossSource(1,:) char {mustBeMember(PathLossSource, {'Free space model'; 'Custom'})} = 'Free space model'
end

properties (SetAccess = private)
    NumInteractions(1,1) double {mustBeNonnegative, mustBeInteger}
end

properties (Dependent)
    PathLoss(1,1) double {mustBeNonnegative} % dB
    PhaseShift(1,1) double {mustBeFinite, mustBeReal} % rad
end

properties (Dependent, Hidden)
    ReflectionLocations
end

properties (SetAccess = private, Hidden)
    NumReflections
    NumDiffractions
end

properties (Access = private)    
    pIsLocked
    pPropDelay
    pAoD
    pAoA
    pPathLoss
    pPhaseShift
    AngleOfIncidence % Angle between the incident ray and facet normal of each interaction
end

properties(Hidden)
    UID
end

methods
    function obj = Ray(varargin)
        coder.extrinsic('rfprop.PropagationModel.initializePropModels');
        obj@comm.internal.ConfigBase(varargin{:});

        % Default values for non-dependent properties
        PVPairs = ...
            {'PathSpecification',   'Locations'; ...
            'CoordinateSystem',    'Cartesian'; ...
            'SystemScale',         1; ...
            'TransmitterLocation', zeros(3, 1); ...
            'ReceiverLocation',    10*ones(3, 1); ...
            'LineOfSight',         true; ...
            'Interactions',        struct('Type','Reflection','Location',[10;10;0],'MaterialName',""); ...
            'Frequency',           1.9e9; ...
            'PathLossSource',      'Free space model'; ...
            'PropagationDistance', sqrt(300); ...
            'NumInteractions',     0; ...
            'NumReflections',      0; ...
            'NumDiffractions',      0; ...
            'pIsLocked',           false; ...
            'pPropDelay',          sqrt(300)/physconst('lightspeed'); ...
            'pAoD',                [45; rad2deg(atan2(10, 10*sqrt(2)))]; ...
            'pAoA',                [-135; -rad2deg(atan2(10, 10*sqrt(2)))]; ...
            'pPathLoss',           0; ...
            'pPhaseShift'          0; ...
            'AngleOfIncidence',    0; ...
            'UID',                 0};

        if isempty(coder.target)
            rfprop.PropagationModel.initializePropModels;
            % Assign default values if not defined yet, except for UID which
            % we don't need in codegen
            for i = coder.unroll(1:size(PVPairs, 1) - 1)
                if isempty(obj.(PVPairs{i,1}))
                    obj.(PVPairs{i,1}) = PVPairs{i,2};
                end
            end
            obj.UID = obj.createUIDNumber;
        else
            coder.const(@feval, 'rfprop.PropagationModel.initializePropModels');
            % Not do this when the object is being const-folded, which
            % typically happpens when it is being assigned to a System object
            % non-tunable property.
            if ~coder.internal.isConstantFolding
                rfprop.PropagationModel.initializePropModels();
            end
            % Assign default values if not defined yet
            for i = coder.unroll(1:size(PVPairs, 1))
                if ~coder.internal.is_defined(obj.(PVPairs{i,1}))
                    obj.(PVPairs{i,1}) = PVPairs{i,2};
                end
            end
        end
    end

    % Set methods for public properties
    function obj = set.CoordinateSystem(obj, val)
        obj.CoordinateSystem = validatestring(val, ...
            {'Cartesian'; 'Geographic'}, 'CoordinateSystem');
    end
    function obj = set.Frequency(obj, val)
        obj.Frequency = double(val);
    end
    function obj = set.Interactions(obj, val)
        arguments
            obj
            val {mustBeA(val, {'comm.internal.RayInteraction'; 'struct'})}
        end

        % Prepare validation function to be compatible with CodeGen
        coder.extrinsic('comm.Ray.validateInteractionsStruct');

        if isstruct(val)
            validateattributes(val, {'struct'}, {'row'});

            % Validate individual structure
            if ~coder.internal.isConstantFolding
                comm.Ray.validateInteractionsStruct(val);
            end
        end

        obj.Interactions = val;
    end
    function obj = set.PathSpecification(obj, val)
        obj.PathSpecification = validatestring(val, ...
            {'Locations'; 'Delay and angles'}, 'PathSpecification');
    end
    function obj = set.ReceiverLocation(obj, val)
        validateattributes(val, {'numeric'}, {'numel',3});
        obj.ReceiverLocation = double(val(:)); % column vector
    end
    function obj = set.SystemScale(obj, val)
        obj.SystemScale = double(val);
    end
    function obj = set.TransmitterLocation(obj, val)
        validateattributes(val, {'numeric'}, {'numel',3});
        obj.TransmitterLocation = double(val(:)); % column vector
    end
end

methods % Set and get methods for dependent properties
  function obj = set.PropagationDelay(obj, propDelay)    
    performReadOnlyCheck(obj, 'PropagationDelay');  
    obj.pPropDelay = propDelay;
  end
  
  function propDelay = get.PropagationDelay(obj)
    if isPropertyReadOnly(obj, 'PropagationDelay')
        path = convertPathToCartesian(obj, obj.propagationPath);
        distance = sum(sqrt(sum(diff(path, 1, 2).^2)));
        c = 299792458; % = physconst('lightspeed') w/o validation
        propDelay = distance/c;
    else
        propDelay = obj.pPropDelay;
    end
  end
  
  function obj = set.AngleOfDeparture(obj, aod)    
    performReadOnlyCheck(obj, 'AngleOfDeparture');
    obj.pAoD = aod;
  end
  
  function aod = get.AngleOfDeparture(obj)
    if isPropertyReadOnly(obj, 'AngleOfDeparture')        
        path = obj.propagationPath;
        aod = localangle(obj, path(:,1), path(:,2));
    else
        aod = obj.pAoD;
    end
  end
  
  function obj = set.AngleOfArrival(obj, aoa)
    performReadOnlyCheck(obj, 'AngleOfArrival');
    obj.pAoA = aoa;
  end
  
  function aoa = get.AngleOfArrival(obj)
    if isPropertyReadOnly(obj, 'AngleOfArrival')
        path = obj.propagationPath;
        aoa = localangle(obj, path(:,end), path(:,end-1));
    else
        aoa = obj.pAoA;
    end
  end  
  
  function obj = set.PathLoss(obj, pl)
    performReadOnlyCheck(obj, 'PathLoss');
    obj.pPathLoss = pl;
  end
  
  function pl = get.PathLoss(obj)
    if isPropertyReadOnly(obj, 'PathLoss')
        c = 299792458; % = physconst('lightspeed') w/o validation
        fc = obj.Frequency;
        pl = fspl(obj.PropagationDelay*c, c/fc);
    else
        pl = obj.pPathLoss; 
    end
  end
  
  function obj = set.PhaseShift(obj, phase)
    performReadOnlyCheck(obj, 'PhaseShift');
    obj.pPhaseShift = phase;
  end

  function phase = get.PhaseShift(obj)
    if isPropertyReadOnly(obj, 'PhaseShift')
        fc = obj.Frequency; 
        phase = mod(2*pi*fc*obj.PropagationDelay, 2*pi);
    else
        phase = obj.pPhaseShift; 
    end
  end
end

methods % Get methods for read-only properties
  function dist = get.PropagationDistance(obj)
    c = 299792458; % = physconst('lightspeed') w/o validation
    dist = obj.PropagationDelay * c;
  end
    
  function numIntact = get.NumInteractions(obj)
    numIntact = (1 - obj.LineOfSight) * numel(obj.Interactions);
  end  
end

methods % Set/get methods for deprecated properties 
  function obj = set.ReflectionLocations(obj, refLoc)
    propName = 'ReflectionLocations';
    validateattributes(refLoc, {'double'}, ...
        {'real','2d','finite','nrows',3}, ...
        [class(obj) '.' propName], propName);
    
    numRef = size(refLoc,2);
    intact = repmat( ...
        struct('Type','Reflection','Location',[10;10;0],'MaterialName',""), 1, numRef);
    for i = 1:numRef
        intact(i).Location = refLoc(:,i);
    end
    obj.Interactions = intact;
  end

  function refLoc = get.ReflectionLocations(obj)
    % This is for backward compatibility. As there was no diffractions in
    % the old code, all the locations retrieved are for reflections.      
    refLoc = [obj.Interactions.Location];
  end
  
  function nor = get.NumReflections(obj)
    % Use [obj.Interactions.Type] instead of {obj.Interactions.Type} for
    % code generation to work 
    nor = (1 - obj.LineOfSight) * ...
        sum('R' == convertStringsToChars([obj.Interactions.Type]), 2);
  end

  function nor = get.NumDiffractions(obj)
    % Use [obj.Interactions.Type] instead of {obj.Interactions.Type} for
    % code generation to work 
    nor = (1 - obj.LineOfSight) * ...
        sum('D' == convertStringsToChars([obj.Interactions.Type]), 2);
  end
end

methods % Get methods for private properties
  function aor = get.AngleOfIncidence(obj)
    if obj.NumInteractions
        path = convertPathToCartesian(obj, obj.propagationPath);
        % Get direction of each segment
        pathDir = diff(path, 1, 2); 
        % Normalize segment direction
        pathDir = pathDir ./ sqrt(sum(pathDir.^2, 1));
        % Dot product between 2 consecutive segments
        aor = acosd(dot(pathDir(:,2:end), -pathDir(:,1:end-1)))/2; 
    else
        aor = NaN;
    end
  end
end

methods (Hidden)
    function UID = createUIDNumber(~)
        %createUIDNumber   Create unique identifier number
        
        persistent uidCounter;
        if isempty(uidCounter)
            uidCounter = 0;
        end
        
        % Set counter value to input or else increment if 'next'
        uidCounter = uidCounter + 1;
        UID = uidCounter;
        UID = sprintf('rayObj%u', UID);
    end
end

methods % Public methods
  function plot(rays, varargin)
      %

    % This method doesn't support code generation
    coder.internal.errorIf(~isempty(coder.target), ...
        'shared_channel:Ray:NoCodegenForPlot');    
    
    % Validate rays
    validateattributes(rays,{'comm.Ray'},{'nonempty'},'plot','',1);
    validateRaysForPlot(rays);
    
    % Process optional name/value pairs
    p = inputParser;
    p.addParameter('Animation', '');
    p.addParameter('EnableWindowLaunch', true);
    p.addParameter('TransmitterSite', txsite);
    p.addParameter('ReceiverSite', rxsite);
    p.addParameter('AssociateWithSites', false);
    p.addParameter('Type', 'pathloss');
    p.addParameter('ColorLimits', []);
    p.addParameter('Colormap', 'jet');
    p.addParameter('ShowLegend', true);
    p.addParameter('Map', []);    
    p.parse(varargin{:});

    % Get Site Viewer and validate web graphics
    viewer = rfprop.internal.Validators.validateMap(p, 'plot', strcmpi(rays(1).CoordinateSystem, 'cartesian'));
    isViewerInitiallyVisible = viewer.Visible;
    [animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'plot');
    
    % Validate Type and and get corresponding default parameter values
    allowedTypes = {'power','pathloss'};
    [type, ~, defaultColorLimits, legendTitle] = ...
        rfprop.internal.Validators.validateType(p, allowedTypes, 'plot');
    isPathlossPlot = strcmp(type,'pathloss');
    
    % Validate and get graphical parameters
    cmap = rfprop.internal.Validators.validateColorMap(p, 'plot');
    if isPathlossPlot
        % Flip colormap to maintain consistency of color meaning (red end of
        % jet means good signal), since path loss is inversely proportional to
        % received power
        cmap = flipud(cmap);
    end
    clim = rfprop.internal.Validators.validateColorLimits(p, defaultColorLimits, 'plot');
    showLegend = rfprop.internal.Validators.validateShowLegend(p, 'plot');
    
    % Validate and get sites
    numRays = numel(rays);
    rayTxs = validateSites('txsite', p.Results.TransmitterSite, numRays);
    rayRxs = validateSites('rxsite', p.Results.ReceiverSite, numRays);
    associateWithSites = p.Results.AssociateWithSites;
    
    % Initialize graphics viewer
    graphicsToRemove = {};
    if (~isempty(viewer.LegendID))
        graphicsToRemove{end+1} = viewer.LegendID;
        viewer.LegendID = '';
        viewer.LegendEntities = strings(0);
    end

    for rayInd = 1:numRays
        if (associateWithSites)
            % Go through the viewer's SiteGraphics and remove all rays on
            % the input antenna sites
            rayTx = rayTxs(rayInd);
            rayRx = rayRxs(rayInd);
            % Remove existing rays
            txRayGraphics = viewer.getGraphic(rayTx.UID, 'rays');
            if (isfield(txRayGraphics, rayRx.UID))
                graphicsToRemove = [graphicsToRemove; txRayGraphics.(rayRx.UID)]; %#ok<AGROW>
                viewer.disassociateSiteGraphics(rayTx.UID, rayRx.UID, 'rays');
            end
            % Remove existing LOS lines
            txLOSGraphics = viewer.getGraphic(rayTx.UID, 'los');
            if (isfield(txLOSGraphics, rayRx.UID))
                graphicsToRemove = [graphicsToRemove; txLOSGraphics.(rayRx.UID)]; %#ok<AGROW>
                viewer.disassociateSiteGraphics(rayTx.UID, rayRx.UID, 'los');
            end
        else
            % Go through Site Graphics and remove all plotted rays on the
            % input ray objects
            rayObj = rays(rayInd);
            rayGraphics = viewer.getGraphic(rayObj.UID, 'rays');
            if (~isempty(rayGraphics))
                graphicsToRemove = [graphicsToRemove; rayGraphics]; %#ok<AGROW>
            end
        end
    end
    
    % Add graphics for each ray
    noPathsMeetThreshold = true;
    
    % Queue the plots so the graphics don't update every frame
    viewer.Visualizer.queuePlots;
    for rayInd = 1:numRays
        ray = rays(rayInd);
        rayTx = rayTxs(rayInd);
        rayRx = rayRxs(rayInd);
        
        % Get ray geometry info
        pathDistance = ray.PropagationDistance;
        numIntact = ray.NumInteractions;
        aod = ray.AngleOfDeparture;
        aoa = ray.AngleOfArrival;
        
        % Compute data value for ray
        pl = ray.PathLoss;
        if isPathlossPlot
            v = pl;
        else
            % Compute signal strength using Friis equation
            fq = rayTx.TransmitterFrequency;
            Ptx_db = 10 * log10(1000*rayTx.TransmitterPower); % Convert W to dBm
            Gtx_db = gain(rayTx, fq, aod(1), aod(2)); % Transmit gain
            Grx_db = gain(rayRx, fq, aoa(1), aoa(2)); % Receive gain
            v = Ptx_db + Gtx_db + Grx_db - pl - rayTx.SystemLoss - rayRx.SystemLoss;
        end
        
        % Discard ray if it does not meet threshold
        if isPathlossPlot
            % Check if path loss is greater than max value
            rayDoesNotMeetThreshold = v > clim(2);
        else
            % Check if power is less than min value
            rayDoesNotMeetThreshold = v < clim(1);
        end
        if rayDoesNotMeetThreshold
            continue
        else
            noPathsMeetThreshold = false;
        end
        
        % Compute color value for ray
        rayColorRGB = rfprop.internal.ColorUtils.colorcode(v, cmap, clim);
        rayColor = rfprop.internal.ColorUtils.rgb2css(rayColorRGB);
        
        % Build infobox description
        powDecPlaces = rfprop.Constants.MaxPowerInfoDecimalPlaces;
        distDecPlaces = rfprop.Constants.MaxDistanceInfoDecimalPlaces;
        phaseDecPlaces = rfprop.Constants.MaxPhaseInfoDecimalPlaces;
        angDecPlaces = rfprop.Constants.MaxAngleInfoDecimalPlaces;
        if isPathlossPlot
            valueDescription = message('shared_channel:rfprop:RaytraceDescriptionPathloss', mat2str(round(v,powDecPlaces))).getString;
        else
            valueDescription = message('shared_channel:rfprop:RaytraceDescriptionPower', mat2str(round(v,powDecPlaces))).getString;
        end
        
        % Get number of reflections and diffractions 
        numRef = (1-ray.LineOfSight)*sum('R' == [ray.Interactions.Type], 2);
        numDiff = ray.NumInteractions - numRef; 
        
        % Formulate the interaction type string if the ray has both
        % reflections and diffractions
        if (numRef > 0) && (numDiff > 0)
            intactType = cellfun(@(x)x(1), {ray.Interactions.Type});
            intactStr = [intactType; repmat(',', 1, length(intactType))];
            intactStr = (intactStr(:))';
            intactStr(end) = [];
            intactDescription = [message('shared_channel:rfprop:RaytraceDescriptionInteractions', intactStr).getString, '<br>'];
        else
            intactDescription = ''; 
        end
        
        rayDesc = [...
            message('shared_channel:rfprop:RaytraceDescriptionNumReflections', mat2str(numRef)).getString, '<br>', ...
            message('shared_channel:rfprop:RaytraceDescriptionNumDiffractions', mat2str(numDiff)).getString, '<br>', ...
            intactDescription, ...
            valueDescription, '<br>', ...
            message('shared_channel:rfprop:RaytraceDescriptionPhaseChange', mat2str(round(ray.PhaseShift,phaseDecPlaces))).getString, '<br>', ...
            message('shared_channel:rfprop:RaytraceDescriptionDistance', mat2str(round(pathDistance,distDecPlaces))).getString, '<br>', ...
            message('shared_channel:rfprop:RaytraceDescriptionAoD', mat2str(round(aod(1),angDecPlaces)), mat2str(round(aod(2),angDecPlaces))).getString, '<br>', ...
            message('shared_channel:rfprop:RaytraceDescriptionAoA', mat2str(round(aoa(1),angDecPlaces)), mat2str(round(aoa(2),angDecPlaces))).getString];
        
        % Fill the ray infobox with the materials if the ray had any
        % interactions
        if ray.NumInteractions > 0
            materialDescription = "Materials: ";
            for interactionIndex = 1:ray.NumInteractions
                materialDescription = materialDescription + ray.Interactions(interactionIndex).MaterialName;
                if interactionIndex < ray.NumInteractions
                    materialDescription = materialDescription + ", ";
                end
            end
            rayDesc = [rayDesc, '<br>', char(materialDescription)]; %#ok<AGROW>
        end

        % Generate viewer-specific ID
        rayID = viewer.getId(1);
        rayID = ['ray' num2str(rayID{1})];
        
        % Create IDs for each raytrace line.
        % Number of interactions dictates number of Lines (and thus, number
        % of lineIDs necessary)
        rayIDToUse = cell(1 + numIntact, 1);
        rayIDToUse{1} = rayID;
        for i = 2:numIntact + 1
            % Line IDs are generated to use the unique rayID + the
            % bounce number
            rayIDToUse{i} = [rayID 'bounce' num2str(i - 1)];
        end
        
        if associateWithSites
            % Associate the two site with each other through the viewer's SiteGraphics
            viewer.associateSiteGraphics(rayTx.UID, rayRx.UID, 'rays', rayIDToUse);
        else
            % Store the plotted graphic ID in SiteGraphics and assocaite it
            % with this ray object using ray.UID
            viewer.setGraphic(ray.UID, rayIDToUse, 'rays');
        end
        
        % Add the ray to the list of items dependent on a legend
        if (showLegend)
            viewer.addToLegendEntities(rayIDToUse);
        end
        
        % Initialize ray vertex locations, which includes the tx and rx
        % locations along with locations of any interactions
        numLocations = 2 + numIntact;
        rayLats = zeros(numLocations,1);
        rayLons = zeros(numLocations,1);
        rayElevs = zeros(numLocations,1);
        
        % Ray path starts with tx location
        rayLats(1) = ray.TransmitterLocation(1);
        rayLons(1) = ray.TransmitterLocation(2);
        rayElevs(1) = ray.TransmitterLocation(3);
        
        % Add interaction locations
        if ~ray.LineOfSight
            intactLocations = [ray.Interactions.Location];
            rayLats(2:numIntact+1) = intactLocations(1,:);
            rayLons(2:numIntact+1) = intactLocations(2,:);
            rayElevs(2:numIntact+1) = intactLocations(3,:);
        end
        
        % Ray path ends with rx location
        rayLats(end) = ray.ReceiverLocation(1);
        rayLons(end) = ray.ReceiverLocation(2);
        rayElevs(end) = ray.ReceiverLocation(3);

        % Add line graphic for ray
        if (associateWithSites)
            linkedGraphics = {rayTx.UID; rayRx.UID};
        else
            linkedGraphics = {};
        end
        
        if viewer.IsCartesian
            args = {"Color", reshape(rayColorRGB, 1, 3), ...
                "Name", message('shared_channel:rfprop:RaytracePropagationPathName').getString, ...
                "Description", rayDesc, "ID", rayIDToUse};
        else
            args = {"Color", {rayColor}, "Description", rayDesc, "ID", rayIDToUse, ...
                "Name", message('shared_channel:rfprop:RaytracePropagationPathName').getString, ...
                "Width", 2, "LinkedGraphics", {linkedGraphics}};
        end
        viewer.Visualizer.line([rayLats rayLons rayElevs], args{:});
    end
    
    % Remove all stale graphics
    viewer.remove(graphicsToRemove);
    
    % Warn if no paths to plot and return early
    if noPathsMeetThreshold
        if isPathlossPlot
            warning(message('shared_channel:rfprop:NoPathsMeetRequiredPathLoss', mat2str(clim(2))));
        else
            warning(message('shared_channel:rfprop:NoPathsMeetRequiredPower', mat2str(clim(1))));
        end
        viewer.Visualizer.submitPlots('Animation', 'none');
        return
    end
    
    % Get legend data
    if showLegend
        [legendColors, legendColorValues] = rfprop.internal.ColorUtils.colormaplegend(cmap, clim);
    else
        legendColors = "";
        legendColorValues = "";
    end
    if isPathlossPlot
        legendColors = fliplr(legendColors);
        legendColorValues = fliplr(legendColorValues);
    end
    if (showLegend)
        legendInd = viewer.getId(1);
        legendID = ['raylegend', num2str(legendInd{1})];
        viewer.Visualizer.legend(legendTitle, legendColors, legendColorValues, "ID", legendID);
        viewer.LegendID = legendID;
    end
    
    if isViewerInitiallyVisible && ~viewer.Visible
        viewer.Visualizer.submitPlots('Animation', 'none');
        return % Abort if Site Viewer has been closed
    end
    
    submitAnim = 'fly';
    % Zoom if the viewer isn't open yet
    if (~isViewerInitiallyVisible && (isempty(animation) || strcmp(animation, 'zoom')))
        submitAnim = 'zoom';
    end
    viewer.Visualizer.submitPlots('Animation', submitAnim);
    % Show the viewer and bring it into focus
    if (~isViewerInitiallyVisible && enableWindowLaunch)
        viewer.Visible = true;
    end
    viewer.bringToFront;
  end
end

methods (Hidden) % Hidden public methods
  function validateConfig(obj)
    for i = 1:numel(obj)
        % Validate location range for geographic 
        coder.internal.errorIf( ...
            ~isInactiveProperty(obj(i), 'TransmitterLocation') &&  ...
            strcmp(obj(i).CoordinateSystem, 'Geographic') && ...
            (obj(i).TransmitterLocation(1) <-90 || ...
            obj(i).TransmitterLocation(1) > 90), ...
            'shared_channel:Ray:InvalidGeoLocation');
    
        coder.internal.errorIf( ...
            ~isInactiveProperty(obj(i), 'ReceiverLocation') &&  ...
            strcmp(obj(i).CoordinateSystem, 'Geographic') && ...
            (obj(i).ReceiverLocation(1) <-90 || ...
            obj(i).ReceiverLocation(1) > 90), ...
            'shared_channel:Ray:InvalidGeoLocation');
        
        if ~isInactiveProperty(obj(i), 'Interactions') &&  ...
            strcmp(obj(i).CoordinateSystem, 'Geographic') 
            intactLoc = [obj(i).Interactions.Location]; 
            coder.internal.errorIf( ...
                (any(intactLoc(1,:) <-90) || any(intactLoc(1,:) > 90)), ...
                'shared_channel:Ray:InvalidGeoLocation');
        end
        
        path = convertPathToCartesian(obj(i), obj(i).propagationPath);
        
        % Check Tx and Rx positions cannot be the same
        coder.internal.errorIf( ...
            norm(path(:,1) - path(:,end)) <= sqrt(eps), ...
            'shared_channel:Ray:SameTxRxPos')

        if obj(i).NumInteractions
            % Check consecutive points on the path cannot be the same
            coder.internal.errorIf(any(sqrt(sum(diff(path, 1, 2).^2)) < sqrt(eps)), ...
                'shared_channel:Ray:SameConsecutivePts');
            
            % If the interaction is reflection, two consecutive segments cannot have the same direction
            pathDir = diff(path, 1, 2); % [3, NI+1]
            pathDir = pathDir ./ sqrt(sum(pathDir.^2, 1));

            for n = 1:obj(i).NumInteractions
                if strcmpi(obj(i).Interactions(n).Type,'Reflection')
                    coder.internal.errorIf(sqrt(sum((pathDir(:,n) - pathDir(:,n+1)).^2)) <= sqrt(eps), ...
                        'shared_channel:Ray:SameConsecutiveDir');
                end
            end
        end
    end
  end
 
  function [polMtx, propMtx] = getPolMatrix(ray, refMtl, txAxes, rxAxes)
    %GETPOLMATRIX Return polarization matrix and propagator matrix
    %   [POLMTX, PROPMTX] = getPolarizationMatrix(RAY, REFMTL, TXAXES,
    %   RXAXES, ISPOL) returns the polarization matrix, POLMTX, in the form
    %   of [HH HV; VH VV] for a ray, RAY, due to reflections, given the
    %   interaction materials, transmit array orientation axes, and receive
    %   array orientation axes. POLMTX is a 2-by-2 matrix. PROPMTX is the
    %   propagator matrix for an electric field i.e. Erx = propMtx*Etx; it
    %   is a 3x3 matrix that operates on an electric field with Cartesian
    %   components as [Ex;Ey;Ez]. exp(-iwt) convention.
    %
    % A note on polarization terminology (relative to local reflecting plane):
    % General       | Electromagnetics          | Physics, Optics | Acoustics
    % -----------------------------------------------------------------------
    % Parallel      | TM, vertical (V), H-Pol   | p-like          | hard
    % -----------------------------------------------------------------------
    % Perpendicular | TE, horizontal (H), E-Pol | s-like          | soft
    % -----------------------------------------------------------------------
    %
    % REFERENCES:
    %   [1] G. Yun, "Polarization Ray Tracing", University of Arizona
    %   dissertation, 2011.  Note that Yun may use different conventions
    %   than common electromagnetics literature.

    % Rotate the receiver's axes for a geographic system
    if strcmp(ray.CoordinateSystem, 'Geographic')
        % For a geographic system, we need an extra rotation from Rx
        % ENU into Tx ENU. Multiplying rxAxes with it gives the "true"
        % rotation from LCS to GCS (which is Tx ENU) for Rx.
        enuRot = enuRotation(ray.TransmitterLocation, ray.ReceiverLocation);
        rxAxes = enuRot{1} * rxAxes;
    end

    % Calculate polarization and propagation matrices. If the ray is direct
    % (line-of-sight), the math is simpler than for indirect ray.
    propMtx = complex(eye(3));
    if ray.LineOfSight
        % If the transmitter's and receiver's axes are the same, then there
        % is no cross-coupling. Otherwise, perform path vector math.
        if isequal(txAxes, rxAxes)
            polMtx = [1 0; 0 1];
        else
            % Determine the ray path unit vectors. The ray path is found by
            % taking end point - start point. Working in Cartesian.
            % Normalize to get unit vector.
            path = convertPathToCartesian(ray, ray.propagationPath); % [3,2]
            pathDir = diff(path, 1, 2); % [3,1]
            pathDir = pathDir/sqrt(sum(pathDir.^2, 1));

            % Get V-pol direction of the TX's E-field.  This uses the fact
            % that the V-pol starts as the ray direction rotated -theta 90
            % degrees.
            incVPolDir = getVPolDir(pathDir, txAxes); % [x;y;z]

            % Get vertical pol. direction of the E field at Rx, using the
            % last segment
            rxVPolDir = getVPolDir(pathDir, rxAxes); % [x;y;z]
            rxHPolDir = cross(rxVPolDir, pathDir);   % Yun convention

            % Project propagated fields against expected receiver polarizations
            propHPol = cross(incVPolDir, pathDir); % Yun convention
            polMtx = [dot(propHPol, rxHPolDir) dot(incVPolDir, rxHPolDir); ...
                dot(propHPol, rxVPolDir) dot(incVPolDir, rxVPolDir)];
        end
    else % Indirect (reflected)
        % Calculate Fresnel reflection coefficients & complex relative perm.
        refCoeff = getReflectance(ray, refMtl); % [H;V]

        % Determine the ray path unit vectors. The ray path is found by
        % taking end point - start point.  Working in Cartesian. Normalize
        % to get unit vector.
        path = convertPathToCartesian(ray, ray.propagationPath);
        pathDir = diff(path, 1, 2);
        pathNorm = sqrt(sum(pathDir.^2, 1));
        pathDir = pathDir./pathNorm;

        % Get V-pol direction of the TX's E-field.  This uses the fact that
        % the V-pol starts as the ray direction rotated -theta 90 degrees.
        incPath = pathDir(:,1);
        incVPolDir = getVPolDir(incPath, txAxes); % [x;y;z]

        % Get pol. directions of the E field at Rx, using the last segment
        rxVPolDir = getVPolDir(pathDir(:,end), rxAxes); % [x;y;z]
        rxHPolDir = cross(rxVPolDir, pathDir(:,end));   % Yun convention

        % Iterate through reflection interactions to derive polarization
        % matrix
        outPath = incPath;
        refCoeffIdx = 1;
        for intIdx = 1:ray.NumInteractions
            % Determine polarization coupling matrix.  This is a repeated
            % product for each interaction.
            incPath = outPath;
            outPath = pathDir(:, intIdx+1);
            % Rotation matrix to local surface for H/V pol.
            % [Yun 2011, (2.2.3)-(2.2.4)]
            if sum(abs(incPath+outPath)) < eps % Normal incidence (outPath = -incPath)
                s_in = cross(incPath, outPath+eps);
            else
                s_in = cross(incPath, outPath);
            end
            s_in = s_in / norm(s_in);
            p_in = cross(incPath, s_in);
            rotMtxIn = [s_in p_in incPath];

            % Rotation matrix back to global coordinates
            % [Yun 2011, (2.2.3)-(2.2.4)]
            p_out = cross(outPath, s_in);
            rotMtxOut = [s_in p_out outPath];

            % Cascade reflection into polarization matrix
            % [Yun 2011, (2.2.6)-(2.2.7), pg. 68]
            propMtx = rotMtxOut * diag([refCoeff(:,refCoeffIdx)' 1]) ...
                * rotMtxIn' * propMtx;

            % Increment reflection coefficient index
            refCoeffIdx = refCoeffIdx + 1;
        end

        % Calculate propagated V & H fields in 3D Cartesian
        propVPol = propMtx * incVPolDir;
        propHPol = propMtx * cross(incVPolDir, pathDir(:,1)); % Yun convention

        % Project propagated fields against expected receiver polarizations
        polMtx = [dot(propHPol, rxHPolDir) dot(propVPol, rxHPolDir); ...
            dot(propHPol, rxVPolDir) dot(propVPol, rxVPolDir)];
    end
  end
end

methods (Access = protected)
  function flag = isInactiveProperty(obj, prop)
    if strcmp(prop, 'SystemScale')
        flag = strcmp(obj.PathSpecification, 'Delay and angles') || ...
            strcmp(obj.CoordinateSystem, 'Geographic');
    elseif any(strcmp(prop, {'CoordinateSystem','SystemScale', ...
        'TransmitterLocation','ReceiverLocation', ...
        'LineOfSight','NumInteractions'}))
        flag = strcmp(obj.PathSpecification, 'Delay and angles');
    elseif any(strcmp(prop, {'Interactions'}))
        flag = strcmp(obj.PathSpecification, 'Delay and angles') || ...
            obj.LineOfSight;
    else
        flag = false;
    end
  end
  
  function groups = getPropertyGroups(obj)
    % Returns two property groups for display: normal and info
    
    % Get all properties and their values
    propName  = properties(obj);

    % Initialization
    propVal = cell(numel(propName), 1);    
    activeIdx = true(size(propName));
    infoIdx = false(size(propName));
    
    % For each property, determine its activeness and read-only status
    for n = 1:numel(propName)
        if isInactiveProperty(obj, propName{n})
            % If it is inactive then do not add it to the property list
            activeIdx(n) = false;
        elseif isUndefinedProperty(obj, propName{n})
            % If it is undefined then set the value to an empty for disp()
            propVal{n} = [];
        else
            propVal{n} = obj.(propName{n});
        end

        % Check read-only or not
        infoIdx(n) = isPropertyReadOnly(obj, propName{n});
    end
    
    % Create two property lists, one of read-only properties and one of
    % normal properties
    normalPropList = cell2struct(propVal(activeIdx&~infoIdx), ...
        propName(activeIdx&~infoIdx));
    readOnlyPropList = cell2struct(propVal(activeIdx&infoIdx),  ...
        propName(activeIdx&infoIdx));
    
    groups = [ ...
        matlab.mixin.util.PropertyGroup(normalPropList) ...
        matlab.mixin.util.PropertyGroup(readOnlyPropList, ...
        'Read-only properties:')];
  end
end

methods (Access = private)
  function performReadOnlyCheck(obj, propName)
    coder.internal.errorIf(isPropertyReadOnly(obj, propName), ...
        'shared_channel:Ray:NoPropertyChangeWhenReadOnly', propName);
  end

  function flag = isPropertyReadOnly(obj, propName)
    if any(strcmp(propName, {'NumInteractions', 'PropagationDelay', ...
        'AngleOfDeparture','AngleOfArrival'}))
        flag = strcmp(obj.PathSpecification, 'Locations');
    elseif any(strcmp(propName, {'PathLoss', 'PhaseShift'}))
        flag = strcmp(obj.PathLossSource, 'Free space model');
    elseif strcmp(propName, 'PropagationDistance') 
        flag = true; 
    else
        flag = false;
    end
  end
    
  function path = propagationPath(obj)
    % Calculate propgation path from tx to rx in one coordinate system
    
    if obj.LineOfSight    
        path = [obj.TransmitterLocation, obj.ReceiverLocation];
    else
        path = [obj.TransmitterLocation, ...
            [obj.Interactions.Location], obj.ReceiverLocation];
    end
  end
  
  function cartPath = convertPathToCartesian(obj, path)
    % For a geographic system, convert every location on the path into the
    % ENU position of the transmitter. No conversion is needed for a
    % Cartesian system.
    if strcmp(obj.CoordinateSystem, 'Geographic')
        [x, y, z] = geodetic2enu(...
            path(1,1),     path(2,1),     path(3,1), ...
            path(1,2:end), path(2,2:end), path(3,2:end));
        cartPath = [0 x; 0 y; 0 z];
    else        
        cartPath = path * obj.SystemScale;
    end
  end  
    
  function ang = localangle(obj, location1, location2)
    % Return [az; el] angle of location 2 relative to location 1 in the ENU
    % system of location 1. 
      
    if strcmp(obj.CoordinateSystem, 'Geographic')
        [X,Y,Z] = geodetic2enu(...
            location1(1), location1(2), location1(3), ...
            location2(1), location2(2), location2(3));
    else
        X = location2(1) - location1(1);
        Y = location2(2) - location1(2);
        Z = location2(3) - location1(3);
    end
       
    % Convert Cartesian position to az/el angles in radians.
    % az is in (-pi, pi] and el is in [-pi/2 pi/2]. 
    [azrad, elrad] = cart2sph(X,Y,Z);
            
    % Convert to degrees
    az = rad2deg(azrad);
    el = rad2deg(elrad);
    
    ang = [az; el];
  end

  function refCoeff = getReflectance(ray, refMtl)
    %GETREFLECTANCE Return reflection coefficient for each material
    %   [REFCOEFF, EPSILON] = getReflectance(RAY, REFMTL) returns the
    %   reflection coefficient, REFCOEFF, in linear scale for horizontal &
    %   vertical polarizations at each reflection for a non-line-of-sight
    %   ray, RAY, given the reflection materials, REFMTL. REFCOEFF is a
    %   2-by-X matrix with the first and second rows representing
    %   horizontal and vertical polarizations respectively, where X is the
    %   number of reflections. exp(-iwt); vertical polarization coefficient
    %   = -Eref/Einc conventions; permeability assumed that of free-space;
    %   ray assumed to propagate in free-space.
    
    if ray.LineOfSight
        refCoeff = complex(ones(2,0));
    else
        fc = ray.Frequency;
        rayAoI = ray.AngleOfIncidence;

        % Calculate number of materials and determine which interactions &
        % materials are for diffractions
        numMtl = ray.NumInteractions;

        % Get complex relative permittivity and angle of incidence
        scale = 1i/(2*pi*fc*8.854187817e-12); % 1i/(2*pi*f*epsilon0)
        epsilon = complex(nan(1,numMtl));
        aoi = ones(1,numMtl);
        for mtlIdx = 1:numMtl
            % Material could be cell aray, char, or numeric array
            if iscell(refMtl)
                thisMtl = refMtl{mtlIdx};
            else
                thisMtl = refMtl(:,mtlIdx);
            end

            % Material value could be NaN, char (material name), or
            % [relative real permittivity; conductivity (S/m)]
            if ischar(thisMtl)
                [~, ~, epsilon(mtlIdx)] = ...
                    comm.internal.factoryMaterialPermittivity(thisMtl, fc); % exp(jwt) convention
                epsilon(mtlIdx) = conj(epsilon(mtlIdx)); % exp(-iwt) convention
                aoi(mtlIdx) = rayAoI(mtlIdx);
            else
                if ~isnan(thisMtl)
                    epsilon(mtlIdx) = thisMtl(1) + thisMtl(2)*scale; % exp(-iwt) convention
                end
                aoi(mtlIdx) = rayAoI(mtlIdx);
            end
        end

        % Calculate reflection coefficients [H;V]. Assume perfect
        % electrical conductor (PEC) ([-1;1])
        refCoeff = repmat(complex([-1;1]), 1, numMtl);

        % For finite permittivity values and non-diffraction interactions,
        % calculate using Fresnel plane-wave reflection formula
        useFresnel = isfinite(epsilon);
        if any(useFresnel)
            AoI = aoi(useFresnel);
            N2 = epsilon(useFresnel); % = permittivity{reflecting medium}/permittivity{incident ray medium = free-space}
            B_H = sqrt(N2 - sind(AoI).^2);
            B_V = B_H./N2;
            cosdAoI = cosd(AoI);
            refCoeff(:, useFresnel) = [ ...
                (cosdAoI - B_H)./ ...
                (cosdAoI + B_H); ... % H
                (cosdAoI - B_V)./ ...
                (cosdAoI + B_V)]; % V
        end
    end
  end

end

methods (Static, Hidden)
    function validateInteractionsStruct(val)
        % Validate Interactions struct

        numFieldsIntact = numel(fields(val));
        if numFieldsIntact == 2
            % Allow rays from previous releases with only 2 fields in the
            % interactions property
            coder.internal.errorIf(~all(isfield(val, {'Type', 'Location'})), ...
                'shared_channel:Ray:InvalidInteractionStruct');
        else
            coder.internal.errorIf((numFieldsIntact ~= 3) || ...
                ~all(isfield(val, {'Type', 'Location', 'MaterialName'})), ...
                'shared_channel:Ray:InvalidInteractionStruct');
        end
        for i = 1:length(val)
            type = val(i).Type;
            coder.internal.errorIf( ...
                (~ischar(type) && ~isStringScalar(type)) || ...
                ~any(strcmp(type, {'Reflection', 'Diffraction'})), ...
                'shared_channel:Ray:InvalidInteractionType', i)
            validateattributes(val(i).Location, {'double'}, ...
                {'real','finite','size',[3 1]}, ...
                'comm.Ray.Interactions.Location', ...
                'Interactions.Location');
        end
    end
end
end

function validateAngleRange(val)
% Validate angle range [azimuth; elevation] in degrees

coder.internal.errorIf( ...
        (val(1) <= -180) || (val(1) > 180) || (abs(val(2)) > 90), ...
        'shared_channel:Ray:InvalidAzEl');

end

function sites = validateSites(type, sites, numRays)

% Sites can be scalar or same-sized as number of rays
if isscalar(sites)
    numSites = 1;
else
    numSites = numRays;
end

validateattributes(sites,{type},{'nonempty','vector','numel',numSites},'plot');

% Perform scalar expansion
if isscalar(sites)
    sites = repmat(sites,numRays,1);
end
end

function validateRaysForPlot(rays)
numRays = numel(rays);
for i = 1:numRays
    coder.internal.errorIf(...
        ~strcmp(rays(i).PathSpecification, "Locations"), ...
        'shared_channel:Ray:InvalidConfigForPlot');
end
end

function VPolDirGCS = getVPolDir(rayDirGCS, LCS2GCSAxes)

% If 'LCS2GCSAxes' is eye(3), then local = global coordinates. Otherwise,
% perform coordinate transformation.
if isequal(LCS2GCSAxes, eye(3))
    % Get ray spherical direction in LCS
    [az, el] = cart2sph(rayDirGCS(1), rayDirGCS(2), rayDirGCS(3));
    % Get V-pol Cartesian direction in LCS
    [x, y, z] = sph2cart(az, el+pi/2, 1);
    VPolDirGCS = [x;y;z];
else
    % Get ray Cartesian direction in LCS
    rayDirLCS = global2localcoord(rayDirGCS, 'rr', zeros(3,1), LCS2GCSAxes);
    % Get ray spherical direction in LCS
    [az, el] = cart2sph(rayDirLCS(1), rayDirLCS(2), rayDirLCS(3));
    % Get V-pol Cartesian direction in LCS
    [x, y, z] = sph2cart(az, el+pi/2, 1);
    % Get V-pol Cartesian direction in GCS
    VPolDirGCS = local2globalcoord([x; y; z], 'rr', zeros(3,1), LCS2GCSAxes);
end

end

function [x,y,z] = geodetic2enu(lat0, lon0, h0, lat, lon, h)

if isempty(coder.target)
    [x, y, z] = rfprop.internal.MapUtils.geodetic2enu(...
        lat0, lon0, h0, lat, lon, h);
else
    a = 6378137;
    f = 1/298.257223563;
    
    s1 = sind(lat0);
    c1 = cosd(lat0);
    
    s2 = sind(lat);
    c2 = cosd(lat);
    
    p1 = c1 .* cosd(lon0);
    p2 = c2 .* cosd(lon);
    
    q1 = c1 .* sind(lon0);
    q2 = c2 .* sind(lon);
    
    e2 = f * (2 - f);
    w1 = 1 ./ sqrt(1 - e2 * s1.^2);
    w2 = 1 ./ sqrt(1 - e2 * s2.^2);
    deltaX =            a * (p2 .* w2 - p1 .* w1) + (h .* p2 - h0 .* p1);
    deltaY =            a * (q2 .* w2 - q1 .* w1) + (h .* q2 - h0 .* q1);
    deltaZ = (1 - e2) * a * (s2 .* w2 - s1 .* w1) + (h .* s2 - h0 .* s1);
    
    cosPhi = cosd(lat0);
    sinPhi = sind(lat0);
    cosLambda = cosd(lon0);
    sinLambda = sind(lon0);
    
    t =  cosLambda .* deltaX + sinLambda .* deltaY;
    x = -sinLambda .* deltaX + cosLambda .* deltaY;    
    z =  cosPhi .* t + sinPhi .* deltaZ;
    y = -sinPhi .* t + cosPhi .* deltaZ;
end

end

function [lat,lon,h] = enu2geodetic(lat0, lon0, h0, pos)

xEast = pos(1,:)';
yNorth = pos(2,:)';
zUp = pos(3,:)';

a = 6378137;
f = 1/298.257223563;

sinphi = sind(lat0);
cosphi = cosd(lat0);

e2 = f * (2 - f);
N  = a ./ sqrt(1 - e2 * sinphi.^2);
rho = (N + h0) .* cosphi;
z0 = (N*(1 - e2) + h0) .* sinphi;
x0 = rho .* cosd(lon0);
y0 = rho .* sind(lon0);

cosPhi = cosd(lat0);
sinPhi = sind(lat0);
cosLambda = cosd(lon0);
sinLambda = sind(lon0);

t = cosPhi .* zUp - sinPhi .* yNorth;
dz = sinPhi .* zUp + cosPhi .* yNorth;
dx = cosLambda .* t - sinLambda .* xEast;
dy = sinLambda .* t + cosLambda .* xEast;

x = x0 + dx;
y = y0 + dy;
z = z0 + dz;

rho = hypot(x,y);
[lat, h] = cylindrical2geodetic(rho, z, a, f);
lon = atan2d(y,x);
   
end

function [phi, h] = cylindrical2geodetic(rho, z, a, f)

b = (1 - f) * a;
e2 = f * (2 - f);
ep2 = e2 / (1 - e2);

beta = atan2d(z, (1 - f) * rho);
phi = atan2d(z + b * ep2 * sind(beta).^3,...
    rho - a * e2 * cosd(beta).^3);

betaNew = atan2d((1 - f)*sind(phi), cosd(phi));
count = 0;
while any(beta(:) ~= betaNew(:)) && count < 5
    beta = betaNew;
    phi = atan2d(z + b * ep2 * sind(beta).^3,...
        rho - a * e2 * cosd(beta).^3);
    betaNew = atan2d((1 - f)*sind(phi), cosd(phi));
    count = count + 1;
end

sinphi = sind(phi);
N = a ./ sqrt(1 - e2 * sinphi.^2);
h = rho .* cosd(phi) + (z + e2 * N .* sinphi) .* sinphi - N;

end

function enuRot = enuRotation(txLoc, rxLoc)

if isempty(coder.target)
    enuRot = rfprop.internal.MapUtils.enuRotation( ...
                txLoc(1), txLoc(2), txLoc(3), ...
                rxLoc(1), rxLoc(2), rxLoc(3));
else
    % Calculate lat/lon/height of axes points ([1;0;0],[0;1;0],[0;0,1]) in
    % ENU of the 2nd point
    [latXYZ, lonXYZ, heightXYZ] = enu2geodetic( ...
        rxLoc(1), rxLoc(2), rxLoc(3), eye(3));

    % Calculate [X;Y;Z] of the 2nd point in the ENU of the 1st point
    [x0, y0, z0] = geodetic2enu( ...
        txLoc(1), txLoc(2), txLoc(3), rxLoc(1), rxLoc(2), rxLoc(3));

    % Calculate [X;Y;Z] of axes points in the ENU of 1st point
    [x1, y1, z1] = geodetic2enu( ...
        txLoc(1), txLoc(2), txLoc(3), latXYZ, lonXYZ, heightXYZ);

    % The difference between the 2nd point's axes and the 2nd point's
    % origin in the ENU of the 1st point gives the rotation matrix.
    enuRot = {[x1'-x0'; y1'-y0'; z1'-z0']};
end

end

function rotMtx = getRotationMatrix(VPolDir1, HPolDir1, VPolDir2)

% Get sin and cos of the rotation angles for this
% reflection. It is the angle between the V-pol (s-pol)
% directions of the incidence ray and the reflected ray.
cosAng = dot(VPolDir1, VPolDir2);

% Precision issue, rounded off to 0 or 1 so the path loss can be "true" inf
cosAng(abs(cosAng - 1) < sqrt(eps)) = 1;
cosAng(abs(cosAng) < sqrt(eps)) = 0;

% Get H-pol direction of the E field for this interaction. It is the cross
% product of ray direction and V-pol direction.
% HPolDir2 = cross(rayDir2, VPolDir2);

% Dot product between H-Pol and next segment's V-pol to decide the sign of
% sin angle
angSign = dot(HPolDir1, VPolDir2);     % [1 NI+1]

% Positive sin if the angSign is positive and negetive sin otherwise
sinAng = (2*(angSign>0)-1) .* sqrt(1-cosAng^2);        % [1 NI+1]

% This following matrix addresses geometric coupling.
rotMtx = [cosAng sinAng; ...
         -sinAng cosAng];

end