function varargout = link(rxs, txs, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.

% Validate sites
validateattributes(rxs,{'rxsite'},{'nonempty'},'link','',1);
validateattributes(txs,{'txsite'},{'nonempty'},'link','',2);
usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(rxs, txs);

% Process optional inputs
p = inputParser;
if nargin > 2 && mod(numel(varargin),2)  % Odd number of inputs
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Animation', '');
p.addParameter('EnableWindowLaunch', true);
p.addParameter('SuccessColor', 'green');
p.addParameter('FailColor', 'red');
p.addParameter('Map', []);
p.parse(varargin{:});

% Get Site Viewer visibility state and validate web graphics
outputRequested = nargout > 0;
if ~outputRequested
    viewer = rfprop.internal.Validators.validateMap(p, 'link', usingCartesian);
    isViewerInitiallyVisible = viewer.Visible;
    if usingCartesian
        [map, ~, ~, mapStruct] = rfprop.internal.Validators.validateCartesianMap(p);
        pm = rfprop.internal.Validators.validateCartesianPropagationModel(p, map, 'link');
        if isprop(pm, 'CoordinateSystem')
            pm.CoordinateSystem = 'cartesian';
        end
    else
        [map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'link');
        rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'geographic');
        pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'link');
    end
else
    isViewerInitiallyVisible = false;
    % Validate and get parameters
    if usingCartesian
        [map, ~, ~, mapStruct] = rfprop.internal.Validators.validateCartesianMap(p);
        rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'cartesian');
        pm = rfprop.internal.Validators.validateCartesianPropagationModel(p, map, 'link');
    else
        [map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'link');
        rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'geographic');
        pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'link');
    end
end
% Validate and get parameters
rfprop.internal.Validators.validateMaxNumReflections(pm, 'link');
[animation, enableWindowLaunch] = rfprop.internal.Validators.validateGraphicsControls(p, isViewerInitiallyVisible, 'link');
successColor = rfprop.internal.Validators.validateLineColor(p, 'SuccessColor', 'link');
failColor = rfprop.internal.Validators.validateLineColor(p, 'FailColor', 'link');

% Get site antenna coordinates
txsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(txs, map);
rxsCoords = rfprop.internal.AntennaSiteCoordinates.createFromAntennaSites(rxs, map);

% Get distances (for infobox) and signal strengths
ds = distance(txs, rxs, 'Map', map,...
    'SourceAntennaSiteCoordinates', txsCoords, ...
    'TargetAntennaSiteCoordinates', rxsCoords);

ss = sigstrength(rxs, txs, pm, 'Type', 'power', ...
    'Map', mapStruct, ...
    'TransmitterAntennaSiteCoordinates', txsCoords, ...
    'ReceiverAntennaSiteCoordinates', rxsCoords);

% Allocate output matrix
numTx = numel(txs);
numRx = numel(rxs);
status = zeros(numTx, numRx, 'logical');

% Compute link status from each transmitter to all receivers
numLinks = numTx * numRx;
startIDs = cell(1,numLinks);
endIDs = startIDs;
startPositions = startIDs;
endPostions = startIDs;
lineColors = startIDs;
lineInfo = startIDs;
graphicsToRemove = {};
if ~outputRequested
    viewer.Visualizer.queuePlots;
    % Show the sites at the start of link to initialize siteGraphics and
    % ensure that stale graphics get cleared before link is plotted.
    show(rxs,'Map',viewer,'Animation','none','EnableWindowLaunch', false, ...
        'AntennaSiteCoordinates',rxsCoords);
    show(txs,'Map',viewer,'Animation','none','EnableWindowLaunch', false, ...
        'AntennaSiteCoordinates',txsCoords);
end
k = 0;
for txInd = 1:numTx
    tx = txs(txInd);
    txCoord = txsCoords.extract(txInd);
    
    for rxInd = 1:numRx
        rx = rxs(rxInd);
        rxCoord = rxsCoords.extract(rxInd);
        if ~outputRequested
            viewer.initializeSiteGraphics(tx.UID);
            viewer.initializeSiteGraphics(rx.UID);
            linkID = viewer.getId(1);
            linkID = ['link' num2str(linkID{1})];
            % Remove existing link lines
            txLinkGraphics = viewer.getGraphic(tx.UID, 'link');
            if isfield(txLinkGraphics, rx.UID)
                graphicsToRemove = [graphicsToRemove; txLinkGraphics.(rx.UID)];
            end
            viewer.disassociateSiteGraphics(tx.UID, rx.UID, 'link');
        else
            linkID = 'dummyID';
        end
       % Compute link status
        Prx = ss(txInd, rxInd);
        sens = rx.ReceiverSensitivity;
        linkIsSuccessful = Prx >= sens;
        
        % Populate output or line data for plot
        if outputRequested
            status(txInd,rxInd) = linkIsSuccessful;
        else
            k = k+1;
            if usingCartesian
                startPositions{k} = txCoord.AntennaPosition;
                endPostions{k} = rxCoord.AntennaPosition;
            else
                startPositions{k} = [txCoord.DisplayLatitudeLongitude, txCoord.SurfaceHeightAboveTerrainReference];
                endPostions{k} = [rxCoord.DisplayLatitudeLongitude, rxCoord.SurfaceHeightAboveTerrainReference];
            end

            startIDs{k} = tx.UID;
            endIDs{k} = rx.UID;

            d = ds(rxInd,txInd)/1000; % Convert from m to km
            
            % Piece together description from parts so <br> can be
            % inserted, as HTML tags would break Message Catalog
            lineInfo{k} = [...
                message('shared_channel:rfprop:LinkDescriptionDistance', mat2str(round(d, rfprop.Constants.MaxDistanceInfoDecimalPlaces))).getString, '<br>', ...
                message('shared_channel:rfprop:LinkDescriptionPower', mat2str(round(Prx, rfprop.Constants.MaxPowerInfoDecimalPlaces))).getString, '<br>', ...
                message('shared_channel:rfprop:LinkDescriptionSensitivity', mat2str(round(sens, rfprop.Constants.MaxPowerInfoDecimalPlaces))).getString];
            if linkIsSuccessful
                if usingCartesian
                    lineColors{k} = successColor;
                else
                    lineColors{k} = rfprop.internal.ColorUtils.rgb2css(successColor);
                end
            else
                if usingCartesian
                    lineColors{k} = failColor;
                else
                    lineColors{k} = rfprop.internal.ColorUtils.rgb2css(failColor);
                end
            end
        end
        if (~outputRequested)
            % Remove stale graphics
            viewer.remove(graphicsToRemove);
            if usingCartesian
                lineNV = {"Color", lineColors{k}, "Width", 0.002, ...
                    "Description", ['<p style="text-align:center;">' message('shared_channel:rfprop:LinkName').getString '</p>' lineInfo{k}], ...
                    "ID", linkID};
                viewer.Visualizer.line([startPositions{k}; endPostions{k}],lineNV{:});
            else
                lineNV = {"Color", lineColors{k}, "Width", 2, "Description", lineInfo{k}, "FollowEllipsoid", true, ...
                    "ID", linkID, "Name", message('shared_channel:rfprop:LinkName').getString, ...
                    "LinkedGraphics", {{tx.UID; rx.UID}}};
                viewer.Visualizer.line([startPositions{k};...
                    endPostions{k}],lineNV{:});
            end
            viewer.associateSiteGraphics(tx.UID, rx.UID, 'link', linkID);
        end
    end
end

if outputRequested
    varargout = {status};
else
    % Show sites (but no animation and do not launch window)
    if isViewerInitiallyVisible && ~viewer.Visible
        return % Abort if Site Viewer has been closed (test before show)
    end
    if isViewerInitiallyVisible && ~viewer.Visible
        return % Abort if Site Viewer has been closed (test before line plot)
    end
    
    if (~isViewerInitiallyVisible && (isempty(animation) || strcmp(animation, 'zoom')))
        animation = 'zoom';
    end
    if isViewerInitiallyVisible
        animation = 'fly';
    end
    viewer.Visualizer.submitPlots("Animation", animation);
    if (enableWindowLaunch)
        if ~viewer.Visible
            viewer.Visible = true;
        end
        viewer.bringToFront;
    end
end
end