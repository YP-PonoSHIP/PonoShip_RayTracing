function [rxlats, rxlons, rxs, data] = radialReceiverLocationsLayoutData(txs, txsCoords, type, ...
    generatingMapMsg, computingDataMsg, pm, mapInput, res, radialResFactor, datarange, ...
    lonWest, lonEast, latSouth, latNorth, rxGain, rxAntennaHeight)
%radialReceiverLocationsLayoutData   Data for radial receiver locations layout
    
%   Copyright 2019-2024 The MathWorks, Inc.

% Check if pre-validated Map struct
if isstruct(mapInput)
    map = mapInput.Map;
else
    map = mapInput;
end

terrainSource = rfprop.internal.Validators.validateTerrainSource(map);

numTx = numel(txs);
isMultipleTx = numTx > 1;

% Compute resolutions and skip factor, which is number of terrain samples
% to advance by to get each receiver site location on the radial
targetradres = res/radialResFactor;
terrainres = min(targetradres, pm.terrainProfileResolution(map));
skipFactor = max(1,floor(targetradres/terrainres));

% For multiple tx, interpolation is required at common locations since
% the data for any particular tx is scattered on radials
txslatlon = txsCoords.LatitudeLongitude;
txslats = txslatlon(:,1);
txslons = txslatlon(:,2);
if isMultipleTx
    % Generate lat/lon data grid using resolution of receiver sites along
    % radials. A frame of data will be generated for each transmitter site
    % using interpolation. Use target resolution separating receiver sites
    % along radials for grid resolution.
    [lonsv, latsv, numGridPts] = rfprop.internal.MapUtils.geovec(targetradres, latNorth, latSouth, lonEast, lonWest);
    data = nan(numTx, numGridPts);
    [datalons,datalats] = meshgrid(lonsv,latsv);
    datalons = datalons(:);
    datalats = datalats(:);
    
    % Find data grid locations that are within range of any txsite
    gcdist = nan(numGridPts,numTx);
    for txInd = 1:numTx
        gcdist(:,txInd) = rfprop.internal.MapUtils.greatCircleDistance(...
            txslats(txInd), txslons(txInd), datalats, datalons);
    end
    isInDataRange = any(gcdist <= datarange,2);
end

% Convert grid resolution to angular resolution radials using chord
% math, where grid res is the chord length. Adjust to next lower value
% that uniformly divides 360 degrees.
txsnumradials = zeros(numTx,1);
for txInd = 1:numTx
    chordang = 2*asind(res/(2*datarange(txInd)));
    txsnumradials(txInd) = ceil(360/chordang);
end

% Launch status dialog mode based on total number of radials
if isa(map,'siteviewer') && ~isempty(generatingMapMsg)
    totalNumRadials = sum(txsnumradials);
    fullProgressDlg = totalNumRadials >= rfprop.Constants.MinNumRadialsFullProgressDialog;
    if fullProgressDlg
        % Use full, three-phase progress dialog
        msg = message('shared_channel:rfprop:ProgressDialogPreparingMapData').getString;
        launchBusyDialog(map, msg)
    else
        % Use single-phase progress dialog
        launchBusyDialog(map, generatingMapMsg)
    end
end

% Get signal strength data for each txsite
for txInd = 1:numTx
    tx = txs(txInd);
    txCoords = txsCoords.extract(txInd);
    txCoords.addCustomData('TxInd',txInd,'NumRadials',txsnumradials')
    txlat = txslats(txInd);
    txlon = txslons(txInd);
    txdatarange = datarange(txInd);
    numradials = txsnumradials(txInd);
    angs = linspace(0, 360, numradials+1);
    angs(end) = []; % Remove 360 since it is a duplicate of 0
    
    % Compute max distance along radials. For multiple txsites, use
    % distance from txsite to furthest data region of any other txsite.
    maxraddist = txdatarange;
    if isMultipleTx
        for otherTxInd = 1:numTx
            otherTx = txs(otherTxInd);
            if ~isequal(tx, otherTx)
                txToTxDist = rfprop.internal.MapUtils.greatCircleDistance(...
                    txlat, txlon, txslats(otherTxInd), txslons(otherTxInd));
                candmaxraddist = txToTxDist + datarange(otherTxInd);
                if candmaxraddist > maxraddist
                    maxraddist = candmaxraddist;
                end
            end
        end
    end
    
    % Compute end location for each radial by creating a vector of
    % candidate end locations (spaced by gridres) and choosing the
    % furthest one that is within range of any site.
    endlats = zeros(numradials,1);
    endlons = zeros(numradials,1);
    raddists = res:res:maxraddist; % Candidate end locations
    numangs = numradials;
    for angInd = 1:numangs
        % Trim candidate end locations to those within range
        [raddistlats, raddistlons] = rfprop.internal.MapUtils.greatCircleForward(...
            txlat, txlon, raddists, angs(angInd));
        [raddistlats, raddistlons] = rfprop.internal.MapUtils.georange(...
            txs, raddistlats, raddistlons, datarange, terrainSource);
        
        if isempty(raddistlats)
            % Skip radial since there is no value on it within range
            numradials = numradials - 1;
        else
            % Select furthest location as radial end location
            endlats(angInd) = raddistlats(end);
            endlons(angInd) = raddistlons(end);
        end
    end
    
    % Compute all terrain sample locations along radials
    endsites = rxsite(...
        'Name', 'internal.radialsites', ...
        'Latitude', endlats(:), ...
        'Longitude', endlons(:));
    [radslats, radslons, actterrainres] = sampleGreatCircle(tx, endsites, terrainres, ...
        'SourceAntennaSiteCoordinates', txCoords, 'Map', map);
    
    % Build lists of all locations and rx site locations along radials
    numlocs = zeros(numradials,1);
    allradslats = [];
    allradslons = [];
    rxlats = [];
    rxlons = [];    
    for k = 1:numradials
        radlats = radslats{k};
        radlons = radslons{k};
        
        % Remove tx from profile as optimization
        radlats = radlats(2:end);
        radlons = radlons(2:end);
        numlocs(k) = numel(radlats);
        
        % Grow list of all locations on radials
        allradslats = [allradslats; radlats]; %#ok<AGROW>
        allradslons = [allradslons; radlons]; %#ok<AGROW>
        
        % Grow list of rx site locations for computing signal strength
        rxlats = [rxlats; radlats(1:skipFactor:end)]; %#ok<AGROW>
        rxlons = [rxlons; radlons(1:skipFactor:end)]; %#ok<AGROW>
    end
    
    % Query elevation of all locations
    Z = rfprop.internal.AntennaSiteCoordinates.querySurfaceHeightAboveGeoid(...
        allradslats, allradslons, map);
    
    % Get terrain for terrain profiles and data locations
    numProfiles = numel(actterrainres);
    Zprofiles = cell(1,numProfiles);
    Zdata = [];
    locInd = 1;
    for k = 1:numProfiles
        % Get terrain of all locations on radial
        locIndEnd = locInd + numlocs(k) - 1;
        Zprofile = Z(locInd:locIndEnd);
        
        % Grow variables for terrain profiles and data locations
        Zprofiles{k} = Zprofile;
        Zdata = [Zdata; Zprofile(1:skipFactor:end)]; %#ok<AGROW>
        
        locInd = locIndEnd + 1;
    end
    
    % Package terrain profiles and receiver locations
    terrainProfiles = {Zprofiles, actterrainres, skipFactor};
    
    % Define rxsite locations
    rxs = rxsite(...
        'Name', 'internal.radialsite', ... % Specify to avoid default site naming
        'Latitude', rxlats, ...
        'Longitude', rxlons, ...
        'AntennaHeight', rxAntennaHeight);
    
    % Compute signal strength at each rxsite
    rxsCoords =  rfprop.internal.AntennaSiteCoordinates([rxlats rxlons Zdata], rxAntennaHeight, map);
    
    % Initiate second phase of full progress dialog. Do this only for first
    % tx if multiple transmitter sites.
    if isa(map,'siteviewer') && ~isempty(computingDataMsg) && fullProgressDlg && txInd == 1
       launchWaitbarDialog(map, computingDataMsg)
    end
    
    ss = sigstrength(rxs, tx, pm, ...
        'Type', type, ...
        'ReceiverGain', rxGain, ...
        'TransmitterAntennaSiteCoordinates', txCoords, ...
        'ReceiverAntennaSiteCoordinates', rxsCoords, ...
        'TerrainProfiles', terrainProfiles, ...
        'Map', mapInput);
    
    % If multiple txsites, then interpolate data at grid locations
    if isMultipleTx
        F = scatteredInterpolant(rxlons,rxlats,ss','natural','none');
        data(txInd,isInDataRange) = F(datalons(isInDataRange),datalats(isInDataRange));
    else
        data = ss;
    end
end

% If multiple txsites, then data corresponds to grid locations
if isMultipleTx
    rxlats = datalats;
    rxlons = datalons;
end
end

function launchBusyDialog(map, msg)

map.showProgressDialog('Message', msg, ...
    'Indeterminate', true, ...
    'Cancelable', false);
end

function launchWaitbarDialog(map, msg)

map.showProgressDialog('Message', msg, ...
    'Value', 0, ...
    'ShowPercentage', true, ...
    'Indeterminate', false, ...
    'Cancelable', true);
end