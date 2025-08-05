function ss = sigstrength(rxs, txs, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.   

% Validate sites
validateattributes(rxs,{'rxsite'},{'nonempty'},'sigstrength','',1);
validateattributes(txs,{'txsite'},{'nonempty'},'sigstrength','',2);


% Allocate output matrix
numTx = numel(txs);
numRx = numel(rxs);
ss = zeros(numTx, numRx);

% Process optional inputs
p = inputParser;
if nargin > 2 && mod(numel(varargin),2)
    % Validator function is necessary for inputParser to allow string
    % option instead of treating it like parameter name
    p.addOptional('PropagationModel', [], @(x)ischar(x)||isstring(x)||isa(x,'rfprop.PropagationModel'));
else
    p.addParameter('PropagationModel', []);
end
p.addParameter('Type', 'power');
p.addParameter('ReceiverGain', []);
p.addParameter('Map', []);
p.addParameter('TransmitterAntennaSiteCoordinates', []);
p.addParameter('ReceiverAntennaSiteCoordinates', []);
p.addParameter('TerrainProfiles', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.TransmitterAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(rxs, txs);
else
    usingCartesian = strcmp(p.Results.TransmitterAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

% Validate and get parameters
if usingCartesian
    [map, ~, ~, mapStruct] = rfprop.internal.Validators.validateCartesianMap(p);
    pm = rfprop.internal.Validators.validateCartesianPropagationModel(p, map, 'sigstrength');
    if isprop(pm, 'CoordinateSystem')
        pm.CoordinateSystem = 'cartesian';
    end
else
    [map, mapStruct] = rfprop.internal.Validators.validateMapTerrainSource(p, 'sigstrength');
    rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, 'geographic');
    pm = rfprop.internal.Validators.validateGeographicPropagationModel(p, map, 'sigstrength');
end
rfprop.internal.Validators.validateMaxNumReflections(pm, 'sigstrength');
isMultipathModel = pm.isMultipathModel;
typeIsEfield = validateType(p);
if typeIsEfield
    % Assuming far-field conditions and reciprocity:  
    %   radiated power density Wt (or St) = Pt*Gt*Lt*Fp^2/(4*pi*r^2);
    %       Pr = Wt*Gr*Lr*Fpol^2*lambda^2/(4*pi) [Friss Transmission Equation];
    %       and Wt = |E x conj(H)|avg = |Erms|^2 / eta. Fp is the
    %       propagation factor accounting for scenario complexity (=1 for
    %       homogenous medium e.g. antennas in free-space) and eta is the
    %       characteristic impedance of the propagation medium.
    % pathloss PL = (4*pi*r/(lambda*Fp))^2
    %   => Fp^2/(4*pi*r^2) = 4*pi/(PL*lambda^2)
    % => Erms = sqrt(eta*Wt) = sqrt(eta*Pt*Gt*4*pi/(PL*lambda^2))
    %   Note:  <ray>.PathLoss includes polarization mismatch (Fpol^2) so
    %   the electric field strength calculated herein is actually that
    %   which would be measured by the specified receive antenna (including
    %   its orientation).
    Z0 = rfprop.Constants.Z0; % propagation medium assumed to be free-space
    EfieldConst = 10*log10(Z0*4*pi) + 120; % dBuV/m offset constant
end

rxGain = validateReceiverGain(p, numTx);

% Get site antenna coordinates
txsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TransmitterAntennaSiteCoordinates, txs, map, 'sigstrength');
rxsCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.ReceiverAntennaSiteCoordinates, rxs, map, 'sigstrength');

% Calculate path loss depending on propagation model type
if isMultipathModel
    rays = raytrace(txs, rxs, pm, Map=mapStruct);
    Lpls = rays;
    rxsAoas = rays;
else
    % Get path loss and propagation path info
    [Lpls_db, info] = pm.pathloss(rxs, txs, 'Map', map, ...
        'TransmitterAntennaSiteCoordinates', txsCoords, ...
        'ReceiverAntennaSiteCoordinates', rxsCoords, ...
        'ComputeAngleOfArrival', isempty(rxGain), ...
        'TerrainProfiles', p.Results.TerrainProfiles);
end

% Optimize pattern computation with re-usable directivity objects when
% transmitters share frequency and all tx antennas or all rx antennas are
% the same Phased Array object.
isSharedTxFrequency = all(arrayfun(@(x)x.TransmitterFrequency==txs(1).TransmitterFrequency, txs(:)));
sharedTxAntenna = txs(1).Antenna;
useSharedTxDirectivity = (numTx > 1) && ~isMultipathModel && isSharedTxFrequency && ...
    rfprop.AntennaSite.isPhasedAntenna(sharedTxAntenna) && ...
    all(arrayfun(@(x)x.Antenna==sharedTxAntenna, txs(:)));
if useSharedTxDirectivity
    txDirectivity = phased.internal.Directivity('Sensor',sharedTxAntenna);
end

sharedRxAntenna = rxs(1).Antenna;
useSharedRxDirectivity = (numRx > 1) && ~isMultipathModel && isSharedTxFrequency && ...
    rfprop.AntennaSite.isPhasedAntenna(sharedRxAntenna) && ...
    all(arrayfun(@(x)x.Antenna==sharedRxAntenna, rxs(:)));
if useSharedRxDirectivity
    rxDirectivity = phased.internal.Directivity('Sensor',sharedRxAntenna);
end

% Compute signal strength from each transmitter to all receivers
for txInd = 1:numTx
    tx = txs(txInd);

    % Compute transmitter constants
    fq = tx.TransmitterFrequency;
    Ptx = tx.TransmitterPower;
    if typeIsEfield
        lambda = rfprop.Constants.LightSpeed/fq;
        txEfieldConst = EfieldConst + 10*log10(Ptx) - 20*log10(lambda); % dBuV/m offset constant
    else
        Ptx_db = 10 * log10(1000*Ptx); % Convert W to dBm (db with reference to mW)
    end
    Ltxsys_db = tx.SystemLoss;

    % Get directivity or gain pattern
    useTxDirectivity = useSharedTxDirectivity || rfprop.AntennaSite.isPhasedAntenna(tx.Antenna);
    useTxGainPattern = rfprop.AntennaSite.isElectromagneticAntenna(tx.Antenna);
    if useTxDirectivity && ~useSharedTxDirectivity
        txDirectivity = phased.internal.Directivity('Sensor',tx.Antenna);
    elseif useTxGainPattern
        [Gtx,Gtxaz,Gtxel] = gainPattern(tx,fq);
    end

    % Get AoD for all propagation paths from tx. For multipath models e.g.
    % ray-trace, also need AoA for all tx-rx pairs and corresponding path
    % loss; path loss will be kept in dB if only 1 ray found per site pair,
    % otherwise path loss will be a unitless field strength.
    if isMultipathModel
        txAods = [];
        numRays = zeros(1,numRx);
        for rxInd = 1:numRx
            rxRays = rays{txInd,rxInd};
            if ~isempty(rxRays)
                rxNumRays = numel(rxRays);
                numRays(rxInd) = numel(rxRays);
                if rxNumRays == 1
                    txAods = [txAods rxRays.AngleOfDeparture]; %#ok<AGROW>
                    rxsAoas{txInd,rxInd} = rxRays.AngleOfArrival; % keep in dB
                    Lpls{txInd,rxInd} = rxRays.PathLoss;
                else
                    rxAoas_temp = zeros(2,rxNumRays);
                    Lpls_temp = zeros(rxNumRays,1);
                    for i = 1:rxNumRays
                        txAods = [txAods rxRays(i).AngleOfDeparture]; %#ok<AGROW>
                        rxAoas_temp(:,i) = rxRays(i).AngleOfArrival;
                        Lpls_temp(i) = exp(1i*rxRays(i).PhaseShift)/ ...
                            10^(rxRays(i).PathLoss/20);
                    end
                    rxsAoas{txInd,rxInd} = rxAoas_temp;
                    Lpls{txInd,rxInd} = Lpls_temp; % unitless field strength
                end
            end
            
        end
    else
        txAods = [info(txInd,:).AngleOfDeparture];
    end

    % Compute gains for all propagation paths from tx
    if ~isempty(txAods)
        txAz = txAods(1,:)';
        txEl = txAods(2,:)';
        if useTxDirectivity
            Gtxrxs_db = gain(tx,fq,txAz,txEl,txDirectivity);
        elseif useTxGainPattern
            Gtxrxs_db = gain(tx,fq,txAz,txEl,Gtx,Gtxaz,Gtxel);
        else
            Gtxrxs_db = gain(tx,fq,txAz,txEl);
        end
    end

    aodStartInd = 1;
    for rxInd = 1:numRx
        % Compute transmitter gain (including system losses) and path loss
        if isMultipathModel
            % Get path loss and check for empty, which corresponds to no
            % propagation paths found
            Lpl_ul = Lpls{txInd,rxInd};
            if isempty(Lpl_ul) % Infinite path loss since no propagation path
                ss(txInd,rxInd) = -inf;
                continue
            end

            % Get gains for this rx
            aodEndInd = aodStartInd + numRays(rxInd) - 1;
            Gtx_db = Gtxrxs_db(aodStartInd:aodEndInd) - Ltxsys_db;
            aodStartInd = aodEndInd + 1;
        else
            % Get path loss
            Lpl_db = Lpls_db(txInd, rxInd);

            Gtx_db = Gtxrxs_db(rxInd) - Ltxsys_db;
        end

        % Compute signal strength depending on if calculating E-Field
        % strength or received power. For E-Field strength, do not need to
        % know rx gain; for received power, do need to know rx gain.
        if typeIsEfield % dBuv/m
            % For multipath model e.g. ray-trace, use coherent phasor sum
            % of the individual path losses.
            if isMultipathModel
                % Perform phasor sum. If only 1 ray was found for this site
                % pair, the path loss will still be in dB, otherwise the
                % path loss will be unitless field strength.
                if numRays(rxInd) == 1 % Note: == 0 case is previously handled
                    ss(txInd,rxInd) = txEfieldConst + Gtx_db - Lpl_ul; % Lpl_ul is still in dB for only 1 ray
                else
                    E = sum(10.^(Gtx_db./20).*Lpl_ul);
                    ss(txInd,rxInd) = txEfieldConst + 20*log10(abs(E));
                end
            else
                ss(txInd,rxInd) = txEfieldConst + Gtx_db - Lpl_db; 
            end
        else % dBm
            % Compute receiver gain, including system losses
            if ~isempty(rxGain)
                Grx_db = repmat(rxGain(txInd), size(Gtx_db));
            else
                rx = rxs(rxInd);

                % Get directivity or gain pattern
                useRxDirectivity = useSharedRxDirectivity || rfprop.AntennaSite.isPhasedAntenna(rx.Antenna);
                useRxGainPattern = rfprop.AntennaSite.isElectromagneticAntenna(rx.Antenna);
                if useRxDirectivity && ~useSharedRxDirectivity
                    rxDirectivity = phased.internal.Directivity('Sensor',rx.Antenna);
                elseif useRxGainPattern
                    [Grx,Grxaz,Grxel] = gainPattern(rx,fq);
                end

                if isMultipathModel
                    rxAoas = rxsAoas{txInd,rxInd};
                else
                    infoStructs = info(txInd,rxInd);
                    rxAoas = [infoStructs.AngleOfArrival];
                end

                rxAz = rxAoas(1,:)';
                rxEl = rxAoas(2,:)';

                if useRxDirectivity
                    Grx_db = gain(rx,fq,rxAz,rxEl,rxDirectivity) - rx.SystemLoss;
                elseif useRxGainPattern
                    Grx_db = gain(rx,fq,rxAz,rxEl,Grx,Grxaz,Grxel) - rx.SystemLoss;
                else
                    Grx_db = gain(rx,fq,rxAz,rxEl) - rx.SystemLoss;
                end
            end

            % Compute signal strength in dBm, using link budget form of
            % Friis equation:
            %
            %  Received Power (dBm) = Transmitted Power (dBm) + Gains (dB) 
            %  - Losses (dB)
            %
            % Applied to tx/rx, this yields:
            %
            %  Prx_db = Ptx_db + Gtx_db + Grx_db - Lpl_db
            %
            % where:
            %  * Prx_db is received power in dBm at receiver input
            %  * Ptx_db is transmitter output power in dBm
            %  * Gtx_db is transmitter system gain in dBi (antenna gain - system loss)
            %  * Grx_db is receiver system gain in dBi (antenna gain - system loss)
            %  * Lpl_db is path loss (dB) as given by propagation model
            %
            % For multipath model e.g. ray-trace, use coherent phasor sum
            % of the individual path losses.
            if isMultipathModel
                % Perform phasor sum. If only 1 ray was found for this site
                % pair, the path loss will still be in dB, otherwise the
                % path loss will be unitless field strength.
                if numRays(rxInd) == 1 % 0 case is already taken care of
                    ss(txInd,rxInd) = Ptx_db + Gtx_db + Grx_db - Lpl_ul; % Lpl_ul is still in dB for only 1 ray
                else
                    E = sum(10.^(Gtx_db./20).*Lpl_ul.*10.^(Grx_db./20));
                    ss(txInd,rxInd) = Ptx_db + 20*log10(abs(E));
                end
            else
                ss(txInd,rxInd) = Ptx_db + Gtx_db + Grx_db - Lpl_db;
            end
        end
    end
end
end

function typeIsEfield = validateType(p)

try
    type = p.Results.Type;
    type = validatestring(type, {'efield','power'}, 'sigstrength','Type');
    typeIsEfield = strcmp(type,"efield");
catch e
    throwAsCaller(e);
end
end

function rxGain = validateReceiverGain(p, numTx)

try
    rxGain = p.Results.ReceiverGain;
    if ~isempty(rxGain)
        validateattributes(rxGain,{'numeric'}, ...
            {'real','finite','nonnan','nonsparse','nonempty'}, 'sigstrength', 'ReceiverGain');
        
        % Expand scalar gain to match length of tx
        if isscalar(rxGain)
            rxGain = repmat(rxGain,1,numTx);
        end
    end
catch e
    throwAsCaller(e);
end
end