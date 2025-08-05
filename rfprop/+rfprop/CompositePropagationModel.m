classdef (Sealed) CompositePropagationModel < rfprop.PropagationModel
%CompositePropagationModel   Composite propagation model
%   PMC = add(PM1,PM2) adds propagation model objects PM1 and PM2 and
%   returns a composite propagation model object PMC which contains PM1 and
%   PM2.
%
%   The path loss computed by PMC is the sum of path losses computed by PM1
%   and PM2. If either PM1 or PM2 is a ray tracing propagation model, then
%   PMC is also a ray tracing propagation model where path losses from
%   rain, gas, or fog models in the composite are added to the path loss
%   computed for each propagation path.
%
%   Notes
%   -----
%   - The syntax 'PM1 + PM2' can be used in place of add.
%   - A composite propagation model cannot contain more than one
%     propagation model object of the same class.
%   - A composite propagation model cannot contain more than one
%     propagation model object which includes effects of free-space loss.
%
%   CompositePropagationModel properties:
%
%      PropagationModels - Propagation model objects
%
%   CompositePropagationModel methods:
%
%      add      - Add propagation models
%      pathloss - Path loss of radio wave propagation
%
%   Example: Calculate signal strength over terrain in rain
%
%   % Specify transmitter and receiver sites
%   tx = txsite("Name","Fenway Park", ...
%       "Latitude",42.3467, ...
%       "Longitude",-71.0972, ...
%       "TransmitterFrequency",6e9);
%   rx = rxsite("Name","Bunker Hill Monument", ...
%       "Latitude",42.3763, ...
%       "Longitude",-71.0611);
%
%   % Calculate signal strength using default Longley-Rice model
%   ss1 = sigstrength(rx,tx)
%
%   % Create composite propagation model with Longley-Rice and specific
%   % atmospheric propagation models
%   pm = propagationModel("longley-rice") + ...
%      propagationModel("gas") + propagationModel("rain");
%
%   % Calculate signal strength using composite propagation model
%   ss2 = sigstrength(rx,tx,pm)

%   Copyright 2018-2021 The MathWorks, Inc.
    
    properties
        %PropagationModels   Propagation model objects
        %   Propagation model objects in the composite propagation model,
        %   specified as as an object array.
        PropagationModels(1,:) rfprop.PropagationModel {validatePropagationModels} = rfprop.PropagationModel.empty
    end
    
    % Only allow explicit construction from propagation models
    methods(Access = ?rfprop.PropagationModel)
        function pm = CompositePropagationModel(varargin)
            if nargin % Add each input model to this
                for modInd = 1:numel(varargin)
                    pm = pm + varargin{modInd};
                end
            end
        end
    end
    
    methods
        function pm = set.PropagationModels(pm, models)
            validateattributes(models, {'rfprop.PropagationModel'}, {'nonempty'},'','PropagationModels');
            validatePropagationModels(models)
            pm.PropagationModels = models;
        end
    end
    
    methods(Hidden)
        function pm = plus(pm, others)
            
            validateattributes(others, {'rfprop.PropagationModel'}, {'scalar'});
            
            % Get array of other models if adding with another composite
            if isa(others, 'rfprop.CompositePropagationModel')
                others = others.PropagationModels;
            end
            
            % Grow list of models in composite
            models = [pm.PropagationModels, others];
            validatePropagationModels(models)
            pm.PropagationModels = models;
        end
        
        function rt = requiresTerrain(pm)
            % Return true if any model requires terrain
            rt = false;
            for model = pm.PropagationModels
                if requiresTerrain(model)
                    rt = true;
                    return
                end
            end
        end
        
        function ismultipath = isMultipathModel(pm)
            % Return true if any model is multipath
            
            ismultipath = false;
            for model = pm.PropagationModels
                if isMultipathModel(model)
                    ismultipath = true;
                    return
                end
            end
        end
        
        function tr = terrainProfileResolution(pm, map)
            % Return lowest resolution
            models = pm.PropagationModels;
            numModels = numel(models);
            trs = zeros(numModels, 1);
            for modInd = 1:numModels
                trs(modInd) = terrainProfileResolution(models(modInd), map);
            end
            tr = min(trs);
        end
        
        function pl = specificPathLoss(pm, rx, tx, ds, angleOfDeparture)
            % Return cumulative specific path loss
            
            pl = 0;
            aod_el = angleOfDeparture(2);
            for model = pm.PropagationModels
                if ~isMultipathModel(model)
                    pl = pl + model.pathlossOverDistance(rx, tx, ds, aod_el);
                end
            end
        end
    end
    
    methods(Access = protected)
        function pm = pathlossSetup(pm, varargin)
            models = pm.PropagationModels;
            for modInd = 1:numel(models)
                models(modInd) = pathlossSetup(models(modInd),varargin{:});
            end
            pm.PropagationModels = models;
        end
        
        function pl = pathlossOverDistance(pm, rxs, tx, d, el)
            % Return sum of path losses for propagation models
            pl = 0;
            for model = pm.PropagationModels
                pl = pl + pathlossOverDistance(model, rxs, tx, d, el);
            end
        end
        
        function [pl,txAngle, rxAngle] = pathlossOverTerrain(pm, rx, tx, res, Z, d, txAngle, rxAngle)
            % Return sum of path losses over terrain or distance for propagation models
            
            pl = 0;
            for model = pm.PropagationModels
                if model.requiresTerrain
                    [modpl, txAngle, rxAngle] = pathlossOverTerrain(model, rx, tx, res, Z, d, txAngle, rxAngle);
                else
                    modpl = pathlossOverDistance(model, rx, tx, d, txAngle);
                end
                pl = pl + modpl;
            end
        end
        
        function pl = pathlossOverVerticalDistance(pm, rx, tx, d, el)
            % Return sum of path losses over vertical distance for propagation models
            
            pl = 0;
            for model = pm.PropagationModels
                pl = pl + pathlossOverVerticalDistance(model, rx, tx, d, el);
            end
        end
        
        function [pls,infos] = pathlossMultipath(pm, rxs, txs, map)
            
            % Find multipath model
            multipathModel = [];
            models = pm.PropagationModels;
            for model = models
                if isMultipathModel(model)
                    multipathModel = model;
                    break
                end
            end

            % Get the raytrace path loss
            [pls,infos] = multipathModel.pathlossMultipath(rxs, txs, map);
            
            % Augment path loss using rest of models
            for txInd = 1:numel(txs)
                tx = txs(txInd);
                for rxInd = 1:numel(rxs)
                    rx = rxs(rxInd);
                    pathsPl = pls{txInd,rxInd};
                    pathsInfo = infos{txInd,rxInd};
                    for pathInd = 1:numel(pathsPl)
                        pathInfo = pathsInfo(pathInd);
                        pathsPl(pathInd) = pathsPl(pathInd) + pm.specificPathLoss(rx, tx, ...
                            pathInfo.PropagationDistance, pathInfo.AngleOfDeparture);
                    end

                    % Filter rays using path loss thresholds
                    if any(isfinite([multipathModel.MaxAbsolutePathLoss multipathModel.MaxRelativePathLoss]))
                        keepRays = pathsPl <= min(multipathModel.MaxAbsolutePathLoss, ...
                            min(pathsPl)+multipathModel.MaxRelativePathLoss);
                        pls{txInd,rxInd} = pathsPl(keepRays);
                        infos{txInd,rxInd} = pathsInfo(keepRays);
                    else
                        pls{txInd,rxInd} = pathsPl;
                    end
                end
            end
        end
        
        function fs = includesFreeSpaceContribution(pm)
            % Return true if any model includes free space contribution
            fs = false;
            for model = pm.PropagationModels
                if includesFreeSpaceContribution(model)
                    fs = true;
                    return
                end
            end
        end
        
        function fqlim = frequencyLimits(pm)
            % Return common frequency limits of all propagation models
            models = pm.PropagationModels;
            numModels = numel(models);
            fqlims = ones(numModels, 2);
            for modInd = 1:numModels
                fqlims(modInd,:) = frequencyLimits(models(modInd));
            end
            fqlim = [max(fqlims(:,1)) min(fqlims(:,2))];
        end
        
        function htlim = antennaHeightLimits(pm)
            % Return merged limits of all propagation models
            models = pm.PropagationModels;
            numModels = numel(models);
            htlims = ones(numModels, 2);
            for modInd = 1:numModels
                htlims(modInd,:) = antennaHeightLimits(models(modInd));
            end
            htlim = [max(htlims(:,1)) min(htlims(:,2))];
        end
    end
end

function validatePropagationModels(models)

try
    % Get info for each model: class name and whether it includes free space
    numModels = numel(models);
    classNames = cell(1,numModels);
    includesFreeSpace = true(1,numModels);
    for modelInd = 1:numModels
        model = models(modelInd);
        classNames{modelInd} = class(model);
        includesFreeSpace(modelInd) = model.includesFreeSpaceContribution;
    end
    
    % Perform validation
    uniqueClassNames = unique(classNames);
    if any(strcmp('rfprop.CompositePropagationModel',uniqueClassNames))
        error(message('shared_channel:rfprop:CompositeInvalidCompositePropagationModel'));
    elseif numel(uniqueClassNames) < numModels
        error(message('shared_channel:rfprop:CompositeInvalidPlusSameClass'));
    elseif numel(find(includesFreeSpace)) > 1
        error(message('shared_channel:rfprop:CompositeInvalidPlusFreeSpaceContributions'));
    end
catch e
    throwAsCaller(e)
end
end
