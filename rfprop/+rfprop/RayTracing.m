classdef RayTracing < rfprop.PropagationModel & matlab.mixin.CustomDisplay
%

% Copyright 2019-2024 The MathWorks, Inc.
    
    properties
        Method(1,1) string = "sbr"
        AngularSeparation = "medium"
        MaxNumReflections(1,1) double {mustBeNonsparse, mustBeInteger, ...
            mustBeInRange(MaxNumReflections, 0, 100)} = 2
        MaxNumDiffractions(1,1) double {mustBeNonsparse, mustBeInteger, ...
            mustBeInRange(MaxNumDiffractions, 0, 2)} = 0
        MaxAbsolutePathLoss(1,1) double {mustBeNonsparse, mustBePositive} = Inf % dB
        MaxRelativePathLoss(1,1) double {mustBeNonsparse, mustBeNonnegative} = 40 % dB
        CoordinateSystem(1,1) string = "geographic"

        BuildingsMaterial(1,1) string = "auto"
        BuildingsMaterialPermittivity(1,1) double {mustBeFinite, mustBeNonsparse, mustBeReal} = 5.31
        BuildingsMaterialConductivity(1,1) double {mustBeNonsparse, mustBeNonnegative} = 0.0548 % S/m

        TerrainMaterial(1,1) string = "concrete"
        TerrainMaterialPermittivity(1,1) double {mustBeFinite, mustBeNonsparse, mustBeReal} = 5.31
        TerrainMaterialConductivity(1,1) double {mustBeNonsparse, mustBeNonnegative} = 0.0548 % S/m
        
        SurfaceMaterial(1,1) string = "auto"
        SurfaceMaterialPermittivity(1,1) double {mustBeFinite, mustBeNonsparse, mustBeReal} = 5.31
        SurfaceMaterialConductivity(1,1) double {mustBeNonsparse, mustBeNonnegative} = 0.0548 % S/m
    end

    properties (Hidden)
        UseGPU(1,1) string = "off"
    end

    properties (Constant, Hidden)
        MethodChoices = ["image", "sbr"]
        AngularSeparationChoices = ["high", "medium", "low"]
        CoordinateSystemChoices = ["geographic", "cartesian"]
        BuildingsMaterialChoices = ["perfect-reflector"
            "custom"
            "concrete"
            "brick"
            "wood"
            "glass"
            "metal"
            "plywood"
            "marble"
            "auto"]
        TerrainMaterialChoices = ["perfect-reflector"
            "custom"
            "concrete"
            "brick"
            "marble"
            "water"
            "vegetation"
            "loam"]
        SurfaceMaterialChoices = ["perfect-reflector"
            "custom"
            "plasterboard"
            "ceiling-board"
            "chipboard"
            "floorboard"
            "concrete"
            "brick"
            "wood"
            "glass"
            "metal"
            "plywood"
            "marble"
            "water"
            "vegetation"
            "loam"
            "auto"]
        UseGPUChoices = ["off", "on", "auto"]
    end
    
    methods        
        function [pl,info] = pathloss(pm,varargin)            
            [pl,info] = pathloss@rfprop.PropagationModel(pm, varargin{:});
        end
        
        function pm3 = add(pm1, pm2)            
            pm3 = add@rfprop.PropagationModel(pm1, pm2);
        end

        % Set methods
        function obj = set.AngularSeparation(obj, val)
            if isnumeric(val)
                validateattributes(val,{'numeric'},{'scalar','>=',0.05,'<=',10},'','AngularSeparation')
                obj.AngularSeparation = double(val);
            else
                obj.AngularSeparation = validatestring(val, obj.AngularSeparationChoices);
            end
        end
        function obj = set.BuildingsMaterial(obj, val)
            obj.BuildingsMaterial = validatestring(val, obj.BuildingsMaterialChoices);
        end
        function obj = set.CoordinateSystem(obj, val)
            obj.CoordinateSystem = validatestring(val, obj.CoordinateSystemChoices);
        end
        function obj = set.MaxNumDiffractions(obj, val)
            try
                validateUseGPU(val, obj.UseGPU);
            catch e
                throwAsCaller(e)
            end

            obj.MaxNumDiffractions = val;
        end
        function obj = set.Method(obj, val)
            obj.Method = validatestring(val, obj.MethodChoices);
        end
        function obj = set.TerrainMaterial(obj, val)
            obj.TerrainMaterial = validatestring(val, obj.TerrainMaterialChoices);
        end
        function obj = set.SurfaceMaterial(obj, val)
            obj.SurfaceMaterial = validatestring(val, obj.SurfaceMaterialChoices);
        end
        function obj = set.UseGPU(obj, val)
            val = validatestring(val, obj.UseGPUChoices);
            try
                validateUseGPU(obj.MaxNumDiffractions, val);
            catch e
                throwAsCaller(e)
            end
            matlab.internal.parallel.resolveUseGPU(val);

            obj.UseGPU = val;
        end
    end
    
    methods(Access = protected)
        function group = getPropertyGroups(pm)
            
            if ~isscalar(pm)
                group = getPropertyGroups@matlab.mixin.CustomDisplay(pm);
            else
                propList = {'CoordinateSystem','Method'};
                if strcmpi(pm.Method, 'sbr')
                    propList = [propList, {'AngularSeparation', ...
                        'MaxNumReflections', 'MaxNumDiffractions', ...
                        'MaxAbsolutePathLoss','MaxRelativePathLoss'}];
                else % image method
                    propList = [propList, {'MaxNumReflections', ...
                        'MaxAbsolutePathLoss','MaxRelativePathLoss'}];
                end
                if strcmpi(pm.CoordinateSystem,'geographic')
                    if strcmpi(pm.BuildingsMaterial,'custom')
                        propList = [propList {'BuildingsMaterial','BuildingsMaterialPermittivity','BuildingsMaterialConductivity'}];
                    else
                        propList = [propList {'BuildingsMaterial'}];
                    end
                    if strcmpi(pm.TerrainMaterial,'custom')
                        propList = [propList {'TerrainMaterial','TerrainMaterialPermittivity','TerrainMaterialConductivity'}];
                    else
                        propList = [propList {'TerrainMaterial'}];
                    end
                else
                    if strcmpi(pm.SurfaceMaterial,'custom')
                        propList = [propList {'SurfaceMaterial','SurfaceMaterialPermittivity','SurfaceMaterialConductivity'}];
                    else
                        propList = [propList {'SurfaceMaterial'}];
                    end
                end
                group = matlab.mixin.util.PropertyGroup(propList);
            end
        end
        
        function [pl,info] = pathlossMultipath(pm, rxs, txs, map)
            
            rays = raytrace(txs, rxs, ...
                'Type', 'pathloss', ...
                'PropagationModel', pm, ...
                'Map', map);
            
            % Outputs are M-by-N cell arrays
            pl = cell(numel(txs), numel(rxs));
            info = cell(numel(txs), numel(rxs));
            
            % Analyze ray objects to compute corresponding outputs
            for txInd = 1:numel(txs)
                for rxInd = 1:numel(rxs)
                    % Get rays, where each ray object corresponds to a propagation
                    % path between the tx and rx (whether LOS or reflected)
                    txrxRays = rays{txInd,rxInd};
                    
                    % Analyze each ray
                    numRays = numel(txrxRays);
                    if numRays == 0
                        continue
                    end
                    
                    % Initialize path loss and info structs entry
                    pls = ones(1,numRays);
                    pathStructs = struct(...
                        "PropagationDistance", cell(1,numRays), ...
                        "AngleOfDeparture", cell(1,numRays), ...
                        "AngleOfArrival", cell(1,numRays), ...
                        "NumReflections", cell(1,numRays), ...
                        "NumDiffractions", cell(1,numRays));
                    
                    % Populate path loss and info structs
                    for rayInd = 1:numRays
                        ray = txrxRays(rayInd);
                        
                        % Populate path loss
                        pls(1,rayInd) = ray.PathLoss;
                        
                        % Populate info struct
                        pathStructs(rayInd) = struct(...
                            "PropagationDistance", ray.PropagationDistance, ...
                            "AngleOfDeparture", ray.AngleOfDeparture, ...
                            "AngleOfArrival", ray.AngleOfArrival, ...
                            "NumReflections", ray.NumReflections, ...
                            "NumDiffractions", ray.NumDiffractions);
                    end
                    
                    pl{txInd,rxInd} = pls;
                    info{txInd,rxInd} = pathStructs;
                end
            end
        end
        
        function lim = frequencyLimits(~)
            lim = [100e6 100e9];
        end
    end
    
    methods(Hidden)
        function rt = isMultipathModel(~)
            rt = true;
        end
    end
end

function validateUseGPU(maxNuMDiff, useGPU)
if maxNuMDiff > 0 && ~strcmp(useGPU, "off")
    error(message('shared_channel:rfprop:UseGPUOffForDiffraction'))
end
end