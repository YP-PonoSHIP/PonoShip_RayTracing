classdef (Abstract) TerrainPropagationModel < rfprop.PropagationModel
    %TerrainPropagationModel   Abstract class for terrain propagation models
    
    %   Copyright 2018-2021 The MathWorks, Inc.
    
    properties(Hidden)
        %TerrainResolution   Resolution of terrain sample locations (m)
        %   Resolution of sample locations used to generate a terrain
        %   profile for computing path loss, specified as a numeric scalar
        %   in meters. The resolution defines the distance between samples
        %   on the great circle path between sites, using using a spherical
        %   Earth model. The default value is 'auto'.
        TerrainResolution = 'auto'
    end
        
    methods(Access = protected)
        function pathlossOverDistance(varargin)
            error(message('shared_channel:rfprop:PropagationModelNoPathlossOverDistance'));
        end
        
        function pl = pathlossOverVerticalDistance(~, ~, tx, d, ~)
            
            % Use free space loss when no terrain profile possible
            lambda = rfprop.Constants.LightSpeed/tx.TransmitterFrequency;
            pl = fspl(d, lambda);
        end
    end
    
    methods(Hidden)
        % Hide "range" method and throw error if called
        function varargout = range(varargin) %#ok<STOUT> 
            error(message('shared_channel:rfprop:PropagationModelNoRange'));
        end
        
        function rt = requiresTerrain(~)
            rt = true;
        end
        
        function resm = terrainProfileResolution(pm, map)
            
            % Return special resolution value if map contains buildings
            % data
            viewer = [];
            if isa(map,'siteviewer')
                viewer = map;
            end
            if ~isempty(viewer) && viewer.HasBuildings
                resm = rfprop.Constants.TerrainProfileResolutionWithBuildings;
                return
            end
            
            % Compute resolution (in meters) from terrain source
            terrainSource = rfprop.internal.Validators.validateTerrainSource(map);
            resm = pm.TerrainResolution;
            if ischar(resm) && strcmpi(resm,'auto')
                resm = terrainSource.IntrinsicResolutionArcLength;
            end
        end
    end
end