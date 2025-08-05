classdef TIREM < rfprop.TerrainPropagationModel
%

%   Copyright 2018-2024 The MathWorks, Inc.

    properties
        AntennaPolarization(1,:) char = 'horizontal'
        GroundConductivity(1,1) double {tirem.internal.Validators.validateGroundConductivity} = 0.005 % S/m
        GroundPermittivity(1,1) double {tirem.internal.Validators.validateGroundPermittivity} = 15
        AtmosphericRefractivity(1,1) double {tirem.internal.Validators.validateAtmosphericRefractivity} = 301 % N-units
        Humidity(1,1) double {tirem.internal.Validators.validateHumidity} = 9 % g/m^3
    end
    
    properties(Access = private)
        PropagationAnalyzer
    end
        
    methods
        function [pl, info] = pathloss(pm, rxs, txs, varargin)            
            [pl, info] = pathloss@rfprop.PropagationModel(pm, rxs, txs, varargin{:});
        end
        
        function pm = set.AntennaPolarization(pm, pol)
            % Validate AntennaPolarization and assign full match value
            pm.AntennaPolarization = tirem.internal.Validators.validateAntennaPolarization(pol);
        end
    end
        
    methods(Access = protected)        
        function pm = pathlossSetup(pm, varargin)
            pm.PropagationAnalyzer = tirem.internal.PropagationAnalyzer;
        end

        function [pl, txAngle, rxAngle] = pathlossOverTerrain(pm, rx, tx, res, Z, ~, txAngle, rxAngle)
            
            % Get distances for terrain profile, noting that profile is
            % evenly spaced by res.
            numProfilePoints = numel(Z);
            ds = 0:res:(res*(numProfilePoints-1));

            % Run analysis
            analyzer = pm.PropagationAnalyzer;
            analyzer.analyze(tx.AntennaHeight, rx.AntennaHeight, tx.TransmitterFrequency, Z, ds, false, ...
                pm.AtmosphericRefractivity, pm.GroundConductivity, pm.GroundPermittivity, pm.Humidity, pm.AntennaPolarization, 0, 0, 0);
            pl = analyzer.TotalLoss;
            
            % Compute output angle. Only update from default if beyond
            % line-of-sight.
            output = analyzer.Outputs;
            propMode = output.PropagationMode;
            if propMode == "DIF"
                txAngle = rad2deg(output.TransmitterHorizonAngle);
                rxAngle = rad2deg(output.ReceiverHorizonAngle);
            elseif propMode == "TRO"
                txAngle = rad2deg(output.TransmitterTroposcatterAngle);
                rxAngle = rad2deg(output.ReceiverTroposcatterAngle);
            end
        end
                
        function lim = frequencyLimits(~)
            lim = [1e6 1000e9]; % 1 MHZ to 1000 GHz
        end
        
        function lim = antennaHeightLimits(~)
            lim = [0 30000];
        end
    end
end