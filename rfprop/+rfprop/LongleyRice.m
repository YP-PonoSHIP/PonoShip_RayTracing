classdef LongleyRice < rfprop.TerrainPropagationModel
%

% Copyright 2018-2024 The MathWorks, Inc.

% Reference: G.A. Hufford, A.G. Longley, W.A. Kissick, A Guide to the Use
% of the ITS Irregular Terrain Model in the Area Prediction Mode, 1982.
    
    properties
        AntennaPolarization(1,:) char {mustBeAntennaPolarization(AntennaPolarization)} = 'horizontal'
        GroundConductivity(1,1) double {mustBeFinite, mustBeNonsparse, mustBeReal, mustBeNonnegative} = 0.005 % S/m
        GroundPermittivity(1,1) double {mustBeFinite, mustBeNonsparse, mustBeReal, mustBeNonnegative} = 15
        AtmosphericRefractivity(1,1) double {mustBeNonsparse, mustBeInRange(AtmosphericRefractivity,[250 400])} = 301 % N-units
        ClimateZone(1,:) char {mustBeClimateZone(ClimateZone)} = 'continental-temperate'
        TimeVariabilityTolerance(1,1) double {mustBeNonsparse, mustBeInRange(TimeVariabilityTolerance,[0.001 0.999])} = 0.5
        SituationVariabilityTolerance(1,1) double {mustBeNonsparse, mustBeInRange(SituationVariabilityTolerance,[0.001 0.999])} = 0.5
    end
       
    properties (Constant, Hidden)
        ClimateZoneChoices = {'equatorial', 'continental-subtropical', ...
           'maritime-subtropical', 'desert', 'continental-temperate', ...
           'maritime-over-land', 'maritime-over-sea'}
        AntennaPolarizationChoices = {'horizontal', 'vertical'}
    end
    
    methods
        function [pl, info] = pathloss(pm, rxs, txs, varargin)            
            [pl, info] = pathloss@rfprop.PropagationModel(pm, rxs, txs, varargin{:});
        end
    end
        
    methods (Access = protected)
        function [pl, txAngle, rxAngle] = pathlossOverTerrain(pm, rx, tx, res, Z, ~, txAngle, rxAngle)
            
            % The elevation profile for the Longley-Rice model must be
            % specified as an array containing one less than the number 
            % of elevation points, the distance between each point,
            % and the series of elevation points.
            elevation = [numel(Z)-1, res, Z];
            
            [pl, ~] = builtin('_longleyricepl', elevation, ...
                pm.ClimateZone , pm.AntennaPolarization, pm.GroundPermittivity, ...
                pm.GroundConductivity, pm.AtmosphericRefractivity, ...
                tx.AntennaHeight, rx.AntennaHeight, ...
                tx.TransmitterFrequency/1e6, pm.SituationVariabilityTolerance, pm.TimeVariabilityTolerance);
        end
        
        function lim = frequencyLimits(~)
            lim = [20e6 20e9];
        end
        
        function lim = antennaHeightLimits(~)
            lim = [0.5 3000];
        end
        
    end
    
end

% AntennaPolarization parameter validation 
function mustBeAntennaPolarization(param)
    mustBeMember(param, rfprop.LongleyRice.AntennaPolarizationChoices);
end

% Climate zone parameter validation
function mustBeClimateZone(param)
    mustBeMember(param, rfprop.LongleyRice.ClimateZoneChoices);
end

% Range checking for parameter validation
function mustBeInRange(param, range)
    mustBeGreaterThanOrEqual(param, range(1));
    mustBeLessThanOrEqual(param, range(2));
end