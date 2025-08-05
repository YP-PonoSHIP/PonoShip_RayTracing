classdef Fog < rfprop.PropagationModel
%

% Copyright 2017-2024 The MathWorks, Inc.  
        
    properties
        Temperature(1,1) double {mustBeFinite, mustBeReal, mustBeNonsparse} = 15 % °C
        WaterDensity(1,1) double {mustBeFinite, mustBeReal, mustBeNonnegative, mustBeNonsparse} = 0.5 % g/m^3
    end
    
    methods(Access = protected)
        function pl = pathlossOverDistance(pm, ~, tx, d, ~)
            pl = fogpl(d, tx.TransmitterFrequency, pm.Temperature, pm.WaterDensity)';
        end
        
        function lim = frequencyLimits(~)
            lim = [10e9 1000e9];
        end
    end
end