classdef Gas < rfprop.PropagationModel
%

% Copyright 2017-2024 The MathWorks, Inc.  
        
    properties
        Temperature(1,1) double {mustBeFinite, mustBeReal, mustBeNonsparse} = 15 % °C
        AirPressure(1,1) double {mustBeFinite, mustBeReal, mustBeNonnegative, mustBeNonsparse} = 101300 % Pa
        WaterDensity(1,1) double {mustBeFinite, mustBeReal, mustBeNonnegative, mustBeNonsparse} = 7.5 % g/m^3
    end
    
    methods(Access = protected)
        function pl = pathlossOverDistance(pm, ~, tx, d, ~)
            pl = gaspl(d, tx.TransmitterFrequency, pm.Temperature, pm.AirPressure, pm.WaterDensity)';
        end
      
        function lim = frequencyLimits(~)
            lim = [1e9 1000e9];
        end
    end
end