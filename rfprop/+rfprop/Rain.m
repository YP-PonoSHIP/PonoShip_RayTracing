classdef Rain < rfprop.PropagationModel
%

% Copyright 2017-2024 The MathWorks, Inc.   

    properties
        RainRate(1,1) double {mustBeFinite, mustBeReal, mustBeNonnegative, mustBeNonsparse} = 16 % mm/h
        Tilt(1,1) double {mustBeFinite, mustBeReal, mustBeNonsparse} = 0 % °
    end
    
    methods(Access = protected)
        function pl = pathlossOverDistance(pm, ~, tx, d, el)
            
            % Pass elevation angle to rainpl, but use Tilt property for
            % polarization tilt angle, even though this is dependent on
            % direction. In the future, this could be automatically
            % computed from the transmitter antenna. Note that default tilt
            % of 0 typically produces "worst case" scenario to maximize
            % path loss.
            pl = rainpl(d,tx.TransmitterFrequency,pm.RainRate,el,pm.Tilt)';
        end
        
        function lim = frequencyLimits(~)
            lim = [1e9 1000e9];
        end
    end
end