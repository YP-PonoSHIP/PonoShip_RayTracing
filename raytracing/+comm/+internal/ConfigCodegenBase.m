classdef (Abstract) ConfigCodegenBase
%ConfigCodegenBase Base object for configuration object code generation

% Copyright 2018-2022 The MathWorks, Inc.

%#codegen

methods
    function obj = ConfigCodegenBase(varargin)
        obj = setProperties(obj,varargin{:});
    end
end

methods (Access = protected)
    function y = validateEnumProperties(obj, prop, value)
        options = obj.([prop, '_Values']);
        y = validatestring(value, options, mfilename, prop);
    end

    function obj = setProperties(obj,varargin)
        coder.internal.errorIf((mod(nargin-1,2) ~= 0),'shared_channel:ConfigBase:InvalidPVPairs');

        for i = 1:2:nargin-1
            obj.(varargin{i}) = varargin{i+1};
        end
    end
end
end

% [EOF]