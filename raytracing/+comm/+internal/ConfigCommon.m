classdef (Abstract, HandleCompatible) ConfigCommon 
%ConfigCommon Base object for configuration objects, for simulation and
%code generation

% Copyright 2018-2022 The MathWorks, Inc.

%#codegen

methods 
  function v = set(obj, prop)
    v = obj.([prop, '_Values']);
  end       
end

methods (Access = protected)  
  % Validate and return a matching string from the enum property options
  function y = validateEnumProperties(obj, prop, value)
    options = set(obj, prop);
    y = validatestring(value, options, mfilename, prop);
  end
  
  function flag = isInactiveProperty(~, ~)
    flag = false;
  end
  
  function [flag, msg] = isUndefinedProperty(~, ~)
    flag = false;
    msg = '';
  end

end

end

% [EOF]