function pm = propagationModel(pmName, varargin)
%

% Copyright 2017-2024 The MathWorks, Inc.  

pmName = rfprop.PropagationModel.validateName(pmName, 'propagationModel');

switch (pmName)
    case 'freespace'
        pm = rfprop.FreeSpace(varargin{:});
    case 'rain'
        pm = rfprop.Rain(varargin{:});
    case 'fog'
        pm = rfprop.Fog(varargin{:});
    case 'gas'
        pm = rfprop.Gas(varargin{:});
    case 'close-in'
        pm = rfprop.CloseIn(varargin{:});
    case 'longley-rice'
        pm = rfprop.LongleyRice(varargin{:});
    case 'tirem'
        pm = rfprop.TIREM(varargin{:});
    case 'raytracing'
        pm = rfprop.RayTracing(varargin{:});
end