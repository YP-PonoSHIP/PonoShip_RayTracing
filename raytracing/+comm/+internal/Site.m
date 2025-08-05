classdef Site < comm.internal.ConfigBase 
%Site Transmitter or receiver site for internal raytrace function
%   SITE = COMM.INTERNAL.SITE creates a transmitter site or receiver site
%   object, SITE. This object specifies the configuration for an antenna
%   site.
% 
%   SITE = COMM.INTERNAL.SITE(TXSITE/RXSITE) creates a site object, SITE,
%   by retrieving relative information from a txsite/rxsite object.
%
%   SITE = COMM.INTERNAL.SITE(Name,Value) creates a site object, SITE, with
%   the specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as (Name1,Value1,
%   ...,NameN,ValueN).
%
%   COMM.INTERNAL.SITE properties:
%
%   Position        - Antenna [x;y;z] position (m)
%   Antenna         - Antenna element or array object
%   OrientationAxes - Antenna orientation axes
%   Frequency       - Transmitter operating frequency (Hz) 
%
%   % Example 1: 
%   %   Create a transmitter site with carrier frequency of 28e9 Hz. 
%
%   site = comm.internal.Site( ...
%       "Position", ones(3,1), ...
%       "OrientationAxes", fliplr(eye(3)), ...
%       "Frequency", 28e9);
%
%   See also txsite, rxsite.

% Copyright 2020-2024 The MathWorks, Inc.

properties
    Position = zeros(3, 1)
    Antenna = 'isotropic'
    OrientationAxes = eye(3)
    Frequency = [];
end

methods
  function obj = Site(varargin)
    if (nargin == 1)
        % For a txsite/rxsite array input, it has been ensured that all the
        % objects in the array have the same CoordinateSystem setting.
         coder.internal.errorIf( ...
            ~isa(varargin{1}, 'rfprop.AntennaSite') || ...
            ~strcmp(varargin{1}(1).CoordinateSystem, 'cartesian'), ...
            'shared_channel:Site:InvalidSiteInput');
        txrxs = varargin{1};
        numSites = numel(txrxs);
        obj(size(txrxs, 1), size(txrxs, 2)) = obj;
        for i = 1:numSites
            obj(i).Position = txrxs(i).AntennaPosition;
            obj(i).Antenna = txrxs(i).Antenna;
            obj(i).OrientationAxes = ...
                rfprop.internal.coordinateTransformationMatrix( ...
                txrxs(i).AntennaAngle);
            if isa(txrxs, 'txsite')
                obj(i).Frequency = txrxs(i).TransmitterFrequency;
            end
        end
    else
        obj = setProperties(obj,varargin{:});
    end
  end
  
  function obj = set.Position(obj, pos)
    validateattributes(pos, {'double'}, ...
        {'real','finite','size',[3 1]}, ...
        [class(obj) '.' 'Position'], 'Position');    
    obj.Position = pos;
  end
  
  function obj = set.Antenna(obj, ant)
    if ~((ischar(ant) || isstring(ant)) && strcmp(ant, 'isotropic'))
        validateattributes(ant, {'arrayConfig', ...
            'phased.internal.AbstractAntennaElement', ...
            'phased.internal.AbstractArray', ...
            'phased.internal.AbstractSubarray', ...
            'em.Antenna', 'em.Array', 'installedAntenna', 'customAntennaStl'}, {'scalar'}, ...
            [class(obj) '.' 'Antenna'], 'Antenna');
    end
    obj.Antenna = ant;
  end
  
  function obj = set.OrientationAxes(obj, axes)
    % Validate orientation axes to be a 3-by-3 unitary matrix
    validateattributes(axes, {'double'}, ...
        {'real','finite','size',[3 3]}, ...
        [class(obj) '.' 'OrientationAxes'], 'OrientationAxes');
    coder.internal.errorIf(...
        max(max(abs(axes'*axes-eye(3)))) > sqrt(eps), ...
        'shared_channel:Site:AxesNotUnitary');
    
    obj.OrientationAxes = axes;
  end

  function obj = set.Frequency(obj, freq)
    if ~isempty(freq)
        validateattributes(freq, {'double'}, ...
            {'real','positive','finite','scalar'}, ...
            [class(obj) '.' 'Frequency'], 'Frequency');
    end
    obj.Frequency = freq;
  end  
end

end