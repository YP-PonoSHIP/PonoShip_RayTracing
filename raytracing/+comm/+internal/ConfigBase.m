classdef (Abstract, HandleCompatible) ConfigBase < comm.internal.ConfigCommon & matlab.mixin.CustomDisplay 
%ConfigBase Base object for configuration objects

% Copyright 2018-2022 The MathWorks, Inc.

%#codegen

methods
  function obj = ConfigBase(varargin)
    obj = setProperties(obj,varargin{:});
  end
end

methods (Access = protected)
  function obj = setProperties(obj,varargin)
    coder.internal.errorIf((mod(nargin-1,2) ~= 0),'shared_channel:ConfigBase:InvalidPVPairs');

    for i = 1:2:nargin-1
      obj.(varargin{i}) = varargin{i+1};
    end
  end

  function groups = getPropertyGroups(obj)
    % Returns three property groups for display; normal, read-only, and
    % constant
    propName  = properties(obj);
    mc = meta.class.fromName(class(obj)); % Get metaclass
    allProperties = {mc.PropertyList.Name};
    % If defined explicitly, consider the order of properties as per
    % 'CustomPropList'
    if any(strcmpi(allProperties,'CustomPropList'))
      propName = obj.CustomPropList;
    end
    propVal = cell(numel(propName),1);
    activeIdx = true(size(propName));
    readOnlyIdx = false(size(propName));
    constIdx = false(size(propName));
    mc = meta.class.fromName(class(obj)); % Get metaclass
    % For each property...
    for n = 1:numel(propName)
      if isInactiveProperty(obj,propName{n})
        % If it is inactive then do not add it to the property list
        activeIdx(n) = false;
      elseif isUndefinedProperty(obj,propName{n})
        % If it is undefined then set the value to an empty for disp()
        propVal{n} = [];
      else
        propVal{n} = obj.(propName{n});
      end
      % Determine if property is read-only or constant from attributes
      id = cellfun(@(x)strcmp(x,propName{n}),{mc.PropertyList.Name});
      % If a read-only property then display it in a separate property group
      readOnlyIdx(n) = strcmp(mc.PropertyList(id).SetAccess,'private') && strcmp(mc.PropertyList(id).GetAccess,'public');
      % If a constant property then display it in a separate property group
      constIdx(n) = mc.PropertyList(id).Constant;
    end
    % Create three property lists, one of read-only properties, one of
    % normal properties, and one of constant properties
    normalPropList = cell2struct(propVal(activeIdx&~readOnlyIdx&~constIdx),propName(activeIdx&~readOnlyIdx&~constIdx));
    % Display non-empty lists as groups
    groups = matlab.mixin.util.PropertyGroup(normalPropList);
    if (any(activeIdx&readOnlyIdx))
      readOnlyPropList = cell2struct(propVal(activeIdx&readOnlyIdx),propName(activeIdx&readOnlyIdx));
      groups = [groups matlab.mixin.util.PropertyGroup(readOnlyPropList,getString(message('shared_channel:ConfigBase:ROProperties')))];
    end
    if (any(activeIdx&constIdx))
      constPropList = cell2struct(propVal(activeIdx&constIdx),propName(activeIdx&constIdx));
      groups = [groups matlab.mixin.util.PropertyGroup(constPropList,getString(message('shared_channel:ConfigBase:ConstProperties')))];
    end
  end
  
  function s = getFooter(obj)
    % If any of the properties are invalid, display the message in the
    % footer. Create a new line for each undefined property
    str = cell(0,1);
    i = 1;
    propName = properties(obj);
    % For all properties...
    for n = 1:numel(propName)
        [isUndefined, msg] = isUndefinedProperty(obj,propName{n});
        if isUndefined
            str{i} = msg;
            i = i+1;
        end
    end
    s = char(str); % Convert the cell array to a character array; each cell is a new line
    if ~isempty(s) % Wrap to command window width if there is something to print
        s = matlab.internal.display.printWrapped(char(str));
    end
  end
  
  function displayNonScalarObject(obj)
    % If an object array, display the size and type of object (with a link
    % to the class doc)
    fprintf('  %s %s %s\n\n',matlab.mixin.CustomDisplay.convertDimensionsToString(obj), ...
        matlab.mixin.CustomDisplay.getClassNameForHeader(obj), ...
        getString(message('shared_channel:ConfigBase:Array')));
  end
  
  function displayEmptyObject(obj)
    % If an empty object array, treat as an object array for display
    displayNonScalarObject(obj);
  end
end

methods(Access = private, Static)
  function name = matlabCodegenRedirect(~)
    name = 'comm.internal.ConfigCodegenBase';
  end
end

end
