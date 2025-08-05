function [remArgs,ParentHandle,ParentType] = parseParentArgs(Args, opts )
    % This function is used to filter the parent name value pair from a
    % cell array of name value pair arguments (varargin).
    % Inputs:
    % Args          : cell array of name value pair arguments.
    % validParents  : cell array of strings with type of valid graphical
    %                 parents. The strings should belong to the following
    %                 list . "figure","axes","uiaxes","panel","container"

    arguments
        Args                (1,:) cell 
        opts.OnlyInput      (1,1) logical = false
        opts.ValidParents   (1,:) cell = {'figure','axes','uiaxes','uipanel','uicontainer'}
    end
    
    OnlyInput = opts.OnlyInput;
    ValidParents = opts.ValidParents;

    if isempty(Args)
        % assign outputs
        remArgs = Args;
        ParentHandle = [];
        ParentType = [];
        return;
    end

    classCell = convertValidParentsToClassName(ValidParents);
    isParentPresent = 0;
    for i = 1:numel(Args)
        if strcmpi(Args{i},'hParent')
            isParentPresent = 1;
            break;
        end
    end
    if isParentPresent
        
        value = Args{i+1};
        if ~isgraphics(value)
            error(message('shared_channel:shared_channel:InvalidHandle'));
        end
        validateattributes(value,classCell,"scalar",'','hParent',i);
        % remove the parent name value pair args
        Args(i:i+1) = [];

        % assign outputs
        remArgs = Args;
        ParentHandle = value;
        ParentType = value.Type;
    else
        if OnlyInput
            % classCell = convertValidParentsToClassName(ValidParents);
            % value = Args{1};
            % validateattributes(value,classCell,"scalar",'','hParent',i);
            % % remove the parent name value pair args
            % Args(1) = [];
            % 
            % % assign outputs
            % remArgs = [];
            % ParentHandle = value;
            % ParentType = value.Type;
        else
            % assign outputs
            remArgs = Args;
            ParentHandle = [];
            ParentType = [];
        end
        
    end
% if OnlyInput
%     p = inputParser;
%     parse(p,Args{:});
% end

end

function classCell = convertValidParentsToClassName(validParents)
    classCell = cell(1,numel(validParents));
    for i = 1:numel(validParents)
        validatestring(validParents{i},{'figure','axes','uiaxes','uipanel','uicontainer'});
        switch validParents{i}
            case 'figure'
                classCell{i} = 'matlab.ui.Figure';
            case 'axes'
                classCell{i} = 'matlab.graphics.axis.Axes';
            case 'uiaxes'
                classCell{i} = 'matlab.ui.control.UIAxes';
            case 'uipanel'
                classCell{i} = 'matlab.ui.container.Panel';
            case 'uicontainer'
                classCell{i} = 'matlab.ui.container.internal.UIContainer';
        end
    end
end
