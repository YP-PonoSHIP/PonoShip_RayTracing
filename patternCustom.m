function varargout =  patternCustom(MagE, theta, phi,varargin)
%PATTERNCUSTOM  Plot 2D and 3D radiation pattern of the antenna in 
%               polar and rectangular coordinate system.
%
% patternCustom(MagE, theta, phi) plots the 2D and 3D radiation pattern
% over the specified phi and theta vectors.
%
% patternCustom(__, Name, Value) plots the 2D and 3D radiation pattern from
% the input data with additional options specified by the Name, Value
% pairs.
%
% hPlot = patternCustom(___) returns handles of the lines or surface in the 
% figure window, using any of the input arguments in the previous syntaxes.
%
% Input Arguments
%
% MagE: It can be either the Magnitude of the quantity plotted specified as
%       a vector or matrix, a phased.CustomAntennaElement object or a
%       phased.IsotropicAntennaElement object.
%       If the quantity is a vector, it should be of the same size
%       as phi and theta. If the quantity is a matrix, it should be of size 
%       phi x theta.
%
% theta: Theta angles in spherical coordinate system specified as a vector.
%        This needs to be provided only when MagE is a vector or matrix.
%
% phi:  Phi angles in the spherical coordinate system specified as a vector.
%       This needs to be provided only when MagE is a vector or matrix.
%       
%
% Below are the list of Name-Value pairs available in the patternCustom
% function:
%
% CoordinateSystem: Coordinate system to visualize the radiation pattern
% specified as polar | rectangular. The default is polar.
%
% Slice: Plane to visualize the 2D data specified as phi | theta. There are
% no defaults for this pair.
%
% SliceValue: Values for the Slice (phi or theta) specified as a vector.
% There are no defaults for this pair.
% 
% PatternOptions: The parameter is used to change the pattern plot properties. 
%                 The input should be a PatternPlotOptions object. When 
%                 PatternOptions are applied for patternCustom, two 
%                 properties of PatternPlotOptions object, Transparency and 
%                 MagnitudeScale can be changed. Other two proprties 
%                 SizeRatio, AntennaOffset which depend on the inset figure
%                 are ignored. It is not applicable for phased.CustomAntennaElement
%                 object and phased.IsotropicAntennaElement object.

%   Copyright 2019-2022 The MathWorks, Inc.


% Parsing through the inputs
parseobj = inputParser;
parseobj.FunctionName = 'patternCustom';

[varargin, hfig, ~] = parseParentArgs(varargin,"ValidParents",{'figure'}, "OnlyInput",false);

% Checking whether the MagE input is a vector or a matrix
if ~isobject(MagE)
    if nargin > 3
        [varargin{:}] = convertStringsToChars(varargin{:});
    end
    if isvector(MagE)
        narginchk(3,11);
        % Making sure the size of theta and phi is same
        if numel(theta) == numel(phi)
            typeValidationMat = @(x) validateattributes(x,{'numeric'},      ...
                {'vector','numel',numel(theta),'nonempty','real'}, 'patternCustom');
        else
            error(message('shared_channel:patternCustom:DimensionMismatch'));
        end
    elseif ismatrix(MagE)
        narginchk(3,11);
        typeValidationMat = @(x) validateattributes(x,{'numeric'},          ...
            {'size',[numel(phi),numel(theta)],'nonempty','real'},'patternCustom');
    end
else
    if nargin>1
        error(message('shared_channel:patternCustom:Support3DVizOnlyForPSTObj'));
    end
    phasedObj = MagE;
    
    if isa(phasedObj,'phased.CustomAntennaElement')||isa(phasedObj,'phased.IsotropicAntennaElement')
        narginchk(1,7);
        if ~isa(phasedObj, 'phased.CustomAntennaElement')
            phasedObj = phased.CustomAntennaElement; % Dummy object to get default isotropic response
        end
%         [MagE,az,el] = pattern(MagE,freqpat); % Dummy frequency
%         MagE = MagE';
        MagE = phasedObj.MagnitudePattern';
        if all(MagE==0)
            MagE = MagE + 1e-3;
        end
        phi = phasedObj.AzimuthAngles';
        theta = 90-phasedObj.ElevationAngles   ;
        typeValidationMat = @(x) validateattributes(x,{'numeric'},          ...
            {'size',[numel(phi),numel(theta)],'nonempty','real'},'patternCustom');
    else
        error(message('shared_channel:patternCustom:UnsupportedObject'));
    end    
end

expectedcoord = {'polar','rectangular'};
expectedslice = {'phi','theta'};

addRequired(parseobj,'MagE',typeValidationMat);
typeValidationVec = @(x) validateattributes(x,{'numeric'},              ...
    {'vector','nonempty', 'real','finite', 'nonnan'},'patternCustom');
addRequired(parseobj,'theta', typeValidationVec);
addRequired(parseobj,'phi',typeValidationVec);
addParameter(parseobj,'CoordinateSystem','polar',                   ...
    @(x) any(validatestring(x,expectedcoord)));
addParameter(parseobj,'Slice','',                   ...
    @(x) any(validatestring(x,expectedslice)));

isAntenaLicense = checkLicense;
if isAntenaLicense
    expectedpatternOptions = PatternPlotOptions;
    addParameter(parseobj,'PatternOptions',expectedpatternOptions,...
        @(x)validateoptions(x));
elseif any(strcmpi(varargin,'PatternOptions'))
    error(message('shared_channel:patternCustom:PatternOptionsReqAntenna'));
end

typeValidationSlice = @(x) validateattributes(x,{'numeric'},          ...
    {'vector','nonempty','real','finite', 'nonnan'}, 'patternCustom');
addParameter(parseobj,'SliceValue',[],typeValidationSlice);
parse(parseobj, MagE, theta, phi, varargin{:});

% Calling the function to get back the rearranged data
[MagEBlocks,theta1,phi1] = antennashared.internal.rearrangeData(MagE,theta,phi);

% Checking for the output argument
if nargout <= 1
    % This loop will be executed if 2D data is to be plotted
    if ((~isequal(parseobj.Results.Slice,'')) &&                        ...
            (~isequal(parseobj.Results.SliceValue,[])))
        
        % Calling the radpattern2Ddata function to get the data in the
        % appropriate format
        [MagStore,AngStore] = antennashared.internal.radpattern2Ddata(MagEBlocks,  ...
            phi1,theta1,1e9,parseobj.Results.Slice,                     ...
            parseobj.Results.SliceValue);
        
        % This loop is executed if the coordinatesystem specified is polar
        if strcmpi(parseobj.Results.CoordinateSystem,'polar')
            if isempty(hfig)
                hfig = gcf;
                ax = gca;
                if ~ishold(ax)
                    cla(ax);
                end
            else
                ax = hfig.CurrentAxes;
                if ~isempty(ax)&& ~ishold(ax)
                    clf(hfig);
                end               
                ax = axes(hfig);
            end
            hPlot = polarpattern(AngStore(:,1),MagStore,                ...
                'DrawGridToOrigin', 1, 'LineWidth', 2, 'GridWidth', 1.5,...
                'AngleResolution', 30, 'Parent',ax);
            if strcmpi(parseobj.Results.Slice,'theta')
                createLabels(hPlot, 'theta=%d#deg',                     ...
                    parseobj.Results.SliceValue(1:end));
            elseif strcmpi(parseobj.Results.Slice,'phi')
                createLabels(hPlot, 'phi=%d#deg',                       ...
                    parseobj.Results.SliceValue(1:end));
            end       
        % coordinatesystem specified is rectangular
        elseif strcmpi(parseobj.Results.CoordinateSystem,'rectangular')
            if isempty(hfig)
                hfig = gcf;
                ax = gca;
                if ~ishold(ax)
                    cla(ax);
                end
            else
                ax = hfig.CurrentAxes;
                if ~isempty(ax)&& ~ishold(ax)
                    clf(hfig);
                end
                ax = axes(hfig);
            end
            str = cell(length(parseobj.Results.SliceValue),1);            
            hPlot = plot(ax,AngStore,MagStore); % gca
            set(hPlot,'LineWidth',2);
            grid(ax, 'on');
            ylabel(ax, 'Magnitude');
            if strcmpi(parseobj.Results.Slice,'theta')
                xlabel(ax, 'Phi (degree)');
                for m = 1:length(parseobj.Results.SliceValue)
                    str{m}= strcat('theta =', num2str(parseobj.Results.SliceValue(m)),' deg');
                end
            elseif strcmpi(parseobj.Results.Slice,'phi')
                xlabel(ax, 'Theta (degree)');
                for m = 1:length(parseobj.Results.SliceValue)
                    str{m}= strcat('phi =', num2str(parseobj.Results.SliceValue(m)),' deg');
                end
            end
            if ishold
                str = {ax.Legend.String{1:end-1},str{:}};
            end
            legend(ax, str);
            legend(ax, 'Location','Best');
        end
        
        % This loop will be executed if 3D data is to be plotted
    elseif ((isequal(parseobj.Results.Slice,'')) &&                     ...
            (isequal(parseobj.Results.SliceValue,[])))        

        if isAntenaLicense
            patternOptions = parseobj.Results.PatternOptions;
            patternOptions.setMagnitude(MagEBlocks);
            patternOptions.setThetaPhi(theta1,phi1);

            % Modify the magnitude based on the magnitude scale from
            % PatternPlotOptions
            MagnitudeScale = patternOptions.MagnitudeScale;
            if ~(isempty(MagnitudeScale))
                minval = MagnitudeScale(1);
                maxval = MagnitudeScale(2);
                MagEBlocks(MagEBlocks<=minval) = minval;
                MagEBlocks(MagEBlocks>=maxval) = maxval;
            end
        end

        % This loop is executed if the coordinatesystem specified is polar
        if isempty(hfig)
            hfig = gcf;
            ax = gca;
            if ~ishold(ax)
                cla(ax);
            end
        else
            ax = hfig.CurrentAxes;
            if ~isempty(ax)&& ~ishold(ax)
                clf(hfig);
            end
            ax = axes(hfig);
        end
        if strcmpi(parseobj.Results.CoordinateSystem,'polar')            
            [hPlot] = antennashared.internal.radiationpattern3D(ax, MagEBlocks,theta1,phi1,'CurrentAxes', 1);
            axesHand = hPlot.Parent;
            axis(axesHand,'normal')
            axesHand.DataAspectRatio = [1 1 1];

            % This loop is executed if the coordinatesystem specified is
            % rectangular
        elseif strcmpi(parseobj.Results.CoordinateSystem,'rectangular')
            [thetaBlocks,phiBlocks] = meshgrid(theta1, phi1);
            [hPlot,axesHand] = antennashared.internal.RectangularPlot3D(ax, thetaBlocks, ...
                phiBlocks,MagEBlocks);
            axis(axesHand,'normal')
            axesHand.DataAspectRatio = [1 1 1];
            %set(axesHand,'Position',[0.1300 0.1100 0.7750 0.8150]);
            xlabel(ax, 'Theta (degree)');
            ylabel(ax, 'Phi (degree)');
            zlabel(ax, 'Magnitude');
        end

        if isAntenaLicense
            if strcmpi(parseobj.Results.CoordinateSystem,'polar')
                % Set transparency
                transparency = patternOptions.Transparency;
                sf = axesHand.Children(1);
                h = hggroup(ax, 'tag','patterngroup');
                set(axesHand.Children(2:end),'Parent',h);
                sf.FaceAlpha = transparency;
                % Apply the MagnitudeScale limit to CLim
                if ~(isempty(patternOptions.MagnitudeScale))
                    axesHand.CLim = patternOptions.MagnitudeScale;
                else
                    % Find the CData of Surface
                    maxlim = max(max(axesHand.Children.Children(end).CData));
                    minlim = min(min(axesHand.Children.Children(end).CData));
                    axesHand.CLim = [minlim maxlim];
                end
            elseif strcmpi(parseobj.Results.CoordinateSystem,'rectangular')
                % Set transparency
                transparency = patternOptions.Transparency;
                % Apply transparency on the new surface if hold on is
                % enabled.
                sf = axesHand.Children(1);
                set(sf,'FaceAlpha',transparency);
            end

            % Set tag to the axis
            set(axesHand,'tag','rectangular');
            %hfig = gcf;

            % Set tag to the plot
            taggrp = hggroup('Tag',num2str(rand*100),'Parent' ,axesHand,...
                'HandleVisibility','off');
            % Set tag to the PatternPlotOptions
            parseobj.Results.PatternOptions.setTag(taggrp);

            % Set current figure to the PatternPlotOptions
            parseobj.Results.PatternOptions.setPlot(hfig);
        end

        % These loops are executed if the number of input arguments are not enough
    elseif ((isequal(parseobj.Results.Slice,'')) && (~isequal(parseobj.Results.SliceValue,[])))
        error(message('shared_channel:patternCustom:UnspecifiedCut'));
    elseif (~(isequal(parseobj.Results.Slice,'')) && (isequal(parseobj.Results.SliceValue,[])))
        error(message('shared_channel:patternCustom:UnspecifiedCutValue'));
    end
    
    if nargout == 1
        varargout{1} = hPlot;
    end
elseif nargout > 1
    error(message('shared_channel:patternCustom:IncorrectNumArguments','output','output','1'));
end

end
function op = validateoptions(options)
    if(strcmpi(class(options) ,'PatternPlotOptions'))&& isscalar(options)
        op = true;
    else
        op = false;
        error(message('shared_channel:patternCustom:InvalidPatternPlotOptions'));
    end
end
function isLicense = checkLicense
if builtin('license', 'test', 'Antenna_Toolbox') && ~isempty(ver('antenna'))
    isLicense = true;
else
    isLicense = false;
end
end