function [yval, xval, legendstr, legendstruct] = radpattern2Ddata(MagE, ...
    phi, theta, frequency, slice, varargin)
%

%   Copyright 2015-2020 The MathWorks, Inc.

if strcmpi(slice, 'frequency') && isscalar(frequency)
    if isscalar(phi)
        yval      = reshape(MagE,numel(theta),1);
        xval      = (90-theta).';
        legendstr = strcat(num2str(phi),' deg');
        legendstruct.type  = 'az';
        legendstruct.value = phi;
        legendstruct.unit  = 'deg';        
        return;
    end
    
    if isscalar(theta)
        yval      = reshape(MagE, numel(phi),1);
        xval      = phi.';
        legendstr = strcat(num2str(90-theta),' deg');
        legendstruct.type  = 'el';
        legendstruct.value = 90-theta;
        legendstruct.unit  = 'deg';        
        return;
    end
end

if strcmpi(slice, 'azimuth')% phi variation
    
    yval      = reshape(MagE, numel(phi), numel(theta)).';
    legendstr = cell(size(phi));
    for m=1:numel(phi)
        legendstr{m}= strcat(num2str(phi(m)),' deg');
    end
    xval      = repmat((90-theta).', [1, length(phi)]);
    legendstruct.type  = 'az';
    legendstruct.value = phi;
    legendstruct.unit  = 'deg';
    return;
    
elseif strcmpi (slice, 'elevation')% theta variation
    
    yval = reshape(MagE, numel(phi), numel(theta));
    legendstr = cell(size(theta));
    for m=1:numel(theta)
        legendstr{m}= strcat(num2str(90-theta(m)),' deg');
    end
    xval     = repmat(phi.', [1, length(theta)]);
    legendstruct.type  = 'el';
    legendstruct.value = 90-theta;
    legendstruct.unit  = 'deg';
    return;
    
elseif strcmpi (slice, 'frequency')
    
    yval = MagE.';
    if isscalar(theta) && isscalar(phi)
        xval = phi;
    elseif isscalar(theta)
        xval     = repmat(phi.', [1, length(frequency)]);
    elseif isscalar(phi)
        xval     = repmat(90-theta.', [1, length(frequency)]);
    end
    legendstr = cell(size(frequency));
    [freqval,~,U]=engunits(frequency);
    for m=1:numel(frequency)
        legendstr{m}= sprintf('%d %sHz',freqval(m),U);%strcat(num2str(freqval(m)),U,' Hz');
    end
    legendstruct.type  = 'freq';
    legendstruct.value = freqval;
    legendstruct.unit  = strcat(U, 'Hz');
    return;
end

legendstr = '';
% This code is for patternCustom() function
if strcmpi(slice,'theta')
    
    SliceValue = varargin{1};
    theta1 = theta;
    phi1 = phi;
    
    thetarep = sort(repmat(theta1,[length(phi1),1]));
    phirep = repmat(phi1,[length(theta1),1]);
    
    MagE = reshape(MagE,[length(phi1)*length(theta1),1]);
    
    yval = zeros(length(phi1),length(SliceValue));
    xval = zeros(length(phi1),length(SliceValue));
    
    for i = 1:length(SliceValue)
        indstore = find(thetarep == SliceValue(i));
        if isempty(indstore)
            error(message('shared_channel:patternCustom:InvalidValueTheta'));
        else
            yval(:,i) = MagE(indstore);
            xval(:,i) = phirep(indstore);
        end
    end
    legendstruct.type  = 'theta';
    legendstruct.value = SliceValue;
    legendstruct.unit  = 'deg';
    
elseif strcmpi(slice,'phi')
    
    SliceValue = varargin{1};
    theta1 = theta;
    phi1 = phi;
    
    thetarep = sort(repmat(theta1,[length(phi1),1]));
    phirep = repmat(phi1,[length(theta1),1]);
    
    MagE = reshape(MagE,[length(phi1)*length(theta1),1]);
    
    yval = zeros(length(theta1),length(SliceValue));
    xval = zeros(length(theta1),length(SliceValue));
    
    for i = 1:length(SliceValue)
        indstore = find(phirep == SliceValue(i));
        if isempty(indstore)
            error(message('shared_channel:patternCustom:InvalidValuePhi'));
        else
            yval(:,i) = MagE(indstore);
            xval(:,i) = thetarep(indstore);
        end
    end
    legendstruct.type  = 'phi';
    legendstruct.value = SliceValue;
    legendstruct.unit  = 'deg';
end
end% of radpattern2Ddata