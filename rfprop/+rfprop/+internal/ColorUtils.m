classdef (Sealed, Hidden) ColorUtils
    %ColorUtils   Reserved for MathWorks internal use only
    
    %   Copyright 2019 The MathWorks, Inc. 
    

    methods(Static)
             
        function colorStrings = colors(flag)
            %colors   Color values
            
            if nargin && strcmp(flag,'long') % Long values only
                colorStrings = {'yellow','magenta','cyan','red','green', ...
                    'blue','white','black'};
            else
                colorStrings = {'y','yellow','m','magenta','c','cyan', ...
                    'r','red','g','green','b','blue','w','white','k','black'};
            end
        end

        function hex = rgb2css(rgb)
            
            
            hex = ['#', ...
                dec2hex(round(255*rgb(1)),2), ...
                dec2hex(round(255*rgb(2)),2), ...
                dec2hex(round(255*rgb(3)),2)];
        end
        
        function rgb = str2rgb(colorspec)
            %str2rgb   RGB value for color specification
            
            switch lower(colorspec)
                case {'y','yellow'}
                    rgb = [1 1 0];
                case {'m','magenta'}
                    rgb = [1 0 1];
                case {'c','cyan'}
                    rgb = [0 1 1];
                case {'r','red'}
                    rgb = [1 0 0];
                case {'g','green'}
                    rgb = [0 1 0];
                case {'b','blue'}
                    rgb = [0 0 1];
                case {'w','white'}
                    rgb = [1 1 1];
                case {'k','black'}
                    rgb = [0 0 0];
                otherwise
                    error(message('shared_channel:rfprop:InvalidColorSpec', colorspec));
            end
        end
        
        function dataRGB = colorcode(data, cmap, clim)
            %colorcode   Color code data using colormap
            
            % Convert data to indices in colormap using color limits. Code
            % adapted from Tip on caxis reference page
            cmin = clim(1);
            cmax = clim(2);
            m = size(cmap, 1);
            colormapInd = fix((data-cmin)/(cmax-cmin)*m)+1;
            
            % Clamp values outside the range [1 m]
            colormapInd(colormapInd < 1) = 1;
            colormapInd(colormapInd > m) = m;
            
            % Get RGB from colormap
            dataRGB = ind2rgb(colormapInd, cmap);
        end
        
        function [legendColors, legendColorValues] = colormaplegend(cmap, clim)
            %legend   Colormap legend colors and values
            
            cmin = clim(1);
            cmax = clim(2);
            m = size(cmap, 1);

            % Compute legend values. The levels in the legend are generated using
            % ColorLimits. The strategy is to show the color limits and
            % intermediate values at a fixed step size, where the step size is
            % chosen to be intuitive and generate at least 10 levels. A default
            % step size of 10 is used. If this fails to produce at least 10 levels,
            % then a step size of 5 is used. If that also fails to produce at least
            % 10 levels, a step size of 2 is used. If that also fails, a final step
            % size of 1 is used.

            % Get rounded color limits
            colorMax = ceil(cmax);
            colorMin = floor(cmin);
            
            % Compute step size for intermediate color values
            colorStrengthRange = colorMax - colorMin;
            if (colorStrengthRange/10) > 20
                colorStep = floor(colorStrengthRange/20);
            elseif (colorStrengthRange/10) >= 9
                colorStep = 10;
            elseif (colorStrengthRange/5) >= 9
                colorStep = 5;
            elseif (colorStrengthRange/2) >= 9
                colorStep = 2;
            else
                colorStep = 1;
            end
            
            % Compute color strengths, ensuring that max is included
            colorStrengths = colorMin:colorStep:colorMax;
            if colorStrengths(end) ~= colorMax
                colorStrengths = [colorStrengths, colorMax];
            end
            
            % Limit legend values to 3 significant digits
            colorStrengths = round(colorStrengths,3,'significant');
            
            % Grow legend values
            numColors = numel(colorStrengths);
            legendColors = strings(1,numColors);
            legendColorValues = strings(1,numColors);
            for colorInd = 1:numColors
                colorStrength = colorStrengths(colorInd);
                legendInd = fix((colorStrength-cmin)/(cmax-cmin)*m)+1;
                rgb = ind2rgb(legendInd, cmap);
                legendColors(colorInd) = rfprop.internal.ColorUtils.rgb2css(rgb(:));
                legendColorValues(colorInd) = mat2str(round(colorStrength));
            end
            legendColors = fliplr(legendColors);
            legendColorValues = fliplr(legendColorValues);
        end
        
        function [dataRGB, legendCSSColors, legendColorValues] = dataColors(data, useColors, colors, cmap, clim, strengths, showLegend)
            %dataColors   Return data and legend colors using color and options
            
            legendCSSColors = string([]);
            legendColorValues = string([]);
            
            % Use Colors if user specified or if single signal level and no Colormap
            % specified
            if useColors
                % Initialize points RGB vectors
                imgR = zeros(size(data));
                imgG = imgR;
                imgB = imgR;
                
                % Color each user-specified level in the points color vector
                numColors = size(colors, 1);
                colorInd = 1;
                for k = 1:numel(strengths)
                    levelInd = (data == strengths(k));
                    
                    % Get level's color and assign in point RGB Vectors
                    color = colors(colorInd,:);
                    imgR(levelInd) = color(1);
                    imgG(levelInd) = color(2);
                    imgB(levelInd) = color(3);
                    
                    % Cycle through colors
                    colorInd = colorInd + 1;
                    if (colorInd > numColors)
                        colorInd = 1;
                    end
                    
                    % Grow legend values
                    if showLegend
                        legendCSSColors(end+1) = rfprop.internal.ColorUtils.rgb2css(color); %#ok<*AGROW>
                        colorStrength = strengths(k);
                        if (floor(colorStrength) == colorStrength)
                            numDigits = 0; % Show integer value
                        else
                            numDigits = 1; % Show one decimal place
                        end
                        legendColorValues(end+1) = mat2str(round(colorStrength,numDigits));
                    end
                end
                
                % Sort legend values so that legend descends from top to bottom
                if showLegend
                    [~,legendInd] = sort(strengths,'descend');
                    legendCSSColors = legendCSSColors(legendInd);
                    legendColorValues = legendColorValues(legendInd);
                end
                
                % Get RGB matrix
                dataRGB = cat(3, imgR, imgG, imgB);
            else
                dataRGB = rfprop.internal.ColorUtils.colorcode(data, cmap, clim);
                if showLegend
                    [legendCSSColors, legendColorValues] = rfprop.internal.ColorUtils.colormaplegend(cmap, clim);
                end
            end
        end
        
        function y = srgb2lin(x)
            gamma = cast(2.4,'like',x);
            a     = cast(1/1.055,'like',x);
            b     = cast(0.055/1.055,'like',x);
            c     = cast(1/12.92,'like',x);
            d     = cast(0.04045,'like',x);
            
            in_sign = -2 * (x < 0) + 1;
            x = abs(x);
            
            lin_range = (x < d);
            gamma_range = ~lin_range;
            
            y = zeros(size(x),'like',x);
            
            y(gamma_range) = exp(gamma .* log(a * x(gamma_range) + b));
            y(lin_range) = c * x(lin_range);
            
            y = y .* in_sign;
        end
    end
end
