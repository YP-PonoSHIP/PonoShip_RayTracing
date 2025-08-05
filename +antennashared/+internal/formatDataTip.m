
function txt = formatDataTip(src, evt, label, data, units, showXYZ)
    % This function can be used to format the data tip to look similar to
    % matlab data tip. The data tip will be rendered as following
    % label : data units
    % if showXYZ is true XYZ coordinates will be added to the list of
    % labels.
    
    % change interpreter to tex
    src.Interpreter = 'tex';
    % if showXYZ add XYZ and position
    if showXYZ
        label = [{'X','Y','Z'},label];
        data = [{evt.Position(1), evt.Position(2), evt.Position(3)}, data];
        units = [{'','',''},units];
    end
    % generate tooltip text based on labels and data
    txt = [];
    for i = 1:numel(label)
        txt = [txt '\color{black}\rm ' label{i} '\color[rgb]{0,0.6,1}\bf ' convertToString(data{i}) '\color{black}\rm ' units{i} newline];
    end
end
% 
function strval = convertToString(val)
    
    if ischar(val)
        strval = val;
    elseif isstring(val)
        strval = char(val);
    elseif isnumeric(val)
        if isscalar(val)
            strval = num2str(val);
        else
            strval = mat2str(val);
        end
    end
end