function tf = figureForwardState(hfig)
plotStruct  = getappdata(hfig,'PlotCustomiZationStructure');
if ~isempty(plotStruct)
    tf          = plotStruct.BringPlotFigureForward;
else
    tf          = true;
end
end