function varargout =  RectangularPlot3D(haxRadPat,X,Y,Z)
%RectangularPlot3D takes the X, Y and Z vectors (the 3 dimensions) and
%plots them in a 3D plot. It returns the axes of the current figure.

%   Copyright 2015-2020 The MathWorks, Inc.


surfHdl = surf(haxRadPat,X,Y,Z, 'FaceColor','interp');
set(surfHdl,'LineStyle','none','FaceAlpha',1.0);

colormap(haxRadPat,jet(256));
grid(haxRadPat,'on');
hold(haxRadPat,'off');
hfig = ancestor(haxRadPat,'figure');

z = zoom(hfig);
z.setAxes3DPanAndZoomStyle(haxRadPat,'camera');

if nargout == 1
    varargout{1} = surfHdl;
elseif nargout == 2
    varargout{1} = surfHdl;
    varargout{2} = haxRadPat;
end
end