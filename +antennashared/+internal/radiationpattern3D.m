function varargout =  radiationpattern3D(ax,MagE1,theta1,phi1, varargin) 
%RADIATIONPATTERN3D plots 3D radiation pattern for the antenna
%
%   The function plots the 3D radiation pattern of the antenna 
%
%   OBJ is the antenna object
%
%   FREQ is a single frequency point at which the radiation pattern is
%   plotted

%   Copyright 2015-2020 The MathWorks, Inc.

% Compute the max and min values of directivity
minval = min(min(MagE1));
maxval = max(max(MagE1));

parserObj          = inputParser();
addParameter(parserObj,'offset',minval);
addParameter(parserObj,'spherical',false);
addParameter(parserObj,'CurrentAxes',0);
addParameter(parserObj,'plottype','');
addParameter(parserObj,'IndBeam',[]);

parse(parserObj,varargin{:});
temp = parserObj.Results.CurrentAxes;
minval = parserObj.Results.offset;
if any(strcmpi(varargin, 'IndBeam'))
    antoff = parserObj.Results.IndBeam{1};
    surfnum = parserObj.Results.IndBeam{2};
    tabn = parserObj.Results.IndBeam{3};
else
    antoff = [0 0 0];
    surfnum = 1;
    tabn = {};
end
if(parserObj.Results.spherical)
    r = ones(length(phi1),length(theta1));
elseif isequal(unique(MagE1),minval)        % The isotropic case
    r = MagE1;
else
    r = MagE1 - minval;
end

if builtin('license', 'test', 'Antenna_Toolbox') && ~isempty(ver('antenna'))
    [~, units] = em.FieldAnalysisWithFeed.getfieldlabels(parserObj.Results.plottype);
else
    units       = '';
end

if ~temp
    % Create the axis for the radiation pattern
    
    haxRadPat = ax;
    
    hfig = ancestor(haxRadPat,'figure'); %hfig = get(haxRadPat, 'Parent');
    if(~matlab.ui.internal.isUIFigure(hfig)) % disabling figure toolbar for webfigure and uifigure 
        % FIXME: switching the toolbar off for uifigure until g3072859 is fixed 
        set(hfig,'toolbar','figure'); % Show standard toolbar %none
    end
%     cameratoolbar(hfig,'hide'); % Hide camera toolbar    
    tempval = 0; 
else
    haxRadPat = ax;
    hfig = ancestor(haxRadPat,'figure');
    tempval = 1;
end

hold(haxRadPat,'on');

add_x_y_z_labels(haxRadPat,antoff);
add_az_el_labels(haxRadPat,antoff);
draw_circle( haxRadPat, 90, 1:1:360,'--mw-graphics-colorSpace-rgb-blue',antoff,'XY_Circle'); % Circle in the x-y plane
draw_circle( haxRadPat, 1:1:360, 0,'--mw-graphics-colorSpace-rgb-green',antoff,'XZ_Circle'); % Circle in the x-z plane
draw_circle( haxRadPat, 1:1:360,90, '--mw-graphics-colorSpace-rgb-red',antoff,'YZ_Circle'); % Circle in the y-z plane

surfHdl = draw_3d_plot(haxRadPat, MagE1,theta1,phi1, r, tempval,units,antoff,surfnum);
dcm = datacursormode(hfig); %gcf

%set(dcm,"UpdateFcn",@(tipObj,evt)createDataTip(tipObj,evt,maxval,minval,units));
z = zoom(hfig);
z.setAxes3DPanAndZoomStyle(haxRadPat,'camera');

% This is for patternCustom function.
if nargout == 1
    varargout{1} = surfHdl;
end

% Add colorbar
% Till geck G1115830 is fixed

cbar = colorbar('peer',haxRadPat);
ylabel(cbar, units);
view(haxRadPat,135,20);
hold(haxRadPat,'off');
try
    set(hfig, 'NextPlot', 'replace');
    if antennashared.internal.figureForwardState(hfig) &&~calc_ind_beam
        figure(hfig); %shg;
    end
catch
end

end% of radiationpattern3D

function surfHdl =  draw_3d_plot(axes1, MagE1,theta1,phi1,r1, val,units,antoff,surfnum)

[theta,phi] = meshgrid(theta1, phi1);
MagE        = reshape(MagE1,length(phi1),length(theta1));
r           = reshape(r1,length(phi1),length(theta1));
[X, Y, Z]   =  antennashared.internal.sph2cart(phi, theta, r./max(max(r)));

X=X+antoff(1);
Y=Y+antoff(2);
Z=Z+antoff(3);

surfHdl = surf(axes1,X,Y,Z,MagE, 'FaceColor','interp');
set(surfHdl,'LineStyle','none','FaceAlpha',1.0,'Tag','3D polar plot');% change to 0.9
createDataTip(surfHdl,X,Y,Z,MagE,units);
if ~val   
    set(axes1, 'Position',[0.2 0.2 0.7 0.7], 'DataAspectRatio',[1 1 1]); %[0.28 0.24 0.7 0.7] %[0.05 0.05 0.9 .85]
    axis(axes1,'vis3d');
    axis(axes1,[-1.2 1.2 -1.2 1.2 -1.2 1.2]);
    axis (axes1,'off');
else
    %set(axes1,'DataAspectRatio',[1 1 1]);
    axis(axes1,'vis3d');   
    axis(axes1,[-1.2 1.2 -1.2 1.2 -1.2 1.2]);
    axis(axes1,'off');
    axis(axes1,'equal');
end
colormap(axes1,jet(256));
surfHdl.UserData.SurfaceNum = surfnum;
end% of draw_3d_plot

function draw_circle (axes1, theta, phi, color,antoff,Tag)

import matlab.graphics.internal.themes.specifyThemePropertyMappings
[theta,phi] = meshgrid(theta, phi);
[X, Y, Z]   = antennashared.internal.sph2cart(phi, theta, 1.1);
X=X+antoff(1);
Y=Y+antoff(2);
Z=Z+antoff(3);
p = plot3(axes1,X,Y,Z,'LineWidth',2,'Tag',Tag);
specifyThemePropertyMappings(p,'Color',color);

end % of draw_circle

function add_x_y_z_labels(axes1,antoff)
% Create pseudo-axes and x/y/z mark ticks

r      = 1.2;
XPos = r;
YPos = r;
ZPos = r;
plot3( axes1, [antoff(1),XPos+antoff(1)],[antoff(2),antoff(2)],[antoff(3),antoff(3)],'r','LineWidth',1.5 );
text(axes1,1.1*XPos+antoff(1),antoff(2),antoff(3), 'x');
plot3( axes1, [antoff(1),antoff(1)],[antoff(2),YPos+antoff(2)],[antoff(3),antoff(3)],'g','LineWidth',1.5 );
text(axes1,antoff(1),1.05*YPos+antoff(2),antoff(3), 'y');
plot3( axes1, [antoff(1),antoff(1)],[antoff(2),antoff(2)],[antoff(3),ZPos+antoff(3)],'b','LineWidth',1.5 );
text(axes1,antoff(1),antoff(2),1.05*ZPos+antoff(3), 'z');
end% of add_x_y_z_labels

function add_az_el_labels(axes1,antoff)
% Display azimuth/elevation

% Create arrows to show azimuth and elevation variation
XPos = 1.15;
%ZPos = 1.15;
draw_arrow(axes1,[XPos+antoff(1) antoff(2)],[XPos+antoff(1) 0.1+antoff(2)],1.5,antoff(3),'xy');
text(axes1,1.2+antoff(1),0.12+antoff(2),0.0+antoff(3), texlabel('az'));
draw_arrow(axes1,[XPos+antoff(1) antoff(3)],[XPos+antoff(1) 0.1+antoff(3)],1.5,antoff(2), 'xz');
text(axes1,1.2+antoff(1),-0.025+antoff(2),0.15+antoff(3), texlabel('el'));
%draw_arrow([0 ZPos],[0.1 ZPos],1.5, 0, 'xz');
%text(0.16,0,1.2, texlabel('el'));

end% of add_az_el_labels



function rtn = createDataTip(surf,X,Y,Z,MagE,units)
    if strcmpi(surf.Tag,'3D polar plot') && strcmpi(surf.Type,'surface')
        [az,el,r] = cart2sph(X,Y,Z);
        az180 = az.*(180/pi);
        el180 = el.*(180/pi);
        magnitude = MagE;
        if ~any(strcmpi(units,{'dBi','deg',''}))
            % if pattern units are not dBi or deg or none donot add eng
            % units
            [magnitude,fact,u] = engunits(magnitude);
            units = [u units];
        end
        
        row = dataTipTextRow('Az',az180,"%.4g deg");
        surf.DataTipTemplate.DataTipRows(1) = row;

        row = dataTipTextRow('El',el180,"%.4g deg");
        surf.DataTipTemplate.DataTipRows(2) = row;

        formatString = "%.4g " + units;
        row = dataTipTextRow('Mag',magnitude,formatString);
        surf.DataTipTemplate.DataTipRows(3) = row;

    end
end
% function rtn = createDataTip(src,evt,maxval,minval,units)
%     if strcmpi(evt.Target.Tag,'3D polar plot') && strcmpi(evt.Target.Type,'surface')
%         Position = evt.Position;
%         [az,el,r] = cart2sph(Position(1),Position(2),Position(3));
%         % magnitude = (maxval - minval)*r + minval;
%         CDataRows = size(evt.Target.CData,1);
%         colIndex = ceil(evt.DataIndex/CDataRows);
%         rowIndex = floor(evt.DataIndex - floor(evt.DataIndex/ CDataRows)*CDataRows);
%         magnitude = evt.Target.CData(rowIndex,colIndex);
%         az180 = az*(180/pi);
%         el180 = el*(180/pi);
% %         rtn = ['Az: ' , num2str(az180), ' deg' , newline , 'El: ' , num2str(el180) , ...
% %              ' deg',newline,'Mag: ' , num2str(round(magnitude,4)) , ' ',units];
%         if ~any(strcmpi(units,{'dBi','deg',''}))
%             % if pattern units are not dBi or deg or none donot add eng
%             % units
%             [magnitude,fact,u] = engunits(magnitude);
%             units = [u units];
%         end
%         rtn = antennashared.internal.formatDataTip(src,evt,{'Az','El','Mag'},...
%             {az180, el180, round(magnitude,4)} , {'deg','deg',units}, 0);
%     else
%         src.delete;
%         rtn = '';
%     end
% end
function draw_arrow(axes1,startpoint,endpoint,headsize, offset, plane)


if nargin == 4
    plane = 'xy';
    offset= 0;
end

v1 = headsize*(startpoint-endpoint)/2.5;

theta      = 22.5*pi/180;
theta1     = -1*22.5*pi/180;
rotMatrix  = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
rotMatrix1 = [cos(theta1) -sin(theta1) ; sin(theta1) cos(theta1)];

v2 = v1*rotMatrix;
v3 = v1*rotMatrix1;
x1 = endpoint;
x2 = x1 + v2;
x3 = x1 + v3;
if strcmpi(plane, 'xy')
    fill3(axes1,[x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],[offset offset offset],'k');    
    plot3(axes1,[startpoint(1) endpoint(1)],[startpoint(2) endpoint(2)],      ...
        [offset offset],'linewidth',1.5,'color','k');
elseif strcmpi(plane,'xz')    
    fill3(axes1,[x1(1) x2(1) x3(1)],[offset offset offset],[x1(2) x2(2) x3(2)],'k');    
    plot3(axes1,[startpoint(1) endpoint(1)],[offset offset],                  ...
        [startpoint(2) endpoint(2)],'linewidth',1.5,'color','k');
elseif strcmpi(plane,'yz')    
    fill3(axes1,[offset offset offset],[x1(1) x2(1) x3(1)],[x1(2) x2(2) x3(2)],'k');    
    plot3(axes1,[offset offset],[startpoint(1) endpoint(1)],                  ...
        [startpoint(2) endpoint(2)],'linewidth',1.5,'color','k');
end
end% of draw arrow