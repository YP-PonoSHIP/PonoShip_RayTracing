function sharpEdges = getSharpEdges(TRI, filterAngle)
%GETSHARPEDGES Flag each edge is sharp or not
%   SHARPEDGES = COMM.INTERNAL.GEOMETRY.GETSHARPEDGES(TRI) finds each edge
%   of every triangle in the triangulation object, TRI, is sharp or not.
%   Any edge whose adjacent triangles have a dihedral angle that deviates
%   from 180 degrees by an angle greater than 2 degrees is considered to
%   be sharp. Edges that are shared by only one triangle, and edges that
%   are shared by more than two triangles are considered to be sharp.
%
%   SHARPEDGES is a logical matrix with the same size as
%   TRI.ConnectiviltyList to indicate sharp or not for the corresponding
%   triangle edges in TRI. The first/second/third column of SHARPEDGES
%   corresponds to the edge from vertex 1/2/3 to vertex 2/3/1 in the
%   corresponding triangle.
%   
%   SHARPEDGES = COMM.INTERNAL.GEOMETRY.GETSHARPEDGES(TRI, FILTERANGLE)
%   specifies the filter angle, FILTERANGLE, in degrees to detect sharp
%   edges. FILTERANGLE is 2 degrees by default.
%
%   See also triangulation/featureEdges

% Copyright 2020-2021 The MathWorks, Inc.

if nargin == 1
    filterAngle = 2; % 2 degrees by default
end

fe = featureEdges(TRI, filterAngle*pi/180);
% Double the feature edges to cover reverse directions
fe = [fe; fliplr(fe)]; 

% Formulate the logical matrix with the same as TRI.ConnectiviltyList to
% indicate each edge of every triangle is sharp or not
%   Column 1 is for edge from vertex 1 to 2 in the corresponding triangle
%   Column 2 is for edge from vertex 2 to 3 in the corresponding triangle
%   Column 3 is for edge from vertex 3 to 1 in the corresponding triangle
connList = TRI.ConnectivityList;    
sharpEdges = false(size(TRI));
sharpEdges(:,1) = ismember([connList(:,1), connList(:,2)], fe, 'rows');
sharpEdges(:,2) = ismember([connList(:,2), connList(:,3)], fe, 'rows');
sharpEdges(:,3) = ismember([connList(:,3), connList(:,1)], fe, 'rows');

end