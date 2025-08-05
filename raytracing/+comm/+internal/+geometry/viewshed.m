function visFacetIdx = viewshed(TRI, RT, pt, varargin)
%VIEWSHED Calculate the viewshed of query point(s)
%   VISFACETIDX = COMM.INTERNAL.GEOMETRY.VIEWSHED(TRI, TR, PT) calculates
%   the indices of the facets that are visible from the query point(s), PT,
%   given a triangulation object, TRI, and the corresponding
%   comm.internal.RayTracer object, RT. PT is a double precision N-by-3
%   matrix with each row representing one point position in the form of [x,
%   y, z]. VISFACETIDX is a N-by-1 cell array. Each cell element is an
%   integer column vector containing the visible facet indices for the
%   corresponding element in PT.
%
%   VISFACETIDX = COMM.INTERNAL.GEOMETRY.VIEWSHED(TRI, TR, PT, ONFACETIDX)
%   specifies the facet indices, ONFACETIDX, on which PT are known to be.
%   This option saves the extra step to check whether or not PT is on a
%   facet. ONFACETIDX is an integer, N-by-1 vector.
%
%   This function is to work around the limitation of the
%   comm.internal.RayTracer.viewshed method when a query point is on a
%   facet.

% Copyright 2019-2020 The MathWorks, Inc.

narginchk(3, 4);

% Tolerance distance to claim two [x y z] points are identical
square_tol = eps;
tol = sqrt(square_tol);

% Distance to push away from facet when a point is on it
pushOutDist = 1e-5; 

% Get the number of query points
numPts = size(pt, 1);

if nargin == 3
    % Check if a query point is on a facet
    [~, onFacetIdx, dist] = closestPoint(RT, pt);
    % Get indices of query points that are on a facet
    isPtOnFacetIdx = (dist < tol);
else % It is known that all the query points are on a facet
    onFacetIdx = varargin{1};
    isPtOnFacetIdx = 1:numPts;
end

% For the points on a facet, push them a little away from the facet so
% the viewshed method won't return empty.
queryPt = pt;
if any(isPtOnFacetIdx)
    queryPt(isPtOnFacetIdx, :) = pt(isPtOnFacetIdx, :) + ...
        pushOutDist * faceNormal(TRI, onFacetIdx(isPtOnFacetIdx));
end

% Iterate through all the points
visFacetIdx = cell(numPts, 1);
for i = 1:numPts
    % The uniform ray-tracing method needs to shoot 101,700 rays with a
    % angular separation of .8 degree. For a triangulation object with <
    % 15,000 facets, we shoot a max of 60,000 rays to determine visibility
    % using the facet check method. 
    if size(TRI, 1) < 15e3
        visFacetIdx{i} = getViewshedAllFacets( ...
            TRI, RT, queryPt(i, :));    
    else        
        visFacetIdx{i} = getViewshedUniformRT( ...
            RT, queryPt(i, :), .8);
    end
     
    if isPtOnFacetIdx(i) 
        % Exclude all the facets that contain or extend to the point. In
        % other words, exclude the facets for which the point and its
        % projection/image are identical.
        prjPts = comm.internal.geometry.project(TRI, ...
            pt(i,:), visFacetIdx{i});
        visFacetIdx{i} = visFacetIdx{i}( ...
            sum((pt(i, :) - prjPts).^2,2) > square_tol);
    end
end
 
end

function facetIdx = getViewshedAllFacets(TRI, RT, pt)
% Calculate viewshed by checking the visibility of each facet. A facet is
% visible when
%   1) Its normal point toward the query point
%   2) Any of the 3 vertices or the centroid is visible to the query point

fn = faceNormal(TRI);   % Facet normal
fc = incenter(TRI);     % Facet center

% Get facets whose normals point to the query point and whose centroids are
% visible to the query point
tol = double(sqrt(eps('single')));
[dir, dist] = matlabshared.internal.segmentToRay(fc, pt);
isHit = hasIntersection(RT, fc, dir, [repmat(tol, size(dist)) dist]);
facetIdx = find(dist <= tol | ~isHit);

% Find all the facets whose normals point toward the query point
visFacetIdxNormal = find(sum(fn .* (pt-fc), 2) > 0);

% Only check the facets that are excluded by the centroid check
visFacetIdxNormal = setdiff(visFacetIdxNormal, facetIdx);
if isempty(visFacetIdxNormal)
    return;
end

% Get all the unique vertices on the candidate facets. So we don't have to
% shoot a ray to a vertex more than once when the vertex belongs to more
% than one facet.
connList  = TRI.ConnectivityList(visFacetIdxNormal, :);
vertexIdx = unique(connList(:));
vertex = TRI.Points(vertexIdx, :);

% Shoot a ray from the query point to every vertex.
[dir, dist] = matlabshared.internal.segmentToRay(pt, vertex);
isHit = hasIntersection(RT, pt, dir, [repmat(tol, size(dist)) dist-tol]);

% If the intersection returns empty, the vertex is visible to the query
% point. We label it.
isVertexVisible = false(1, size(TRI.Points, 1));
isVertexVisible(vertexIdx) = ~isHit;

% Get all the candidate facets that have at least one vertex visible but
% centroid invisible
visFacetIdxVertex = visFacetIdxNormal(any(isVertexVisible(connList), 2));

% Aggregate two visible facet candidates
facetIdx = [facetIdx; visFacetIdxVertex];

end

function facetIdx = getViewshedUniformRT(RT, pt, angSep)
% To get visible facets (viewshed), we shoot many rays at uniformed
% distributed az & el angles and gather all the facets being first
% intersected. ANGSEP is the uniform angle separation for az and el. 

% Uniform az and el angles in radians
az = (-180+angSep:angSep:180)/180*pi;
el = (-90:angSep:90)/180*pi;

% Get all az & el combinations
[az, el] = meshgrid(az, el);

% Calculate directions for all rays to be shot
[x, y, z] = sph2cart(az(:), el(:), 1);

% Shoot rays
tol = double(sqrt(eps('single')));
[~, intFacets] = firstIntersection(RT, pt, [x y z], ...
    repmat([tol Inf], length(x), 1));

% Gather facets being first intersected by any one ray
facetIdx = unique(intFacets(~isnan(intFacets)));

end

% [EOF]