function prjPts = project(TR, pt, varargin)
%PROJECTION Calculate the projection points on facets
%   PRJPTS = COMM.INTERNAL.GEOMETRY.PROJECT(TR, PT) calculates the
%   projections of a point PT on all the facets of a triangulation object
%   TR. PT is a double precision 1-by-3 vector in the form of [x, y, z].
%   PRJPTS is a double precision, N-by-3 matrix with each row representing
%   the projection on the corresponding facet in TR, where N is the number
%   of facets of TR, i.e., size(TR, 1).
%
%   PRJPTS = COMM.INTERNAL.GEOMETRY.PROJECT(TR, PT, FACETIDX) calculates
%   the projections of a point PT against select facets as specified by
%   FACETIDX. FACETIDX is an M-length column vector with integer values
%   between 1 and N, inclusive. Correspondingly, IMGPTS is a M-by-3 matrix.

% Copyright 2018-2019 The MathWorks, Inc.

% Math: If p is an arbitrary point in a plane whose normal is n, the
% projection of point x on that plane is y = x - <x-p, n> * n.

narginchk(2,3);

if nargin == 3
    facetIdx = varargin{1};
else    
    facetIdx = (1:size(TR, 1))';
end

% Get normal for the specified facets
fn = faceNormal(TR, facetIdx);

% Any one point from the facet will work. We can use a vertex of the facet
% but turns out the array indexing is very slow. 
% facetPt = TR.Points(TR.ConnectivityList(facetIdx,1), :);

% So, instead, we use the center point of the facet. 
facetPt = incenter(TR, facetIdx);

% Calculate projection points
prjPts = pt - (sum((pt - facetPt).* fn, 2) .* fn);

end

% [EOF]