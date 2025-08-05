function [firstPt, firstFacetIdx, isNLOS, firstDist, isOnSharpEdge, ...
    isOnEdge] = firstIntersect(env, txPos, rxPosOrDir, rayType, varargin)
%FIRSTINTERSECT Calculate the first intersection point
%   FIRSTPT = COMM.INTERNAL.GEOMETRY.FIRSTINTERSECT(ENV, TXPOS, RXPOS,
%   'SEGMENT') calculates the first intersection point(s) from the
%   transmitter (Tx) position(s), TXPOS, to the receiver (Rx) position(s),
%   RXPOS, given an an environmental description, ENV. TXPOS/RXPOS is a
%   double precision, 1-by-3 vector or N-by-3 matrix with each row
%   representing one Tx/Rx position in the form of [x y z].
%
%   ENV is a structure with three fields: 
%     Triangulation:  Triangulation object for the environmental geometry
%     RayTracer:      Ray tracer object built from the triangulation object
%     SharpEdgeFlags: Logical matrix of the same size as
%                     Triangulation.ConnnectivityList indicating each edge
%                     of every triangle is sharp or not
%
%   The output FIRSTPT is a double precision, 1-by-3 vector when both TXPOS
%   and RXPOS are 1-by-3 vectors and N-by-3 matrix otherwise. Each row of
%   FIRSTPT is the first intersection point in the form of [x y z] between
%   the corresponding Tx and Rx positions. When Tx and Rx have
%   line-of-sight (LOS), the corresponding row in FIRSTPT is filled with
%   NaN. 
%   
%   FIRSTPT = COMM.INTERNAL.GEOMETRY.FIRSTINTERSECT(ENV, TXPOS, DIR, 'RAY')
%   calculates the first intersection point(s) from the Tx position(s),
%   TXPOS, along the direction(s), DIR. DIR is a double precision, 1-by-3
%   vector or N-by-3 matrix with each row representing one ray direction in
%   the form of [x y z].
%
%   [FIRSTPT, FIRSTFACETIDX, ISNLOS, FIRSTDIST, isOnSharpEdge, isOnEdge] =
%   COMM.INTERNAL.GEOMETRY.FIRSTINTERSECT(... ) returns the intersecting
%   facet index, FIRSTFACETIDX, logical ISNLOS to indicate whether or not
%   it is non-line-of-sight (NLOS) between Tx and Rx, distance from Tx to
%   the first intersection point, FIRSTDIST, logical ISONSHARPEDGE to
%   indicate whether or not the intersection point is on a sharp edge, and
%   logical ISONEDGE to indicate whether or not the intersection point is
%   on an edge of the intersecting facet. FIRSTFACETIDX, ISNLOS, FIRSTDIST,
%   and ISONSHARPEDGE are N-by-1 vectors. ISONEDGE is a N-by-3 matrix for
%   the edges from vertex 1/2/3 to vertex 2/3/1 on the corresponding
%   intersecting facets. 
% 
%   FIRSTFACETIDX is NaN, FIRSTDIST is inf, ISONSHARPEDGE is false, and
%   ISONEDGE is [false, false, false] when the corresponding Tx and Rx have
%   LOS.
%   
%   [...] = COMM.INTERNAL.GEOMETRY.FIRSTINTERSECT(...,
%   EXCLUDESHARPEDGEPOINT) changes the row of FIRSTPT, FIRSTFACETIDX and
%   FIRSTDIST to NaN if the first intersection point is on a sharp edge,
%   when EXCLUDESHARPEDGEPOINT is set to true. The corresponding Tx and Rx
%   are still considered NLOS, i.e., ISNLOS is independent of
%   EXCLUDESHARPEDGEPOINT. EXCLUDESHARPEDGEPOINT is false by default.

% Copyright 2019-2021 The MathWorks, Inc.

% Known limitation: The function may not work when Tx and Rx are on two
% different facets that belong to the same flat surface.

narginchk(4, 5);

excludeSharpEdgePt = (nargin == 5) && varargin{1};

% Perform ray tracing to find the first intersection point using Embree
RT = env.RayTracer;
sqrtTol = double(sqrt(eps('single')));
if strcmpi(rayType, 'segment')
    [dir, dist] = matlabshared.internal.segmentToRay(txPos, rxPosOrDir);
    [firstPt, firstFacetIdx, firstDist, firstUV] = firstIntersection( ...
        RT, txPos, dir, [repmat(sqrtTol, size(dist)) dist-sqrtTol]);
else % ray
    numRays = max(size(txPos,1), size(rxPosOrDir,1));
    [firstPt, firstFacetIdx, firstDist, firstUV] = firstIntersection( ...
        RT, txPos, rxPosOrDir, repmat([sqrtTol inf], numRays, 1));    
end

isNLOS = ~isnan(firstFacetIdx);
firstDist(~isNLOS) = inf; % distance is defined as inf for LOS rays

% Use barycentric coordinates to detect the first intersection point is on
% a sharp/knife edge or not:
%   * v = 0     if the point is on the edge from vertex 1 to 2
%   * 1-u-v = 0 if the point is on the edge from vertex 2 to 3
%   * u = 0     if the point is on the edge from vertex 3 to 1
baryC = [firstUV(:, 2) 1-sum(firstUV, 2) firstUV(:, 1)];
isOnEdge = (abs(baryC) < sqrtTol);     
    
% Assign a dummy facet number for LOS rays. The isOnEdge must be false for
% LOS rays because firstUV is [NaN, NaN].
tempFacetIdx = firstFacetIdx;
tempFacetIdx(~isNLOS) = 1; 
isOnSharpEdge = any(env.SharpEdgeFlags(tempFacetIdx,:) & isOnEdge, 2);

% Return NaN when the intersection point is on a sharp edge if
% EXCLUDESHARPEDGEPT is set to true.
if excludeSharpEdgePt && any(isOnSharpEdge)
    firstPt(isOnSharpEdge, :) = nan;
    firstFacetIdx(isOnSharpEdge) = nan;
    firstDist(isOnSharpEdge) = nan;
end

end

% [EOF]