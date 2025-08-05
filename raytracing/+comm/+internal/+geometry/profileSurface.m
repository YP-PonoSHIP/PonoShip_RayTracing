function surfacePos = profileSurface(TRI, queryPos, varargin)
%PROFILESURFACE Find the surface (outermost) point
%   SURFACEPOS = COMM.INTERNAL.GEOMETRY.PROFILESURFACE(TRI, QUERYPOS)
%   calculates the geometric surface position(s), SURFACEPOS, vertically
%   from the query position(s), QUERYPOS, given a triangulation or
%   matlabshared.internal.StaticSceneRayTracer object, TRI. SURFACEPOS can
%   be above or below QUERYPOS in z-axis.
%   
%   QUERYPOS is a double precision, N-by-3 matrix with each row
%   representing one query position in the form of [x, y, z]. SURFACEPOS is
%   of the same data type and dimension as QUERYPOS. When there is no
%   surface point, the corresponding row in QUERYPOS is filled with NaN.
%   
%   SURFACEPOS = COMM.INTERNAL.GEOMETRY.PROFILESURFACE(TRI, QUERYPOS,
%   QUERYDIR) calculates the surface position(s) from the query position(s)
%   on the direction(s), QUERYDIR. QUERYDIR is a double precision, 1-by-3
%   vector applying to all query positions or N-by-3 matrix with each row
%   applying to one query position. Each row of QUERYDIR specifies the
%   direction in the form of [x, y, z] where z must be non-zero.

% Copyright 2019-2021 The MathWorks, Inc.

% Algorithm: For each query position, find a point on the positive/negative
% query direction that is above max Z-axis in the geometry (so definitely
% not on a facet). Then shot a ray from that point downward (the inverse of
% query direction). The first intersection point is the function return. 

narginchk(2, 3);

if isa(TRI, 'triangulation') % Build AABB Tree
    RT = matlabshared.internal.StaticSceneRayTracer(TRI);
    % Ray origin is 10 meters above the max Z value
    rayOriginZ = max(TRI.Points(:, 3)) + 10; 
else % matlabshared.internal.StaticSceneRayTracer object
    RT = TRI;
    % Ray origin at 100,000 meters above the ground should be good for RF
    % Prop applications. 
    rayOriginZ = 1e5; 
end

% Get query direction
if ~isempty(varargin)
    queryDir = varargin{1};
else
    queryDir = [0 0 1];
end

% Get the origin for each ray to be shot
t = (rayOriginZ - queryPos(:,3))./queryDir(:, 3); 
rayOrigin = queryPos + t .* queryDir; 

% Perform ray-tracing to find the surface positions (first intersection
% point for each ray)
dir = matlabshared.internal.segmentToRay(rayOrigin, queryPos);
surfacePos = firstIntersection(RT, rayOrigin, dir);

end

% [EOF]