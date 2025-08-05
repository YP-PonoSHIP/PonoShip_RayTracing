function imgPts = mirror(TR, point, varargin)
%MIRROR Calculate the imaging points against facets
%   IMGPTS = COMM.INTERNAL.GEOMETRY.MIRROR(TR, PT) calculates the images of
%   a point PT against all the facets of a triangulation object TR. PT is a
%   double precision 1-by-3 vector in the form of [x, y, z]. IMGPTS is a
%   double precision, N-by-3 matrix with each row representing the image
%   against the corresponding facet in TR, where N is the number of facets
%   of TR, i.e., size(TR, 1).
%
%   IMGPTS = COMM.INTERNAL.GEOMETRY.MIRROR(TR, PT, FACETIDX) calculates the
%   images of a point PT against select facets as specified by FACETIDX.
%   FACETIDX is an M-length column vector with integer values between 1 and
%   N, inclusive. Correspondingly, IMGPTS is a M-by-3 matrix.

% Copyright 2018-2019 The MathWorks, Inc.

if nargin == 3
    facetIdx = varargin{1};
    imgPts = 2*comm.internal.geometry.project(TR, point, facetIdx) - point;    
else    
    imgPts = 2*comm.internal.geometry.project(TR, point) - point;    
end

end

% [EOF]