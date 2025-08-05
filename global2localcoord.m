function lclCoord = global2localcoord(gCoord,optionArg,localOriginArg,localAxesArg)
%global2localcoord Global to local coordinates conversion
%   LCOORD = global2localcoord(GCOORD, OPTION) converts the global
%   coordinate GCOORD to its local counterpart LCOORD.  The global
%   coordinates' origin is located at [0; 0; 0] and its three axes are in
%   the directions of [1; 0; 0], [0; 1; 0] and [0; 0; 1], respectively.
%   OPTION specifies the form of GCOORD and LCOORD using one of 'rr' | 'rs'
%   | 'ss' | 'sr', where the default is 'rr'. The 'rr' option converts
%   global rectangular coordinates to local rectangular coordinates. The
%   'rs' option converts global rectangular coordinates to local spherical
%   coordinates. The 'sr' option converts global spherical coordinates to
%   local rectangular coordinates. The 'ss' option converts global
%   spherical coordinates to local spherical coordinates. GCOORD and LCOORD
%   are 3xN matrices where each column is coordinates in either rectangular
%   or spherical coordinate form.
%
%   If the coordinates are in rectangular form, the three elements
%   represent (x,y,z) in meters. If the coordinates are in spherical form,
%   the three elements represent (az, el, r), where azimuth (az, in
%   degrees) is measured from x axis toward y axis, elevation (el, in
%   degrees) is measured from x-y plane toward z axis and r (in meters) is
%   the radius. If GCOORD is a matrix, the resulting LCOORD is also a
%   matrix with the same dimensions. Each column of LCOORD contains the
%   three coordinates that define the converted point.
%
%   LCOORD = global2localcoord(...,LORIGIN) specifies the origin, LORIGIN,
%   of the local coordinate system. LORIGIN is a 3xN matrix whose columns
%   contain the rectangular coordinates of the local coordinate system
%   origins with respect to the global coordinate system. If LORIGIN is a
%   column vector, the origin is shared by all N local coordinate systems.
%   The default for LORIGIN is [0; 0; 0].
%
%   LCOORD = global2localcoord(...,LAXES) specifies the axes, in LAXES, of
%   the local coordinate system. LAXES is a 3x3xN array whose pages are
%   axes of the local coordinate systems. Within each page, each column
%   specifies the direction of local x, y and z axis in the global
%   coordinate system. If LAXES is a matrix, the local coordinate system is
%   the same for all conversions. The default for LAXES is [1 0 0;0 1 0;0 0
%   1]. In general, LAXES is orthonormal and you can use isorthonorm()
%   function to check if the specified LAXES is orthonormal.
%
%   % Example:
%   %   A point is located at (2,3,0) in the global coordinate system.
%   %   Determine the coordinates of this point in the local coordinate 
%   %   system whose origin is at (1,1,0) and whose axes are [0; 1; 0], 
%   %   [1; 0; 0] and [0; 0; -1].
%
%   lclcoord = global2localcoord([2; 3; 0],'rr',[1; 1; 0],...
%               [0 1 0;1 0 0;0 0 -1])
%
%   See also phased, local2globalcoord, rangeangle.

%   Copyright 2008-2019 The MathWorks, Inc.

%   Reference
%   [1] J. D. Foley et al. Computer Graphics: Principles and Practice in C,
%       2nd Ed., Addison-Wesley, 1995

%#codegen 
%#ok<*EMCA>

narginchk(1,4);

if nargin < 4 || isempty(localAxesArg)
    localAxes = eye(3);
else
    localAxes = localAxesArg;
end
if nargin < 3 || isempty(localOriginArg)
    localOrigin = [0; 0; 0];
else
    localOrigin = localOriginArg;
end
if nargin < 2 || isempty(optionArg)
    option = 'rr';
else
    option = optionArg;
end

coder.internal.assertFixedSize( ...
    2:nargin, gCoord,option,localOrigin,localAxes);

option = validatestring(option,{'rs','rr','sr','ss'},...
    'global2localcoord','OPTION');
validateattributes(gCoord, {'double'}, ...
    {'finite','nonnan','nonempty','real','2d','nrows',3}, ...
    'global2localcoord','GCOORD');
validateattributes(localOrigin, {'double'}, ...
    {'finite','nonnan','nonempty','real','2d','nrows',3}, ...
    'global2localcoord','LORIGIN');
validateattributes(localAxes,{'double'},{'finite','nonnan','nonempty','real',...
    '3d','nrows',3,'ncols',3},'global2localcoord','LAXES');

numcoord = size(gCoord,2);
numorigin = size(localOrigin,2);
numaxes = size(localAxes,3);
cond = (numcoord~=numorigin) && (numcoord>1) && (numorigin>1);
if cond
    coder.internal.errorIf(cond, ...
        'shared_channel:shared_channel:invalidColumnNumbers', ...
        'LORIGIN',numcoord);
end
cond = (numaxes~=numorigin) && (numaxes>1) && (numorigin>1);
if cond
    coder.internal.errorIf(cond, ...
        'shared_channel:shared_channel:invalidPageNumbers', ...
        'LAXES',numorigin);
end

localAxes = bsxfun(@rdivide,localAxes,sqrt(sum(localAxes.^2)));

lclCoord = phased.internal.global2localcoord(gCoord,option,...
    localOrigin,localAxes);

% [EOF]