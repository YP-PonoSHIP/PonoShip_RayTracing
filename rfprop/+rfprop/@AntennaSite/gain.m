function G = gain(site, fq, azs, els, varargin)
%gain   Antenna gain
%   G = gain(SITE, FQ) returns the peak gain in dBi of SITE's antenna.
%
%   G = gain(SITE, FQ, AZ, EL) returns the gain in dBi of SITE's antenna
%   at the specified frequency FQ and angles AZ and EL.
%
%   G = gain(SITE, FQ, AZ, EL, GPATTERN, GAZ, GEL) returns the gain of
%   SITE's antenna at the specified frequency FQ and angles AZ and EL given
%   antenna gain pattern GPATTERN and corresponding angles GAZ and GEL.
%
%   The inputs AZ and EL can be scalars or vectors. If vectors, the gain G
%   is a vector of the same length corresponding to the gain of each AZ,EL
%   pair.

%   Copyright 2017-2021 The MathWorks, Inc.

% Implement isotropic antenna
ant = site.Antenna;
if ischar(ant)
    if nargin < 3
        G = 0; % Peak gain
    else
        G = zeros(numel(azs),1);
    end
    return
end

% If no angles specified, compute peak gain in any direction
if nargin < 3    
    Gpattern = pattern(ant, fq);
    G = max(Gpattern(:));
    return
end

% Transform angle to local coordinates defined by site's AntennaAngle
[azs, els] = rfprop.internal.transformAngle(site.AntennaAngle, azs, els);

if nargin == 5
    % Use Phased Array directivity object
    phasedDirectivity = varargin{1};
    G = step(phasedDirectivity, fq, [azs';els']);
elseif nargin > 4 || ~rfprop.AntennaSite.isCommAntenna(ant)
    % Get antenna gain pattern
    if nargin == 4
        [Gpattern,Gaz,Gel] = gainPattern(site,fq);
    else
        Gpattern = varargin{1};
        Gaz = varargin{2};
        Gel = varargin{3};
    end

    % Query antenna pattern for input angles
    if sum(any(~isfinite(Gpattern))) > 0
        interpMethod = "linear";
    else
        interpMethod = "spline";
    end
    G = interp2(Gaz, Gel, Gpattern, azs, els, interpMethod);
else
    % Compute gain for each angle from pattern. This approach is used for
    % arrayConfig, which models a basic phased array. Computing gain for
    % each angle produces equivalent results as for phased arrays with
    % directivity patterns, e.g., phased.ULA objects.
    numAngles = numel(azs);
    G = zeros(numAngles,1);
    for k = 1:numAngles
        G(k) = pattern(ant, fq, azs(k), els(k));
    end
end