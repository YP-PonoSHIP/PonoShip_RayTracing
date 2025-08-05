function [X, Y, Z]= sph2cart(phi, theta, r)
%

%   Copyright 2015-2020 The MathWorks, Inc.

Z  = r.*cosd(theta);
X  = r.*sind(theta).*cosd(phi);
Y  = r.*sind(theta).*sind(phi);