function L = fogpl(R,f,T,den)
%fogpl    Path loss due to fog and cloud
%   L = fogpl(R,F,T,DEN) returns the path loss, L (in dB), due to fog or
%   cloud with a given liquid water density in the cloud or fog. R is an
%   length-M vector whose entries represent the propagation distances (in
%   meters). F is an length-N vector whose entries represent the signal
%   carrier frequency (in Hz). T is a scalar indicating the temperature (in
%   Celsius) and DEN is a scalar specifying the liquid water density (in
%   g/m^3). The liquid water density in fog is typically about 0.05 g/m^3
%   for medium fog (visibility of the order of 300m) and 0.5 g/m^3 for
%   thick fog (visibility of the order of 50 m).
%
%   L is an MxN matrix whose elements represents the path loss of each
%   propagation path under corresponding frequency.
%
%   The ITU fog model is valid between 10 to 1000 GHz.
%
%   % Examples:
%
%   % Example 1:
%   %   Compute the path loss of a 3 km link at 10 GHz due to a medium fog.
%   %   Assume the temperature is 15 degrees Celsius. 
%   
%   L = fogpl(3000,10e9,15,0.05)
%
%   % Example 2:
%   %   Plot the attenuation model for a thick fog from 10 GHz to 1000 GHz.
%
%   freq = (10:1000)*1e9;
%   L = fogpl(1e3,freq,15,0.5);
%   loglog(freq/1e9,L); grid on;
%   xlabel('Frequency (GHz)'); ylabel('Fog attenuation (dB/km)')
%
%   See also cranerainpl, fspl, gaspl, rainpl.

%   Copyright 2015-2016 The MathWorks, Inc.

% References
% [1] Recommendation ITU-R P.840-6, 2013
% [2] John Seybold, Introduction to RF Propagation, Wiley, 2005

%#codegen

validateattributes(R,{'double'}, {'nonnan','nonempty','real', ...
    'nonnegative','vector','finite'}, 'fogpl', 'R');
validateattributes(f,{'double'}, {'finite','nonempty','real', ...
    'vector','positive','>=',10e9,'<=',1000e9}, 'fogpl', 'F');
validateattributes(T,{'double'}, {'finite','nonempty','real', ...
    'scalar'}, 'fogpl', 'T');
validateattributes(den,{'double'}, {'finite','nonempty','real', ...
    'scalar','nonnegative'}, 'fogpl', 'DEN');

gamma = fogatt(f(:).',T,den);
Rkm = R(:)/1e3;
L = bsxfun(@times,Rkm,gamma);

end

function gamma = fogatt(f,Tc,den)

fGHz = f/1e9;
T = Tc + 273.15; % Convert Celsius to Kelvin 

theta = 300/T;
epsilon0 = 77.66+103.3*(theta-1);
epsilon1 = 0.0671*epsilon0;
epsilon2 = 3.52;

fpGHz = 20.20-146*(theta-1)+316*(theta-1).^2;
fsGHz = 39.8*fpGHz;

epsilon2p = fGHz.*(epsilon0-epsilon1)./(fpGHz.*(1+(fGHz./fpGHz).^2)) + ...
    fGHz.*(epsilon1-epsilon2)./(fsGHz.*(1+(fGHz./fsGHz).^2));
epsilon1p = (epsilon0-epsilon1)./(1+(fGHz./fpGHz).^2) + ...
    (epsilon1-epsilon2)./(1+(fGHz./fsGHz).^2) + epsilon2;

eta = (2+epsilon1p)./epsilon2p;
Kl = 0.819*fGHz./(epsilon2p.*(1+eta.*eta));

gamma = Kl*den;

end

% [EOF]