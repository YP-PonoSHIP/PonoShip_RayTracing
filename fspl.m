function L = fspl(R,lambda)
%fspl     Free space path loss
%   L = fspl(R,LAMBDA) returns the free space path loss L (in dB) suffered
%   by a signal with wavelength LAMBDA (in meters) when it is propagated in
%   free space for a distance of R (in meters). R can be a length-M vector
%   and LAMBDA is a length-N vector. L is an MxN matrix whose elements are
%   the free space path loss for the corresponding propagation distance
%   specified in R at the wavelength specified in LAMBDA. When LAMBDA is a
%   scalar, L has the same dimensions as R.
%
%   Note that the best case is lossless so the loss is always greater than
%   or equal to 0 dB.
%
%   % Examples:
%
%   % Example 1:
%   %   Calculate the free space loss for a signal whose wavelength is 30
%   %   cm. The signal is propagated for 1 km.
%
%   L = fspl(1000,0.3)
%
%   % Example 2:
%   %   Plot the attenuation for a 1 km free space propagation from 1 GHz
%   %   to 1000 GHz.
%
%   freq = (1:1000)*1e9;
%   c = physconst('lightspeed');
%   L = fspl(1e3,c./freq);
%   loglog(freq/1e9,L); grid on;
%   xlabel('Frequency (GHz)'); ylabel('Free space attenuation (dB)')
%
%   See also cranerainpl, fogpl, gaspl, rainpl.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] John Proakis, Digital Communications, 4th Ed., McGraw-Hill, 2001
%   [2] Recommendation ITU-R P.525-2, 1994

%#codegen 

validateattributes(R, {'double'}, {'nonnan','nonempty','real', ...
    'vector','nonnegative'}, 'fspl', 'R');
validateattributes(lambda, {'double'}, {'nonnan','nonempty','real', ...
    'vector','positive'}, 'fspl', 'LAMBDA');

if isscalar(lambda)
    L = 4*pi*R/lambda;
else
    L = 4*pi*R(:)*(1./(lambda(:).'));
end

L = max(L, 1);
L = 20*log10(L); % Convert to dB

end

% [EOF]