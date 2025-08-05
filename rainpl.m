function L = rainpl(R,f,rr,el_in,tau_in,pct)
%rainpl    Path loss due to rain using the ITU model
%   L = rainpl(R,F,RR) returns the path loss, L (in dB), due to a long-term
%   statistical rain rate. R is a length-M vector whose entries represent
%   the propagation distances (in meters). F is a length-N vector whose
%   entries represent the signal carrier frequency (in Hz). RR is a scalar
%   indicating the long-term statistical 0.01% rain rate (in mm/h), which
%   is the rain rate that is exceeded 0.01% of the time.
%
%   L is an MxN matrix whose elements represents the path loss of each
%   propagation path under the corresponding frequency.
%
%   L = rainpl(R,F,RR,EL) specifies the elevation angle (in degrees) of the
%   propagation path in EL. EL can be either a scalar or a vector. If EL is
%   a scalar, all propagation paths have the same elevation angle. If EL is
%   a vector, its length must be M so each element in EL is the elevation
%   angle of the corresponding propagation path in R. The default value of
%   EL is 0.
%
%   L = rainpl(R,F,RR,EL,TAU) specifies the polarization tilt angle (in
%   degrees) of the signal in TAU. TAU can be either a scalar or a vector.
%   If TAU is a scalar, all propagation paths have the same tilt angle. If
%   TAU is a vector, its length must be M so each element in TAU is the
%   tilt angle of the signal of the corresponding propagation path. The
%   default value of TAU is 0, representing horizontal polarization.
%
%   L = rainpl(R,F,RR,EL,TAU,PCT) returns the attenuation L for the
%   specified percentage of time PCT (%). PCT is a scalar in the range of
%   0.001 to 1, inclusive. The loss L in this case is deduced from a power
%   law using the long-term statistical 0.01% rain rate (in mm/h).
%
%   The ITU rain model is valid between 1 to 1000 GHz.
%
%   % Examples:
%
%   % Example 1:
%   %   Compute the path loss of a 3 km link at 1 GHz due to a rain rate of
%   %   95 mm/h. 
%   
%   L = rainpl(3000,1e9,95)
%
%   % Example 2:
%   %   Plot the attenuation model for a rain rate of 95 mm/h from 1 GHz to
%   %   1000 GHz over a range equal to 5 km.
%
%   freq = (1:1000)*1e9;
%   L = rainpl(5e3,freq,95);
%   loglog(freq/1e9,L)
%   grid on
%   xlabel('Frequency (GHz)')
%   ylabel('Rain attenuation (dB)')
%
%   See also cranerainpl, fogpl, fspl, gaspl.

%   Copyright 2015-2019 The MathWorks, Inc.

% References
% [1] Recommendation ITU-R P.530-17, 2017
% [2] Recommendation ITU-R P.838-3, 2005

%#codegen

narginchk(3,6);

validateattributes(R, {'double'}, {'nonnan','nonempty','real', ...
    'nonnegative','vector','finite'}, 'rainpl', 'R');
validateattributes(f, {'double'}, {'finite','nonempty','real', ...
    'vector','positive','>=',1e9,'<=',1000e9}, 'rainpl', 'F');
validateattributes(rr, {'double'}, {'finite','nonnan','nonempty','real'...
    'nonnegative','scalar'}, 'rainpl', 'RR');

M = numel(R);

if nargin >= 4 
    if isscalar(el_in)
        el = repmat(el_in,1,M);
    else
        el = el_in;
    end
    validateattributes(el, {'double'}, {'finite','nonnan','nonempty', ...
        'real','vector','numel',M,'>=',-90,'<=',90}, 'rainpl', 'EL');
else
    el = zeros(1, M);
end

if nargin >= 5
    if isscalar(tau_in)
        tau = repmat(tau_in,1,M);
    else
        tau = tau_in;
    end    
    validateattributes(tau, {'double'}, {'finite','nonnan','nonempty', ...
        'real','vector','numel',M,'>=',-90,'<=',90}, 'rainpl', 'TAU');
else
    tau = zeros(1, M);
end

if nargin == 6
    validateattributes(pct, {'double'}, {'finite','nonempty','real', ...
    'scalar','positive','>=',0.001,'<=',1}, 'rainpl', 'PCT');
end

% Calculate specific attenuation gamma
[gamma,alpha] = phased.internal.rainatt(f(:).',rr,el(:),tau(:));

% Calculate distance factor r
d = R(:)/1e3; % Actual path length (km)
fGHz = (f(:)./1e9).'; 

% Calculate denominator of effective range equation
% Use bsxfun to handle implicit expansion
% denom = 0.477.*d.^0.633.*rr.^(0.073.*alpha).*fGHz.^0.123 - 10.579*(1 - exp(-0.024.*d));
denom = bsxfun(@times,0.477.*d.^0.633,rr.^(0.073.*alpha)); 
denom = bsxfun(@times,denom,fGHz.^0.123);
denom = bsxfun(@minus,denom,10.579*(1 - exp(-0.024.*d)));

% Calculate effective range
r = 1./denom;
r(denom < 0.4) = 2.5; 

% Calculate effective path length (km)
deff = bsxfun(@times,r,d); % bsxfun to handle implicit expansion

L = bsxfun(@times,deff,gamma);

% Calculate loss for other percentages in the range of 0.001 - 1
if nargin == 6 && pct~=0.01
    C0 = zeros(1,numel(fGHz));
    C0(fGHz>=10) = 0.12 + 0.4*log10((fGHz(fGHz>=10)/10).^0.8);
    C0(fGHz<10) = 0.12; 
    C1 = (0.07.^C0).*(0.12.^(1 - C0));
    C2 = 0.855*C0 + 0.546*(1 - C0); 
    C3 = 0.139*C0 + 0.043*(1 - C0); 
    PL = C1.*pct.^(-(C2 + C3.*log10(pct))); 
    L = bsxfun(@times,L,PL); 
end

end

% LocalWords:  ITU cranerainpl fogpl fspl gaspl denom rr
