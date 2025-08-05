function L = cranerainpl(R,f,rr,el_in,tau_in)
%cranerainpl    Path loss due to rain using the Crane model
%   L = cranerainpl(R,F,RR) returns the path loss, L (in dB), due to a
%   long-term statistical rain rate using the Crane attenuation model. R is
%   a length-M vector whose entries represent the propagation distances (in
%   meters). F is a length-N vector whose entries represent the signal
%   carrier frequency (in Hz). RR is a scalar indicating the long-term
%   statistical rainfall rate (in mm/hr). The Crane model is validated up
%   to 22.5 km.
%
%   L is an MxN matrix whose elements represent the path loss of each
%   propagation path under the corresponding frequency.
%
%   L = cranerainpl(R,F,RR,EL) specifies the elevation angle (in deg) of
%   the propagation path in EL. EL can be either a scalar or a vector. If
%   EL is a scalar, all propagation paths have the same elevation angle. If
%   EL is a vector, its length must be M so each element in EL is the
%   elevation angle of the corresponding propagation path in R. The default
%   value of EL is 0.
%
%   L = cranerainpl(R,F,RR,EL,TAU) specifies the polarization tilt angle
%   (in deg) of the signal in TAU. TAU can be either a scalar or a vector.
%   If TAU is a scalar, all propagation paths have the same tilt angle. If
%   TAU is a vector, its length must be M so each element in TAU is the
%   tilt angle of the signal of the corresponding propagation path. The
%   default value of TAU is 0, representing horizontal polarization.
%
%   The Crane global attenuation model was developed for Earth-space and
%   terrestrial paths and was based on geophysical observations of rain
%   rate, rain structure, and the vertical variation of temperature in the
%   atmosphere.
%
%   The k and alpha coefficients that compose the specific attenuation
%   factor gamma given by Crane and ITU-R P.838 (International
%   Telecommunication Union) are identical and are valid from 1 to 1000
%   GHz.
%
%   The 0.01% rainfall rate tables provided by Crane and the ITU are
%   different. While both models are international, the Crane rain regions
%   are more specific to North America. Comparing results, the ITU rain
%   attenuation model generally computes lower rain attenuation values than
%   Crane.
%
%   % Example 1:
%   %   Compute the path loss of a 3 km link at 1 GHz due to a rain rate
%   %   of 95 mm/h.
%   Lcrane = cranerainpl(3000,1e9,95)
% 
%   % Example 2:
%   %   Plot the attenuation model for a rain rate of 95 mm/h from 1 GHz to
%   %   1000 GHz.
%   freq = (1:1000)*1e9; % Hz
%   Lcrane = cranerainpl(1e3,freq,95);
%   loglog(freq/1e9,Lcrane)
%   grid on
%   xlabel('Frequency (GHz)')
%   ylabel('Rain attenuation (dB/km)')
%
%   See also fogpl, fspl, gaspl, rainpl.

%   Copyright 2015-2023 The MathWorks, Inc.

% References
%   [1] Crane, R. K. Electromagnetic Wave Propagation through Rain, New
%       York: John Wiley, 1996.
%   [2] International Telecommunication Union (ITU). "Specific Attenuation
%       Model for Rain for Use in Prediction Methods." Recommendation ITU-R
%       P.838-3, P Series, Radiowave Propagation, 2005.

%#codegen
narginchk(3,5);

validateattributes(R, {'double'}, {'nonnan','nonempty','real', ...
    'nonnegative','vector','finite','<=',22.5e3}, 'cranerainpl', 'R');
validateattributes(f, {'double'}, {'finite','nonempty','real', ...
    'vector','positive','>=',1e9,'<=',1000e9}, 'cranerainpl', 'F');
validateattributes(rr, {'double'}, {'finite','nonempty','real'...
    'nonnegative','scalar'}, 'cranerainpl', 'RR');

M = numel(R);
N = numel(f); 

if nargin >= 4 
    if isscalar(el_in)
        el = repmat(el_in,1,M);
    else
        el = el_in;
    end
    validateattributes(el, {'double'}, {'finite','nonempty', ...
        'real','vector','numel',M,'>=',-90,'<=',90}, 'cranerainpl', 'EL');
else
    el = zeros(1, M);
end

if nargin == 5
    if isscalar(tau_in)
        tau = repmat(tau_in,1,M);
    else
        tau = tau_in;
    end    
    validateattributes(tau, {'double'}, {'finite','nonempty', ...
        'real','vector','numel',M,'>=',-90,'<=',90}, 'cranerainpl', 'TAU');
else
    tau = zeros(1, M);
end

[gamma,alpha] = phased.internal.rainatt(f(:).',rr,el(:),tau(:));
Rkm = R(:)/1e3;

L = zeros(M,N);
if rr > 0
    % Parameters
    b = 2.3*rr.^-0.17;
    c = 0.026 - 0.03.*log(rr);
    delta = 3.8 - 0.6.*log(rr);
    mu = log(b*exp(c*delta))/delta;
    y = mu*alpha;
    z = c*alpha;
    
    % First region
    idx = Rkm<=delta;
    if any(idx)
        L(idx,:) = gamma(idx,:).*(exp(bsxfun(@times,y(idx,:),Rkm(idx))) - 1)./y(idx,:); % bsxfun to handle implicit expansion
    end
    
    % Second region
    idx = ~idx;
    if any(idx)
        p1 = (exp(y(idx,:).*delta) - 1)./y(idx,:);
        p2 = (b.^alpha(idx,:)).*exp(z(idx,:).*delta)./z(idx,:);
        p3 = (b.^alpha(idx,:)).*exp(bsxfun(@times,z(idx,:),Rkm(idx)))./z(idx,:); % bsxfun to handle implicit expansion
        L(idx,:) = gamma(idx,:).*(p1 - p2 + p3);
    end
end

end

% LocalWords:  ITU Lcrane fogpl fspl gaspl rainpl Radiowave
