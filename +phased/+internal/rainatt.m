function [gamma,alpha,k] = rainatt(f,rr,el,tau)
%   This class is for internal use only. It may be removed in the future.

% rainatt    Calculate specific attenuation and coefficients
%   GAMMA = rainatt(F,RR,EL,TAU) calculates the specific attenuation GAMMA
%   as defined by ITU-R P.838-3. F is an length-N vector whose entries
%   represent the signal carrier frequency (in Hz). RR is a scalar
%   indicating the rain rate (in mm/h). EL specifies the elevation angle
%   (in degrees) of the propagation path. EL can be either a scalar or a
%   vector. If EL is a scalar, all propagation paths have the same
%   elevation angle. If EL is a vector, its length must be M so each
%   element in EL is the elevation angle of the corresponding propagation
%   path in R. TAU specifies the polarization tilt angle (in degrees). TAU
%   can be either a scalar or a vector. If TAU is a scalar, all propagation
%   paths have the same tilt angle. If TAU is a vector, its length must be
%   M so each element in TAU is the tilt angle of the signal of the
%   corresponding propagation path.
%
%   [GAMMA,ALPHA,K] = rainatt(F,RR,EL,TAU) outputs the attenuation
%   coefficients ALPHA and K.
%
%   Copyright 2019 The MathWorks, Inc.

%   Reference 
%   [1] International Telecommunication Union (ITU). "Specific Attenuation
%       Model for Rain for Use in Prediction Methods." Recommendation ITU-R
%       P.838-3, P Series, Radiowave Propagation, Mar. 2005.

%#codegen

[kH,kV,alphaH,alphaV] = rainattcoeff(f);

k = bsxfun(@plus,kH+kV,bsxfun(@times,(kH-kV),cosd(el).^2.*cosd(2*tau)))/2;
alpha = bsxfun(@plus,kH.*alphaH+kV.*alphaV,...
    bsxfun(@times,(kH.*alphaH-kV.*alphaV),cosd(el).^2.*cosd(2*tau)))./(2*k);

gamma = k.*rr.^alpha;
end

function [kH,kV,alphaH,alphaV] = rainattcoeff(f)

fGHz = f/1e9;  % convert to GHz

kHtab = [-5.33980 -0.10008 1.13098; ...
         -0.35351 1.26970 0.45400; ...
         -0.23789 0.86036 0.15354; ...
         -0.94158 0.64552 0.16817];
kHm = -0.18961;
kHc = 0.71147;

kVtab = [-3.80595 0.56934 0.81061; ...
         -3.44965 -0.22911 0.51059; ...
         -0.39902 0.73042 0.11899; ...
         0.50167 1.07319 0.27195];
kVm = -0.16398;
kVc = 0.63297;

alphaHtab = [-0.14318 1.82442 -0.55187; ...
             0.29591 0.77564 0.19822; ...
             0.32177 0.63773 0.13164; ...
             -5.37610 -0.96230 1.47828; ...
             16.1721 -3.29980 3.43990];
alphaHm = 0.67849;
alphaHc = -1.95537;

alphaVtab = [-0.07771 2.33840 -0.76284; ...
             0.56727 0.95545 0.54039; ...
             -0.20238 1.14520 0.26809; ...
             -48.2991 0.791669 0.116226; ...
             48.5833 0.791459 0.116479];
alphaVm = -0.053739;
alphaVc = 0.83433;

tempkH = bsxfun(@rdivide,bsxfun(@minus,log10(fGHz),kHtab(:,2)),kHtab(:,3));
log10kH = sum(bsxfun(@times,kHtab(:,1),exp(-tempkH.^2)))+...
    kHm.*log10(fGHz)+kHc;
kH = 10.^log10kH;

tempkV = bsxfun(@rdivide,bsxfun(@minus,log10(fGHz),kVtab(:,2)),kVtab(:,3));
log10kV = sum(bsxfun(@times,kVtab(:,1),exp(-tempkV.^2)))+...
    kVm.*log10(fGHz)+kVc;
kV = 10.^log10kV;

tempalphaH = bsxfun(@rdivide,bsxfun(@minus,log10(fGHz),alphaHtab(:,2)),alphaHtab(:,3));
alphaH = sum(bsxfun(@times,alphaHtab(:,1),exp(-tempalphaH.^2)))+...
    alphaHm.*log10(fGHz)+alphaHc;

tempalphaV = bsxfun(@rdivide,bsxfun(@minus,log10(fGHz),alphaVtab(:,2)),alphaVtab(:,3));
alphaV = sum(bsxfun(@times,alphaVtab(:,1),exp(-tempalphaV.^2)))+...
    alphaVm.*log10(fGHz)+alphaVc;

end

