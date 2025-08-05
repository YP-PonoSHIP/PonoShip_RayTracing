function Edbu = pow2efield(Prx_db, fq, Grx_dbi)
%pow2efield   Convert power quantity to electric field quantity
%   E = rfprop.pow2efield(P,FQ,GRX) returns electric field quantity in
%   dBuV/m given power quantity P in dBm, transmitter frequency FQ, and
%   receiver gain GRX in dBi.

%   Copyright 2017-2019 The MathWorks, Inc.   

% Convert dB quantities to numeric quantities
Prx_W = 10^((Prx_db-30)/10);
Grx = 10^(Grx_dbi/10);

% Define constants
Z0 = rfprop.Constants.Z0;
lambda = rfprop.Constants.LightSpeed/fq;

% Compute signal strength in dBuV/m, using field strength equation derived 
% from following equations:
%  1. Prx = S * Ae,  (antenna efficiency equation)
%  2. Grx = 4 * pi * Ae /(lambda^2) (antenna gain equation)
%  3. S = E^2/Z0, Z0 = 120*pi (power flux density equation for plane wave since far-field)
%  4. Prx = Ptx * Gtx * Grx * Gpl (Friis equation in numeric form)
%
% Using substitution and solving for E yields:
%
%  Erx = sqrt(4 * pi * Z0 * Ptx * Gtx * Gpl) / lambda
%
%   or substituting Ptx * Gtx * Gpl = Prx / Grx:
%
%  Erx = sqrt(4 * pi * Z0 * Prx / Grx) / lambda
%
% where:
%  * Erx is magnitude of electric field incident on receiver antenna
%  * Ptx is transmitter output power (W)
%  * Gtx is numeric transmitter gain (including antenna gain and system loss)
%  * Gpl is numeric path loss gain
%
% Note that substituting Gpl = (lambda/4*pi*d)^2 for free space
% yields the more familiar equation E = sqrt(30*Ptx*Gtx)/d
E = sqrt(4*pi*Z0*Prx_W/Grx)/lambda;

% Convert V/m to dBuV/m
Edbu = 20*log10(E) + 120;