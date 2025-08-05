function Prx_db = efield2pow(Edbu, fq, Grx_dbi)
%efield2pow   Convert electric field quantity to power quantity
%   P = rfprop.efield2pow(E,FQ,GRX) returns power quantity in dBm given
%   electric field quantity E in dBuV/m, transmitter frequency FQ, and
%   receiver gain GRX in dBi.

%   Copyright 2017 The MathWorks, Inc.   

% Convert dB quantities to numeric quantities
E = 10^((Edbu-120)/20);
Grx = 10^(Grx_dbi/10);

% Define constants
Z0 = rfprop.Constants.Z0;
lambda = rfprop.Constants.LightSpeed/fq;

% Use equation for E field given power received and gain (see
% sigstrength), but solve for Prx
Prx = ((E^2)*(lambda^2)*Grx)/(4*pi*Z0);

% Convert W to dBm
Prx_db = 10*log10(Prx) + 30;