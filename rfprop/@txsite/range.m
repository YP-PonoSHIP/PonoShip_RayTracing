function maxrange = range(tx, ss, Grx_db, type, pm)
%

% Copyright 2017-2014 The MathWorks, Inc.

% Convert electric field strength quantity to power quantity
if strcmp(type, 'efield')
    ss = rfprop.efield2pow(ss, tx.TransmitterFrequency, Grx_db);
end

% Calculate path loss to reach minimum signal level. Use Friis equation as
% defined in sigstrength and solve for Lpl:
%
% Lpl_db = Ptx_db + Gtx_db + Grx_db - Prx_db
%
% where:
%  * Prx_db is received power in dBm at receiver input
%  * Ptx_db is transmitter output power in dBm
%  * Gtx_db is transmitter gain in dBi (antenna gain - system loss)
%  * Grx_db is receiver gain in dBi (antenna gain - system loss)
%  * Lpl_db is path loss (dB) as given by propagation model

Prx_db = ss;
Ptx = tx.TransmitterPower; % Units: W
Ptx_db = 10*log10(1000 * Ptx); % Convert W to dBm
Gtx_db = gain(tx, tx.TransmitterFrequency) - tx.SystemLoss; % Units: dBi
pl = Ptx_db + Gtx_db + Grx_db - Prx_db;

% Get maximum range from propagation model, which computes distance to
% reach path loss. Any greater than this distance will result in path loss
% that produces less than minimum signal level.
maxrange = pm.range(tx, pl);

% Saturate maxrange at diameter of Earth
earthDiameter = 2*rfprop.Constants.SphericalEarthRadius;
if maxrange > earthDiameter
    maxrange = earthDiameter;
end