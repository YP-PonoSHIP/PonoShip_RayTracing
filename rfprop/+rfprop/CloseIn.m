classdef CloseIn < rfprop.PropagationModel
%

% Copyright 2017-2024 The MathWorks, Inc. 

% Reference: S. Sun, T.S. Rappaport, T. Thomas, A. Ghosh, H. Nguyen, I.
% Kovacs, I. Rodriguez, O. Koymen, A. Partyka, "Investigation of prediction
% accuracy, sensitivity, and parameter stability of large-scale propagation
% path loss models for 5G wireless communications," IEEE Transactions on
% Vehicular Technology, vol. 65, no. 5, pp. 2843 - 2860, May 2016.

% Default property values obtained from Table III ibid.
    properties
        ReferenceDistance(1,1) double {mustBeFinite, mustBeReal, mustBePositive, mustBeNonsparse} = 1 % m
        PathLossExponent(1,1) double {mustBeFinite, mustBeReal, mustBeNonsparse} = 2.9
        Sigma(1,1) double {mustBeFinite, mustBeReal, mustBePositive, mustBeNonsparse} = 5.7 % dB
        NumDataPoints(1,1) double {mustBeFinite, mustBeReal, mustBePositive, mustBeInteger, mustBeNonsparse} = 1869
    end
    
    methods(Access = protected)
        function pl = pathlossOverDistance(pm, ~, tx, d, ~)
            
            % Get CI model parameters
            fq = tx.TransmitterFrequency;
            d0 = pm.ReferenceDistance;
            n = pm.PathLossExponent;
            sigma = pm.Sigma;
            N = pm.NumDataPoints;
            
            % Compute path loss
            if d < d0
                pl = 0;
            else
                pl = cipl(fq,d0,d,n,sigma,N);
            end
        end
    end
end

function PL = cipl(f, d0, d, n, sigma, N)

c = rfprop.Constants.LightSpeed;
lambda = c/f;

% The terms
term1 = fspl(d0,lambda);
term2 = 10*n*log10(d./d0);

% Gaussian RV - sigma is voltage, var is power
sigma_lin = 10^(sigma/20);
Xi = sigma_lin.*randn(N,1);
Xi = Xi.^2;

PL_rv = term1 +term2 + 10*log10(Xi);
PL = mean(PL_rv);
end