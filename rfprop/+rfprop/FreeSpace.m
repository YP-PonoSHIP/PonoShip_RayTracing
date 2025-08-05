classdef FreeSpace < rfprop.PropagationModel
%

% Copyright 2017-2024 The MathWorks, Inc.   

    methods(Access = protected)
        function pl = pathlossOverDistance(~, ~, tx, d, ~)
            lambda = rfprop.Constants.LightSpeed/tx.TransmitterFrequency;
            pl = fspl(d,lambda);
        end
    end
    
    methods
        function r = range(pm, txs, pl)
            %
            
            % Override range since model has equation to compute instead of
            % using fzero to find answer
            
            % Validate inputs
            validateattributes(pm,{'rfprop.FreeSpace'}, ...
                {'scalar'},'range','',1);
            validateattributes(txs,{'txsite'},{'nonempty'},'range','',2);
            validateattributes(pl,{'numeric'}, ...
                {'real','finite','nonnan','nonsparse','scalar'},'range','',3);
            
            numTx = numel(txs);
            r = zeros(numTx,1);
            for txInd = 1:numel(txs)
                tx = txs(txInd);
                maxrange = pm.fsrange(tx.TransmitterFrequency, pl);
                r(txInd) = validateRange(pm, maxrange);
            end
        end
    end
    
    % Use static method for range computation so it can be called by
    % rfprop.PropagationModel.range
    methods(Static, Hidden)
        function d = fsrange(fq,pl)
            lambda = rfprop.Constants.LightSpeed/fq;
            d = (lambda*(10^(pl/20)))/(4*pi);
        end
    end
end