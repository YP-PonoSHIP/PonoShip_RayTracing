function [MagEBlocks,theta1,phi1] = rearrangeData(MagE,theta,phi)
%rerrangeData takes the MagE, theta and phi in the user input form and
%rearranges this data so that it can be plotted. This functions takes care
%of the folllowing 4 cases:
%
% 1) Input data is in correct format (unscrambled and arranged)
% 2) Input data is in scrambled form
% 3) There are some missing data points in the input data vectors
% 4) There are some extra data points for repeated combinations of theta
%    and phi.

% Both theta and phi need to be column vectors

%   Copyright 2015-2020 The MathWorks, Inc.
if isrow(theta)
    theta = theta.';
end

if isrow(phi)
    phi = phi.';
end

% Checking if the data is in vector form with correct dimensions
if (numel(MagE) == numel(phi)) && (numel(MagE) == numel(theta))

    % Force the inputs to be column vectors
    MagE = reshape(MagE,[length(MagE),1]);
    phi = reshape(phi,[length(phi),1]);
    theta = reshape(theta,[length(theta),1]);

    % Taking the unique values for theta and phi
    [theta1] = unique(theta,'sorted');
    [phi1] = unique(phi,'sorted');
    sizeThetaOrder = length(theta1);
    sizePhiOrder = length(phi1);

    % Determine the total size of the MagE vector
    TotalSize = sizeThetaOrder*sizePhiOrder;

    % Test for repeated values
    phithetacomb = [phi, theta];
    [~,ic,~] = unique(phithetacomb,'rows');
    dupeinds = setdiff(1:size(phithetacomb,1),ic);
    dupemat = phithetacomb(dupeinds,:);

    % Finding if the MagE values are same for these repeated points
    [~,ib,~] = intersect(phithetacomb,dupemat,'rows','stable');
    Mag1 = MagE(dupeinds);
    Mag2 = MagE(ib);

    [~,ig] = setdiff(Mag1,Mag2,'stable');
    diffpoint = dupeinds(ig);

    [~,ih] = setdiff(dupeinds,diffpoint,'stable');
    samepoint = dupeinds(ih);

    % Fixing the input vectors based on the results
    MagE(samepoint) = [];
    theta(samepoint) = [];
    phi(samepoint) = [];
    dupemat(ih,:) = [];

    % Erroring out if different values of MagE
    if ~isempty(dupemat)
        error(message('shared_channel:patternCustom:RepeatedPoints',dupemat(1,1),dupemat(1,2)));
    end

    % If the sizes match, this loop is executed
    if (TotalSize == length(MagE))
        ThetaBlocks = reshape(theta,[sizePhiOrder, sizeThetaOrder]);
        PhiBlocks = reshape(phi,[sizePhiOrder, sizeThetaOrder]);

        % Producing Phi and Theta for comparison
        ThetaBlocks1 = (repmat(theta1,1,sizePhiOrder))';
        PhiBlocks1 = repmat(phi1,1,sizeThetaOrder);

        % If they are same, then this loop is exectued
        if (isequal(ThetaBlocks,ThetaBlocks1) && isequal(PhiBlocks,PhiBlocks1))
            MagEBlocks = reshape(MagE,[sizePhiOrder, sizeThetaOrder]);

            % This loop is executed, if they are not same
        else
            phithetacomb = [phi, theta];
            [~,idx] = sortrows(phithetacomb,2);
            MagEBlocks = MagE(idx);
            MagEBlocks = reshape(MagEBlocks,[sizePhiOrder, sizeThetaOrder]);
        end

    else
        % Initializing the MagE matrix to be of appropriate size
        MagEBlocks = NaN(TotalSize,1);

        % Putting the data in the arranged format
        phithetacomb = [phi, theta];

        thetarep = sort(repmat(theta1,[sizePhiOrder,1]));
        phirep = repmat(phi1,[sizeThetaOrder,1]);

        [~,idx] = sortrows(phithetacomb,2);
        MagE = MagE(idx);

        combophitheta = [phirep, thetarep];
        [~,ia,~] = intersect(combophitheta,phithetacomb,'rows','stable');
        MagEBlocks(ia) = MagE;

        MagEBlocks = reshape(MagEBlocks,[sizePhiOrder, sizeThetaOrder]);
    end
else
    % Data is already in intended form
    % Ensure the data is increasing
    [theta1, Itheta] = sort(theta);
    [phi1, Iphi]     = sort(phi);
    MagEBlocks       = MagE(Iphi, Itheta);
end