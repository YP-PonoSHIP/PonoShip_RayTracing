function [refPts, refFacetIdx] = calcImageReflectionPath(env, txPos, ...
    visFacetTx, rxPos, currentOrder, maxOrder)
%CALCREFLECTIONPATH Calculate reflection paths
%   [REFPTS, REFFACETIDX] = COMM.INTERNAL.GEOMETRY.CALCREFLECTIONPATH(ENV,
%   RT, TXPOS, VISFACETTX, RXPOS, CURRENTORDER, MAXORDER) calculates the
%   reflection paths with order of CURRENTORDER, from the transmitter
%   position, TXPOS, to the receiver position, RXPOS, given an an
%   environmental description, ENV. TXPOS/RXPOS is a double precision,
%   1-by-3 vector representing Tx/Rx position in the form of [x y z].
%   VISFACETTX is an integer column vector containing the indices of the
%   facets which are visible from TXPOS. MAXORDER (>= CURRENTORDER) is the
%   order of reflection being search by a ray-tracer.
%   
%   ENV is a structure with three fields: 
%     Triangulation:  Triangulation object for the environmental geometry
%     RayTracer:      Ray tracer object built from the triangulation object
%     SharpEdgeFlags: Logical matrix of the same size as
%                     Triangulation.ConnnectivityList indicating each edge
%                     of every triangle is sharp or not
%
%   REFPTS is a double precision, NRAY-by-3-by-CURRENTORDER array, where
%   NRAY is the number of rays being found. The 2nd dimension is the
%   reflection point positions in the form of [x y z]. REFFACETIDX is an
%   integer, NRAY-by-CURRENTORDER matrix.
%
%   This function is recursively called to find all the reflection paths.
%   When CURRENTORDER < MAXORDER, it returns all the "candidate" rays that
%   have collected the last CURRENTORDER many of reflection points before
%   Rx. When CURRENTORDER == MAXORDER, it returns all the reflection paths
%   from Tx to Rx. 

% Copyright 2019-2022 The MathWorks, Inc.

% Parse env input
TRI = env.Triangulation;
RT  = env.RayTracer;
excludeSharpEdges = true;

% Get Tx images against all visible facets
txImages = comm.internal.geometry.mirror(TRI, txPos, visFacetTx);

if currentOrder == 1 
    % This is the last leg of the rays connecting Rx. Shoot a ray from Rx
    % to all Tx images. No intersection point on a sharp edge. 
    [sPts, sFacetIdx] = comm.internal.geometry.firstIntersect( ...
        env, rxPos, txImages, 'segment', excludeSharpEdges);

    % Find the candidate rays whose last intersection facet is the same as
    % the facet the Tx is imaged against. 
    candRayIdx = (sFacetIdx == visFacetTx); 
    
    if any(candRayIdx) % Extract candidate rays
        % Get the last reflection point for candidate rays. It is of size
        % [numRay, 3, numRef] with numRef = 1
        refPts = sPts(candRayIdx, :); 
        % Get the last reflection facet index for candidate rays. It is of
        % size [numRay, numRef] with numRef = 1 here
        refFacetIdx = sFacetIdx(candRayIdx);
    else % No candidate ray
        refPts = zeros(0, 3, 1);
        refFacetIdx = zeros(0, 1); 
    end
else % A middle leg of the ray
    if currentOrder > 2 
        visFacetTxImage = comm.internal.geometry.viewshed(TRI, RT, ...
            incenter(TRI, visFacetTx), visFacetTx);
        visFacetIdx = 1:length(visFacetTx);
    else % currentOrder == 2 
        % For the second last leg of rays, we can use the single viewshed
        % of Rx instead of multiple viewsheds of all Tx's
        % (RefOrder-1)-order images.
        visFacetTxImage = comm.internal.geometry.viewshed(TRI, RT, rxPos);
        visFacetIdx = ones(1, length(visFacetTx));
    end
    
    % Initialize returns
    refPts = zeros(0, 3, currentOrder);
    refFacetIdx = zeros(0, currentOrder); 

    % Iterate through all the visible facets from Tx
    for i = 1:length(visFacetTx)
        if ~isempty(visFacetTxImage{visFacetIdx(i)})
            % Get all candiate rays from this Tx image to Rx. thisRefPts is
            % of size [numRay, 3, currentOrder-1] and thisRefFacetIdx is of
            % size [numRay, currentOrder-1].
            [thisRefPts, thisRefFacetIdx] = ... 
                comm.internal.geometry.calcImageReflectionPath(...
                env, txImages(i,:), visFacetTxImage{visFacetIdx(i)}, ...
                rxPos, currentOrder-1, maxOrder);
            
            if isempty(thisRefPts)
                continue; 
            end

            % Shoot a ray from the last reflection point of all candiate
            % rays to the Tx image. No intersection point on a sharp edge. 
            [sPts, sFacetIdx] = comm.internal.geometry.firstIntersect( ...
                env, thisRefPts(:,:,1), txImages(i,:), 'segment', ...
                excludeSharpEdges);

            % The "surviving" candidate rays should be those that have the
            % intersection facet the same as what the Tx image is against.
            candRayIdx = (sFacetIdx == visFacetTx(i, :)); 
                        
            if ~any(candRayIdx) 
                continue;  
            end
            
            % For the "surviving" candidate rays, get the intersection
            % point and facet index. 
            newRefPts = sPts(candRayIdx, :);        % [numRay, 3, 1]
            newRefFacetIdx = sFacetIdx(candRayIdx); % [numRay, 1]
        
            % Update the last reflection points and facet indices for the
            % "surviving" candidate rays. The point is of size [numRay, 3,
            % currentOrder]; The facet index is of size [numRay,
            % currentOrder].
            thisRefPts = cat(3, newRefPts, thisRefPts(candRayIdx,:,:));
            thisRefFacetIdx = [newRefFacetIdx, thisRefFacetIdx(candRayIdx, :)];
            
            % Add candiate rays from this iteration to the overall
            % candidate ray repository.
            refPts = cat(1, refPts, thisRefPts);
            refFacetIdx = cat(1, refFacetIdx, thisRefFacetIdx);
        end
    end
end    

if (currentOrder == maxOrder) && ~isempty(refPts)
    % This is the first leg of the rays connecting Tx. Shoot a ray from Tx
    % to the last reflection point of all candidate Rays. If it is LOS, the
    % ray is "surviving".
    [~, occFacetIdx, isNLOS] = comm.internal.geometry.firstIntersect(env, ...
        txPos, refPts(:,:,1), 'segment');
    noOcclusion = ~isNLOS | refFacetIdx(:,1) == occFacetIdx;

    % Only return those "surviving" rays
    refPts = refPts(noOcclusion,:,:);
    refFacetIdx = refFacetIdx(noOcclusion,:,:);
end

end

% [EOF]