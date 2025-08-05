classdef RayInteraction < comm.internal.ConfigBase 
%RAYINTERACTION Ray-surface interaction
%   RAYINT = COMM.INTERNAL.RAYINTERACTION creates a ray-surface
%   interaction object, RAYINT, for ray tracing applications. It describes
%   one interaction of a ray with the environmental surface.
% 
%   RAYINT = COMM.INTERNAL.RAYINTERACTION(Name,Value) creates a ray-surface
%   interaction object, RAYINT, with the specified property Name set to the
%   specified value. You can specify additional name-value pair arguments
%   in any order as (Name1,Value1, ...,NameN,ValueN).
%    
%   COMM.INTERNAL.RAYINTERACTION properties:
%
%   Type             - Interaction type
%   Location         - Interaction location
%   Material         - Interaction material
%   AngleOfIncidence - Angle of incidence
%
%   References: 
%   [1] D.A. McNamara, C.W.I. Pistorius, and J.A.G. Malherbe,
%   Introduction to the Uniform Geometrical Theory of Diffraction, Artech
%   House, 1990.

% Copyright 2021-2024 The MathWorks, Inc.
%codegen

properties
    %Type Interaction type
    %   Specify the ray-surface interaction type as one of "Reflection" |
    %   "Diffraction". The default value of this property is "Reflection".
    Type = 'Reflection'
    %Location Interaction location
    %   Specify the ray-surface interaction location as a double precision,
    %   3-by-1 vector. In a geographic coordinate system defined with
    %   reference to Earth ellipsoid model WGS-84, this property represents
    %   the location in the form of [latitude; longitude; height]. The
    %   latitude coordinate in degrees must be in the range [-90, 90]. In a
    %   Cartesian coordinate system, this property represents the location
    %   in the form of [x; y; z]. The default value of this property is
    %   [10; 10; 0] in a Cartesian system.
    Location = [10; 10; 0]
    %Material Interaction material
    %   Specify the ray-surface interaction material as a 2x1 vector for a
    %   reflection and 2x2 matrix for a diffraction. The first row is for
    %   permittivity and the second row is for conductivity. The first and
    %   second columns are for the o-face and n-face respectively for a
    %   diffraction. The default value of this property is [NaN; NaN] for
    %   an unknown material.
    Material = [nan; nan];
    %AngleOfIncidence Angle of incidence 
    %   Specify the angle of incidence in degrees as a scalar for a
    %   reflection and 1-by-2 vector for a diffraction. The first and
    %   second elements are the angle of incidence upon the o-face and
    %   n-face respectively for a diffraction.
    AngleOfIncidence
    %DiffractionEdgeDirection Diffraction edge direction
    %   Specify the unit edge direction as a 3-by-1 vector for a
    %   diffraction. This properyt applies when Type is set to
    %   'Diffraction'.
    DiffractionEdgeDirection = ones(3, 1)
    %DiffractionFaceNormal Diffraction face normal
    %   Specify the normals of the two faces that form the diffraction edge
    %   as a 3-by-2 matrix. The first and second columns are for the o-face
    %   and n-face respectively. This properyt applies when Type is set to
    %   'Diffraction'.
    DiffractionFaceNormal = ones(3, 2)  
    %DiffractionFaceCenter Diffraction face center
    %   Specify the in-center of the two faces that form the diffraction
    %   edge as a 3-by-2 matrix. The first and second columns are for the
    %   o-face and n-face respectively. This properyt applies when Type is
    %   set to 'Diffraction'.
    DiffractionFaceCenter = ones(3, 2)    
    %PrimitiveID Primitive ID
    %   Specify the primitive ID of the interacting surface as a positive,
    %   integer scalar for a reflection and 1-by-2 vector for a
    %   diffraction. The default value of this property is 1.
    PrimitiveID = 1
end

properties (Hidden)
    MaterialName
end

properties(Constant, Hidden)
    Type_Values  = {'Transmission', 'Reflection', 'Diffraction'};
end

methods
    function obj = RayInteraction(varargin)
        obj@comm.internal.ConfigBase(varargin{:});
    end
end

methods (Access = protected)
    function flag = isInactiveProperty(obj, prop)
        if any(strcmp(prop, {'DiffractionEdgeDirection','DiffractionFaceNormal', ...
                'DiffractionFaceCenter'}))
            flag = ~strcmp(obj.Type, 'Diffraction');
        else
            flag = false;
        end
    end
end

methods(Static, Hidden)
    function obj = createRefPath(refLoc, refMtl, srcLoc, dstLoc)
        numRef = size(refLoc, 2);

        % Calculate angle of incidence for each reflection
        path = [srcLoc refLoc dstLoc];
        pathDir = diff(path, 1, 2);
        pathDir =  pathDir./ sqrt(sum(pathDir.^2, 1));
        aoi = acosd(dot(pathDir(:,2:end), -pathDir(:,1:end-1)))/2;

        % Find segments that have 0 incident angle
        idx = (sum((pathDir(:,2:end) + pathDir(:,1:end-1)).^2) < eps);
        if any(idx)
            aoi(idx) = 0;
        end

        % Create and assign ray interaction objects
        obj(1, numRef) = comm.internal.RayInteraction;
        for i = 1:numRef
            obj(1, i).Location = refLoc(:, i);
            if size(refMtl,2) > 1
                obj(1, i).Material = refMtl(:, i);
            else
                obj(1, i).Material = refMtl(:, 1);
            end
            obj(1,i).AngleOfIncidence = aoi(i);
        end
    end
end

end