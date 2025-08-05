classdef GeographicSceneModel < handle
%GeographicSceneModel Data container for handling map data import and query
%
%   FOR INTERNAL USE ONLY -- This class is intentionally undocumented
%     and is intended for use only within other toolbox classes and
%     functions. Its behavior may change, or the class itself may be
%     removed in a future release.
%
%   GeographicSceneModel is a handle class meant to allow the import of
%   terrain data and buildings data. It enables querying of map data such
%   as extracting a mesh from a region or extracting profiles over a line.
%   Currently it handles materials-related operations, but this is only
%   temporary and will likely be pulled out in favor of a separate
%   materials management implementation.
%
%   GeographicSceneModel properties:
%       Terrain - The name of the terrain data used (i.e "gmted2010")
%       Materials - A cell array list of chars representing materials
%       Objects - A cell array of structs representing 3D object files added
%       OverridePropagationModelMaterial - Determines if the prop model's
%           material will be overridden by the materials found on the scene's mesh
%       DefaultMaterial - The default material of all meshes added
%       Buildings - An array of OSMBuildings objects in the scene
%
%   GeographicSceneModel methods:
%       GeographicSceneModel - Constructor
%       setTerrainMaterial - Sets the material of the entire terrain data
%       getBuilding - Returns OSMBuilding whose footprint contains input
%       addObject - Adds 3D object to location from input file
%       sceneMesh - Extract Geomtry.mesh from region
%       inBuildings - Returns true if input is within building
%       terrainProfile - Returns terrain heights over a line
%       surfaceProfile - Returns surface heights over a line
%       terrainHeight - Returns terrain height at input point
%       surfaceHeight - Returns surface height at input point

% Copyright 2022-2023 The MathWorks, Inc.

    properties
        Terrain = 'none'
        Objects = {}
        OverridePropagationModelMaterial = false;
        DefaultMaterial = "concrete"
        Buildings = []
    end

    properties(Hidden)
        BuildingsModel = []
        BuildingsMeshData
        BuildingsTriangulation
        MeshWithoutTerrain
        TerrainMaterial = "concrete"
        TerrainSource = 'none'
    end

    methods
        function gsc = GeographicSceneModel(NameValueArgs)
            arguments
                NameValueArgs.Buildings = []
                NameValueArgs.Terrain = 'none'
            end
            if ~strcmp(NameValueArgs.Terrain, 'none')
                gsc.setTerrain(NameValueArgs.Terrain);
            end
            if ~isempty(NameValueArgs.Buildings)
                gsc.setBuildings(NameValueArgs.Buildings);
            end
        end

        function setTerrainMaterial(gsc, material)
            gsc.TerrainMaterial = string(material);
        end

        function bldg = getBuilding(gsc, lat, lon)
            % Iterate over all buildings and get the building whose
            % footprint contains the input
            bldgsArray = gsc.Buildings;
            numBldgs = numel(bldgsArray);
            bldg = [];
            for k = 1:numBldgs
                [isIn, isOn] = isinterior(bldgsArray(k).Footprint, lat, lon);
                if isIn || isOn
                    bldg = bldgsArray(k);
                    break;
                end
            end
        end

        function addObject(gsc, objFile, location, arguments)
            arguments
                gsc
                objFile
                location
                arguments.Material = gsc.DefaultMaterial
                arguments.Scale = 1
                arguments.Rotation = [0, 0, 0]
            end
            [~, ~, ext] = fileparts(objFile);
            objMesh = Geometry.meshRead(objFile, string(arguments.Material));
            objMesh = Geometry.scale(objMesh, arguments.Scale);
            objStruct = struct( ...
                "Mesh", objMesh, ...
                "Location", location, ...
                "Rotation", arguments.Rotation, ...
                "FileName", objFile, ...
                "Extension", ext, ...
                "Scale", arguments.Scale);
            gsc.Objects{end+1} = objStruct;
        end

        function [scMesh, localOrigin] = sceneMesh(gsc, latlim, lonlim, localOrigin, NameValueArgs)
            arguments
                gsc
                latlim = gsc.BuildingsModel.BuildingsLimits(1:2)
                lonlim = gsc.BuildingsModel.BuildingsLimits(3:4)
                localOrigin = [(latlim(1) + latlim(2))/2, (lonlim(1) + lonlim(2))/2, 0]
                NameValueArgs.TerrainResolution = terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution
                NameValueArgs.Materials = dictionary("", gsc.DefaultMaterial)
            end
            latlonlim = [latlim, lonlim];
            % Create meshes from buildings and terrain
            terrainResolution = NameValueArgs.TerrainResolution;
            bldgMesh = gsc.buildingsMesh(latlonlim, localOrigin, NameValueArgs.Materials);
            terrainMesh = gsc.terrainMesh(latlonlim, terrainResolution, localOrigin);

            % Collect any other objects in the scene
            numObjects = numel(gsc.Objects);
            otherMesh = Geometry.mesh;
            if numObjects > 0
                for k = 1:numel(gsc.Objects)
                    objRotation = gsc.Objects{k}.Rotation;
                    % Rotate in ZYX order
                    objMesh = gsc.Objects{k}.Mesh;
                    objMesh = Geometry.rotate(objMesh, [0, 0, 1], objRotation(3)); % Z
                    objMesh = Geometry.rotate(objMesh, [0, 1, 0], objRotation(2)); % Y
                    objMesh = Geometry.rotate(objMesh, [1, 0, 0], objRotation(1)); % X
                    % Translate the mesh to the proper location in the
                    % scene
                    objLoc = gsc.Objects{k}.Location;
                    [x,y,z] = rfprop.internal.MapUtils.geodetic2enu(localOrigin(1), localOrigin(2), localOrigin(3), objLoc(1), objLoc(2), objLoc(3));
                    objMesh = Geometry.translate(objMesh, x, y, z);
                    if strcmp(gsc.Objects{k}.Extension, ".glb")
                        otherMesh = Geometry.join(otherMesh, objMesh);
                    end
                end
            end

            % Combine the meshes into one scene mesh. Save the
            % MeshWithoutTerrain to facilitate operations that do not
            % involve terrain, such as visualization or 3D model
            % generation.
            meshWithoutTerrain = Geometry.join(otherMesh, bldgMesh);
            scMesh = Geometry.join(terrainMesh, meshWithoutTerrain);
            gsc.MeshWithoutTerrain = bldgMesh;
        end

        function isWithinBldg = inbuildings(gsc, lats, lons)
            bldgs = gsc.Buildings;
            % Check if locations are within footprint of any buildings
            isWithinBldg = false(numel(lons),1);
            for bldgInd = 1:numel(bldgs)
                [inp, onp] = bldgs(bldgInd).Footprint.isinterior(lats,lons);
                isWithinBldg = isWithinBldg | inp | onp;
            end
        end

        function [heights, lats, lons] = terrainProfile(gsc, startLatLon, endLatLon, NameValueArgs)
            arguments
                gsc
                startLatLon
                endLatLon
                NameValueArgs.NumQueryPoints = 10
                NameValueArgs.TerrainResolution = terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution
                NameValueArgs.HeightReference = "geoid"
            end
            % Compute intermediary points using LocationResolution and
            % then pass the rest of the arguments to surfaceHeight

            % Compute the distance between sample points using the distance
            % between the input points and the number of query points
            sampleDist = rfprop.internal.MapUtils.greatCircleDistance( ...
                startLatLon(1), startLatLon(2), endLatLon(1), endLatLon(2)) / NameValueArgs.NumQueryPoints;
            [lats, lons] = rfprop.internal.MapUtils.sampleGreatCircle( ...
                startLatLon(1), startLatLon(2), endLatLon(1), endLatLon(2), sampleDist);

            heights = terrainHeight(gsc, lats, lons, "HeightReference", NameValueArgs.HeightReference);
            lats = lats(:);
            lons = lons(:);
        end

        function [heights, lats, lons] = surfaceProfile(gsc, startLatLon, endLatLon, NameValueArgs)
            arguments
                gsc
                startLatLon
                endLatLon
                NameValueArgs.NumQueryPoints = 10
                NameValueArgs.TerrainResolution = terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution
                NameValueArgs.HeightReference = "geoid"
            end
            % Compute intermediary points using LocationResolution and
            % then pass the rest of the arguments to surfaceHeight
            sampleDist = rfprop.internal.MapUtils.greatCircleDistance( ...
                startLatLon(1), startLatLon(2), endLatLon(1), endLatLon(2)) / NameValueArgs.NumQueryPoints;
            [lats, lons] = rfprop.internal.MapUtils.sampleGreatCircle( ...
                startLatLon(1), startLatLon(2), endLatLon(1), endLatLon(2), sampleDist);
            
            % Compute the distance between sample points using the distance
            % between the input points and the number of query points
            heights = surfaceHeight(gsc, lats, lons, "HeightReference", NameValueArgs.HeightReference);
            lats = lats(:);
            lons = lons(:);
        end

        function heights = surfaceHeight(gsc, lats, lons, NameValueArgs)
            arguments
                gsc
                lats
                lons
                NameValueArgs.TerrainResolution = terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution
                NameValueArgs.HeightReference = "geoid"
            end
            % The height of any points within the buildings limits must be
            % computed using the triangulation of that portion of the
            % scene, rather than simply querying the terrain.
            % Because of this, the points within the buildings limits must be
            % determined first.
            bldgsLimits = gsc.BuildingsModel.BuildingsLimits;
            lats = lats(:);
            lons = lons(:);
            indicesOfPointsWithinBldgs = lats >= bldgsLimits(1) & lons >= bldgsLimits(3) & lats <= bldgsLimits(2) & lons <= bldgsLimits(4);

            % Instead of logically indexing, we want to preserve the
            % ordering of all points so we can return the heights in the
            % order that the user gave the points.
            pointsWithinBldgs = [lats .* indicesOfPointsWithinBldgs  lons .* indicesOfPointsWithinBldgs];
            pointsOutsideBldgs = [lats .* ~indicesOfPointsWithinBldgs lons.* ~indicesOfPointsWithinBldgs];
            pointsWithinBldgs(pointsWithinBldgs ==0) = NaN;
            pointsOutsideBldgs(pointsOutsideBldgs==0) = NaN;
            
            % Compute heights for points within buildings region using the
            % triangulation. Start by computing the enu positions.
            bldgCenter = gsc.BuildingsModel.BuildingsCenter;
            [x,y,z] = rfprop.internal.MapUtils.geodetic2enu(bldgCenter(1), bldgCenter(2), bldgCenter(3), pointsWithinBldgs(:,1), pointsWithinBldgs(:,2), 0);

            % Compute normal vectors at locations
            [xp,yp,zp] = rfprop.internal.MapUtils.geodetic2enu(bldgCenter(1), bldgCenter(2), bldgCenter(3), pointsWithinBldgs(:,1), pointsWithinBldgs(:,2), 1);
            xnorm = xp - x;
            ynorm = yp - y;
            znorm = zp - z;

            % Query triangulation surface points
            if ~all(isnan(pointsWithinBldgs))
                mesh = gsc.sceneMesh(bldgsLimits);
                tri = Geometry.meshToTriangulation(Geometry.translate(mesh, 0, 0, -bldgCenter(3)));
                pos = comm.internal.geometry.profileSurface(tri, ...
                    [x y z], [xnorm ynorm znorm]);

                % Compute elevation of surface point using the triangulation.
                [~,~,heightsInBuildings] = rfprop.internal.MapUtils.enu2geodetic(...
                    bldgCenter(1), bldgCenter(2), bldgCenter(3), pos');
            else
                heightsInBuildings = 0;
            end

            
            % The elevation above is with respect to the reference specified by
            % the terrain source. Transform it if different from the
            % requested output height reference.
            terrainSource = gsc.TerrainSource;
            if ~strcmp(terrainSource, "none")
                if strcmpi(terrainSource.HeightReference,'ellipsoid') && strcmpi(NameValueArgs.HeightReference,'geoid')
                    heightsInBuildings = terrain.internal.HeightTransformation.ellipsoidalToOrthometric(heightsInBuildings, lats, lons);
                    heightsInBuildings = terrain.internal.TerrainSource.snapMeanSeaLevel(heightsInBuildings);
                elseif strcmpi(terrainSource.HeightReference,'geoid') && strcmpi(NameValueArgs.HeightReference,'ellipsoid')
                    heightsInBuildings = terrain.internal.HeightTransformation.orthometricToEllipsoidal(heightsInBuildings, lats, lons);
                end
            end
            
            % For all other points that are outside the buildings region,
            % the height can simply be the height of the terrain.
            heightsInTerrain = terrainHeight(gsc, pointsOutsideBldgs(:,1), pointsOutsideBldgs(:,2), "HeightReference", NameValueArgs.HeightReference);
            
            % Since we maintained the order of the points from the first
            % step, all that needs to be done is to interlace all non-nan
            % points together in their existing order.
            heights = heightsInBuildings;
            heights(isnan(heights)) = heightsInTerrain(~isnan(heightsInTerrain));
        end

        function heights = terrainHeight(gsc, lats, lons, NameValueArgs)
            arguments
                gsc
                lats
                lons
                NameValueArgs.HeightReference = "geoid"
            end
            heights = terrain.internal.TerrainSource.queryTerrain(gsc.TerrainSource, [lats(:), lons(:)], NameValueArgs.HeightReference);
        end
    end

    methods(Hidden)
        % This function is used by raytrace to set the material of all
        % buildings if OverridePropagationModelMaterial is false.
        function setBuildingsMaterial(gsc, material)
            material = string(material);
            for k = 1:numel(gsc.Buildings)
                gsc.Buildings(k).Material = material;
            end
        end
        % These functions will be used by siteviewer to add buildings and
        % terrain one after the other during construction.
        function setBuildings(gsc, buildings)
            buildings = globe.internal.Buildings3DModel(buildings, gsc.TerrainSource);
            gsc.BuildingsTriangulation = buildings.Model;
            gsc.BuildingsModel = buildings;
            gsc.Buildings = extractBuildings(gsc.BuildingsModel.BuildingsWithinTerrain);
            numBuildings = numel(gsc.Buildings);
            % This structure is used to filter buildings when specifying a
            % region of interest and apply materials to buildings before
            % constructing a scene mesh
            gsc.BuildingsMeshData = struct("Centroid", zeros(numBuildings, 3), "Mesh", cell(numBuildings, 1));
            for k = 1:numBuildings
                gsc.BuildingsMeshData(k).Centroid = gsc.Buildings(k).Centroid;
                gsc.BuildingsMeshData(k).Mesh = Geometry.mesh(gsc.BuildingsModel.BuildingsTriangulations(k), gsc.DefaultMaterial);
            end
        end

        function setTerrain(gsc, terrainName, terrainObj)
            if nargin < 3
                terrainObj = terrain.internal.TerrainSource.createFromSettings(terrainName);
            end
            gsc.Terrain = terrainName;
            gsc.TerrainSource = terrainObj;
        end
        
        % This function is used by sceneMesh to extract the buildings
        % within the region specified by the user.
        function bldgMesh = buildingsMesh(gsc, latlonlim, regionCenter, materials)
            arguments
                gsc
                latlonlim = gsc.BuildingsModel.BuildingsLimits
                regionCenter = gsc.BuildingsModel.BuildingsCenter
                materials = []
            end
            if isempty(materials)
                % By default, our dictionary will map just the default
                % material
                materials = dictionary("", gsc.DefaultMaterial);
            end
            bldgMesh = Geometry.mesh;
            if ~isempty(gsc.BuildingsModel)
                bldgsLimits = gsc.BuildingsModel.BuildingsLimits;
                % Return an empty mesh if our region doesn't contain
                % the buildings.
                regionExcludesBuildings = (bldgsLimits(1) > latlonlim(2)) || (latlonlim(1) > bldgsLimits(2)) || ...
                    (bldgsLimits(3) > latlonlim(4)) || (latlonlim(3) > bldgsLimits(4));
                if regionExcludesBuildings
                    return
                end

                if materials.numEntries == 1
                    % If only 1 material exists, then instead of creating
                    % the entire buildings mesh from scratch we can re-use
                    % the buildings triangulation as a whole and apply a
                    % single material to the entire mesh.
                    bldgMesh = Geometry.mesh(gsc.BuildingsTriangulation, materials(materials.keys));
                else
                    % If the region includes the buildings, generate the mesh
                    % from each individual building
                    bldgMeshes = gsc.BuildingsMeshData;
                    numMeshes = numel(bldgMeshes);
                    for k = 1:numMeshes
                        % Currently we don't have the ability to 'slice' the
                        % mesh along the rectangle using the Geometry library.
                        % Instead, we must filter the buildings out based on
                        % whether their centroid is contained by the region
                        % limits.

                        % Check if the buildings centroid is within our limits
                        if bldgMeshes(k).Centroid(1) > latlonlim(1) && bldgMeshes(k).Centroid(1) < latlonlim(2) ...
                                && bldgMeshes(k).Centroid(2) > latlonlim(3) && bldgMeshes(k).Centroid(2) < latlonlim(4)
                            % If the buildings material is contained in the
                            % materials dictionary, then it was found in our
                            % rf prop catalog. Otherwise, use the default
                            % material in the dictionary.
                            if materials.isKey(gsc.Buildings(k).Material)
                                material = materials(gsc.Buildings(k).Material);
                            else
                                material = materials("");
                            end
                            % Remove the existing material on the mesh to avoid
                            % previous runs affecting the total materials on
                            % the mesh
                            RemoveMaterial(bldgMeshes(k).Mesh, bldgMeshes(k).Mesh.GetMaterial(1).GetName);
                            AddMaterial(bldgMeshes(k).Mesh, string(material));
                            Geometry.setMaterial(bldgMeshes(k).Mesh, string(material));
                            bldgMesh = Geometry.join(bldgMesh, bldgMeshes(k).Mesh);
                        end
                    end
                end

                % Translate the building mesh up by the buildingscenter
                % amount, as that is the true height of the buildings.
                bldgMesh = Geometry.translate(bldgMesh, 0, 0, gsc.BuildingsModel.BuildingsCenter(3) - regionCenter(3));
                % Translate the building mesh over by the region center
                bldgCenter = gsc.BuildingsModel.BuildingsCenter;
                [dX, dY, ~] = rfprop.internal.MapUtils.geodetic2enu(bldgCenter(1), bldgCenter(2), bldgCenter(3), regionCenter(1), regionCenter(2), regionCenter(3));
                bldgMesh = Geometry.translate(bldgMesh, -dX, -dY, 0);
                % The above transformations ensure that the region center
                % specified by the user is honored by the resulting mesh.
            end
        end

        function terrainMesh = terrainMesh(gsc, latlonlim, resolution, regionCenter)
            arguments
                gsc
                latlonlim = gsc.BuildingsModel.BuildingsLimits;
                resolution = terrain.internal.TerrainSource.MaxBuildingsTerrainTriangulationResolution;
                regionCenter = gsc.BuildingsModel.BuildingsCenter
            end
            terrainTri = terrain.internal.TerrainSource.terrainGridTriangulation(...
                gsc.TerrainSource, ...
                regionCenter(1), regionCenter(2), regionCenter(3), ...
                latlonlim(1:2), latlonlim(3:4), resolution);
            terrainMesh = Geometry.mesh(terrainTri, string(gsc.TerrainMaterial));
        end

    end
end

% Recursively extract buildings from an array of containing potentially
% nested buildings. This is used by the setBuildings function.
function bldgs = extractBuildings(inputBuildings)
bldgs = matlabshared.maps.internal.OSMPrism.empty;
for k = 1:numel(inputBuildings)
    if isprop(inputBuildings(k), "Children")
        bldgs = [bldgs extractBuildings(inputBuildings(k).Children)];
    else
        bldgs(end+1) = inputBuildings(k);
    end
end
end
