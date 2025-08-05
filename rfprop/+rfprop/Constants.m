classdef Constants
%Constants   RF propagation constants

%   Copyright 2017-2021 The MathWorks, Inc.
    
    properties(Constant)
        LightSpeed = physconst('lightspeed');
        SphericalEarthRadius = earthRadius; % Earth radius using spherical model
        MinorEarthRadius = 6356752; % Minor Earth radius of WGS-84 ellipsoid model
        Z0 = 120 * pi; % Impedance of free space
        MaxNumSitesAutoResolutionGrid = 160000; % Maximum number of sites in grid using automatic resolution (same size as 400-by-400)
        MaxNumSitesShow = 25000; % Maximum number of sites in an array to show
        RequestTimeout = 60; % Maximum time in seconds to wait for Site Viewer requests
        MaxPropagationDistance = earthRadius; % Maximum allowed distance for propagation modeling
        MaxPropagationDistanceUsingTerrain = maxPropagationDistanceOverTerrain; % Maximum great circle distance (m) for terrain propagation modeling
        % Maximum distance to use terrain sampling to compute los instead
        % of distance-to-horizon approximation
        MaxLOSTerrainSampleDistance = earthRadius;
        MaxLOSNumTerrainSamples = 10000; % Maximum number of terrain samples 
        DefaultLOSTerrainResolution = 30; % Default terrain sampling distance
        MaxLatLonInfoDecimalPlaces = 6; % Maximum number of decimal places to show for info of lat/lon values
        MaxDistanceInfoDecimalPlaces = 2; % Maximum number of decimal places to show for info of distance values
        MaxPhaseInfoDecimalPlaces = 2; % Maximum number of decimal places to show for info of phase change values
        MaxAngleInfoDecimalPlaces = 1; % Maximum number of decimal places to show for info of angle values
        MaxPowerInfoDecimalPlaces = 1; % Maximum number of decimal places to show for info of power values
        MapTileCacheSize = 1000; % Cache size of image and terrain tiles
        DefaultMaxRange = 30000; % Default MaxRange distance (m)
        DefaultMultipathModelMaxRange = 500; % Default MaxRange distance when multipath model is used (m)
        TerrainProfileResolutionWithBuildings = 10; % Resolution for terrain profiles when Site Viewer contains buildings (m)
        DefaultMaxImageSize = 16384; % Default MaxImageSize for use without Site Viewer
        ProgressDialogUpdateSkipFactor = 3; % Skip factor for updating progress dialog
        MinNumRadialsFullProgressDialog = 600; % Minimum number of terrain radials for showing full progress dialog
        PhasedPatternResolution = struct("low", 5, "medium", 3, "high", 1);
        AntennaPatternResolution = struct("low", 15, "medium", 7.5, "high", 'auto');
        AntennaPatternInterpolation = struct("low", 10, "medium", 5, "high", 1);
        IsotropicPatternResolution = struct("low", 10, "medium", 5, "high", 1);
        GainAntennaPatternResolution = 1;
        ShowLocationBuildingSideThreshold = 1;
        DefaultRaytracingBuildingsMaterial = "concrete";
        DefaultRaytracingSurfaceMaterial = "concrete";
    end
end

function d = maxPropagationDistanceOverTerrain

s = settings;
d = s.shared.channel.rfprop.MaxPropagationDistanceOverTerrain.ActiveValue;
end