function [X,Y,Z] = position(sourceSites, targetSites, varargin)
%position   Position between sites
%   [X,Y,Z] = position(SITE1,SITE2) returns the position in meters of SITE2
%   relative to antenna site SITE1.  The position [X,Y,Z] is expressed in 
%   the local coordinate system of SITE1. For sites with CoordinateSystem 
%   set to 'geographic', the local coordinate system is defined using the
%   East-North-Up (ENU) frame.
%
%   The inputs SITE1 and SITE2 can be scalars or arrays. The outputs are
%   arrays of size M-by-N, where M is the number of sites in SITE2 and N is
%   the number of sites in SITE1.

%   Copyright 2017-2020 The MathWorks, Inc.

% Process optional name/value pairs
p = inputParser;
p.addParameter('Map', []);
p.addParameter('SourceAntennaSiteCoordinates', []);
p.addParameter('TargetAntennaSiteCoordinates', []);
p.parse(varargin{:});

% Get usingCartesian from CoordinateSystem validation or from pre-validated 
% AntennaSiteCoordinates
if isempty(p.Results.SourceAntennaSiteCoordinates)
    usingCartesian = rfprop.internal.Validators.validateCoordinateSystem(sourceSites, targetSites);
else
    usingCartesian = strcmp(p.Results.SourceAntennaSiteCoordinates.CoordinateSystem,'cartesian');
end

% Get antenna site coordinates
if usingCartesian
    map = rfprop.internal.Validators.validateCartesianMap(p);
    coordSys = 'cartesian';
else
    map = rfprop.internal.Validators.validateMapTerrainSource(p, 'position');
    coordSys = 'geographic';
end
rfprop.internal.Validators.validateMapCoordinateSystemMatch(map, coordSys);
sourceCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.SourceAntennaSiteCoordinates, sourceSites, map, 'position');
targetCoords = rfprop.internal.Validators.validateAntennaSiteCoordinates(...
    p.Results.TargetAntennaSiteCoordinates, targetSites, map, 'position');

% Compute ENU coordinates of targets relative to sources
[X,Y,Z] = targetCoords.position(sourceCoords);

