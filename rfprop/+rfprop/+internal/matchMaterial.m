function matchedMaterial = matchMaterial(material)
%MATCHMATERIAL Match material string with RF Propagation catalog materials
%   MATCHEDMATERIAL = RFPROP.INTERNAL.MATCHMATERIAL(MATERIAL) matches input
%   material text MATERIAL with materials in the RF Propagation catalog.
%   The matching algorithm is insensitive to capitalization, spaces,
%   underscores, and hyphens, and will also look for matches to similar
%   spellings and materials. If the matching algorithm does not find a
%   match, it will attempt two pattern matches including ignoring
%   non-alphabetic characters in MATERIAL. Empty string "" will be used if
%   no unambiguous match is found.
%   
%   MATERIAL must be a character vector, string array, or cell array of
%   character vectors.
%
%   MATCHEDMATERIAL is a string array with size of string(MATERIAL).
%
%   The RF Propagation catalog of materials consists of the materials
%   accepted in the ray tracing propagation model.

% Copyright 2023 The MathWorks, Inc.

% Validate inputs
arguments
    material {mustBeText}
end

% Cast input as string
material = string(material);

% Catalog of materials & substitutes
% Material          Substitute(s)
catalog = {"brick"	[];
"ceiling-board"	    "ceiling panel";
"chipboard"	        ["low-density fiberboard"; "LDF"; "particle board"];
"concrete"	        ["asphalt"; "bitumen"; "blacktop"; "cement"; "pavement"; "tarmac"];
"floorboard"	    [];
"glass"	            [];
"loam"	            ["dirt"; "earth"; "ground"; "soil"];
"marble"            ["granite", "gravel", "slate", "stone", "rock"];
"metal"	            ["aluminium"; "aluminum"; "brass"; "bronze"; "copper"; "gold"; "iron"; "lead"; "nickel"; "silver"; "steel"; "tungsten"; "wire"; "zinc"];
"perfect-reflector"	"pec";
"plasterboard"	    ["drywall"; "gypsum"; "sheetrock"];
"plywood"           [];
"vegetation"	    ["bush"; "flower"; "forest"; "garden"; "grass"; "hedge"; "jungle"; "plant"; "tree"; "woods"];
"water"	            ["lake"; "pond"; "river"];
"wood"	            ["lumber"; "plank"; "timber"]};

% Determine number of substitutes per catalog material
numCatalogMaterials = size(catalog, 1);
numSubstitutePerIdx = zeros(numCatalogMaterials,1);
for mtlIdx = 1:numCatalogMaterials
    numSubstitutePerIdx(mtlIdx) = numel(catalog{mtlIdx,2});
end

% Create list of base catalog materials + substitutes
numCombinedCatalogNames = numCatalogMaterials + sum(numSubstitutePerIdx);
combinedCatalogNames = strings(numCombinedCatalogNames, 1);
combinedCatalogRowIdx = ones(numCombinedCatalogNames, 1);
addIdx = numCatalogMaterials;
for mtlIdx = 1:numCatalogMaterials
    % Insert base material name and its index
    combinedCatalogNames(mtlIdx) = catalog{mtlIdx,1};
    combinedCatalogRowIdx(mtlIdx) = mtlIdx;

    % Add substitute names after all the primary materials, but with a row
    % number referencing the base material's
    for substituteIdx = 1:numSubstitutePerIdx(mtlIdx)
        addIdx = addIdx + 1;
        combinedCatalogNames(addIdx) = catalog{mtlIdx,2}(substituteIdx);
        combinedCatalogRowIdx(addIdx) = mtlIdx;
    end
end

% Strip catalog and input materials' names of capitalization, spaces,
% underscores, and hyphens
combinedCatalogNames = arrayfun(@(x) lower(regexprep(x, [" "; "-"; "_"], '')), combinedCatalogNames);
material = arrayfun(@(x) lower(regexprep(x, [" "; "-"; "_"], '')), material);

% For each input material, find matching catalog material
matchedMaterial = strings(size(material));
for mtlIdx = 1:numel(material)
    thisMaterial = material(mtlIdx);

    % Skip if empty string
    if strcmp(thisMaterial, "")
        continue
    end

    % Check for exact match first, then pattern matching. Ambiguous pattern
    % matches (more than 1 match) will be discarded.
    checkMatch = strcmp(thisMaterial, combinedCatalogNames);
    matchedCatalogRowIdx = unique(combinedCatalogRowIdx(checkMatch));
    numFound = numel(matchedCatalogRowIdx);
    if numFound == 1
        matchedMaterial(mtlIdx) = catalog{matchedCatalogRowIdx, 1};
    else
        % Perform pattern match for material input into catalog list e.g.
        % "concr" in "concrete". Must start at first character.
        checkMatch = arrayfun(@(x) any(strfind(x, thisMaterial) == 1), combinedCatalogNames);
        matchedCatalogRowIdx = unique(combinedCatalogRowIdx(checkMatch));
        numFound = numel(matchedCatalogRowIdx);

        if numFound == 1
            matchedMaterial(mtlIdx) = catalog{matchedCatalogRowIdx, 1};
        elseif numFound == 0
            % Strip input material of non-alphabetic characters
            thisMaterial = regexprep(thisMaterial, '[^a-z]', '');

            % Perform pattern match using catalog list into material input
            % e.g. "concrete" in "concrete123".
            checkMatch = arrayfun(@(x) any(strfind(thisMaterial, x)), combinedCatalogNames);
            matchedCatalogRowIdx = unique(combinedCatalogRowIdx(checkMatch));
            numFound = numel(matchedCatalogRowIdx);

            if numFound == 1
                matchedMaterial(mtlIdx) = catalog{matchedCatalogRowIdx, 1};
            end
        end
    end
end
end

% [EOF]