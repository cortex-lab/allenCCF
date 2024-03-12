function boundary_mask = mask_boundaries(av, margin_radius_um)
% mask_boundaries will return a 3D array of logical which defines the mask
% for structural boundaries
%
% SYNTAX
% mask_boundaries(av, margin_radius_um)
%
% INPUT ARGUMENTS
% av          the 3D array provided from Allen Brain Institute as Allen CCF
%             Grid interval is 10 um (micrometers)
%
% margin_radius_um
%             positive number
%             Boundary margin as a radius in um (micrometers)
%
%
% OUTPUT ARGUMENTS
% boundary_mask  
%             logical
%             The same size as av.
%             Voxels within margin_radius_um from boundaries are labelled
%             as true. Otherwise the values are false.
%
%
% EXAMPLE
% boundary_mask = mask_boundaries(av, 20)
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 16-Aug-2023 11:22:57
%
% See also
% doc


arguments
    av (:, :, :) uint16
    margin_radius_um (1,1) {mustBePositive}
end

% Initialize boundary_mask
boundary_mask = false(size(av));

% Convert margin_um to grid units
margin_units = round(margin_radius_um / 10);

% Iterate through each voxel in av
for x = 1:size(av, 1)
    tic
    for y = 1:size(av, 2)
        for z = 1:size(av, 3)
            % Check neighbors in all three dimensions
            neighbors = av(max(1, x-margin_units):min(size(av, 1), x+margin_units), ...
                          max(1, y-margin_units):min(size(av, 2), y+margin_units), ...
                          max(1, z-margin_units):min(size(av, 3), z+margin_units));
            
            % If any neighbor has a different value, mark as boundary
            if any(neighbors(:) ~= av(x, y, z))
                boundary_mask(x, y, z) = true;
            end
        end
    end
    toc
end