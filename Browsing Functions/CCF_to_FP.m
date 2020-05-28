function [ann_FP, name_FP, acr_FP] = CCF_to_FP(ccf_ap, ccf_dv, ccf_ml, CCFtoFPtable, FPtable, chon_images_loc)
% convert Allen CCF coordinates to Franklin-Paxinos coordinates, using data
% from Chon et al. Enhanced and unified anatomical labeling for a
% common mouse brain atlas (2020).


% argmin to get nearest slice
[m, idx] = min(abs(ccf_ap - CCFtoFPtable.slice_num));

% get the label image name
curr_label_file = char(CCFtoFPtable.label_file(idx));

% load the image
curr_label_image = imread(fullfile(chon_images_loc, strcat(curr_label_file, '.tif')));

% testing - show the image
% imshow(curr_label_image);

% find the label at the particular coordinates
curr_label = curr_label_image(ccf_dv, ccf_ml);

% if it misses the brain -- try second closest slice
if ~curr_label
    % if FP slice is after (more posterior), get the previous slice
    if (CCFtoFPtable.slice_num(idx) - ccf_ap) >= 0
        idx = idx - 1;
    % otherwise get the next (anterior) slice
    else
        idx = idx + 1;
    end

    % get the label image name
    curr_label_file = char(CCFtoFPtable.label_file(idx));

    % load the image
    curr_label_image = imread(fullfile(chon_images_loc, strcat(curr_label_file, '.tif')));

    % find the label at the particular coordinates
    curr_label = curr_label_image(ccf_dv, ccf_ml);

end
    
% get idx from the FP table
curr_idx = find(FPtable.Structural_ID==curr_label);

try
    % get annotation label
    ann_FP = curr_label;

    % get full name
    name_FP = FPtable.Franklin_Paxinos_Full_name{curr_idx};

    % get abbreviated name
    acr_FP = FPtable.Franklin_Paxinos_abbreviation{curr_idx};
catch
    % get annotation label
    ann_FP = curr_label;

    % get full name
    name_FP = 'not found';

    % get abbreviated name
    acr_FP = 'not found';
end


