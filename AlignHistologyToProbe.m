% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME

    
% directory of histology
processed_images_folder = 'C:\Drive\Histology\for tutorial - sample data\Richards_done\processed_3'; 

% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = 'orioll'; % 'compare_angle' 'compare_neither'

% directory of reference atlas files
annotation_volume_location = 'C:\Drive\Histology\for tutorial\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Drive\Histology\for tutorial\structure_tree_safe_2017.csv';
template_volume_location = 'C:\Drive\Histology\for tutorial\template_volume_10um.npy';




%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% create Atals viewer figure
% f = figure('Name','Atlas Viewer'); 
% 
% % show histology in Slice Viewer
% try; figure(slice_figure_browser); title('');
% catch; slice_figure_browser = figure('Name','Slice Viewer'); end
% sliceBrowser(slice_figure_browser, processed_images_folder, f);
% 
% 
% % use application in Atlas Transform Viewer
% f = AtlasTransformBrowser(f, tv,av,st, slice_figure_browser, processed_images_folder, probe_save_name_suffix); % use this function if you a processed_images_folder with appropriately processed .tif histology images


% use the simpler version, which does not interface with processed slice images
% just run these two lines instead of the previous 5 lines of code

save_location = processed_images_folder;
f = allenAtlasBrowser(tv,av,st, save_location, probe_save_name_suffix); % use this function if you lack a processed_images_folder

