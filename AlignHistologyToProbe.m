% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME


% directory of histology
processed_images_folder = 'C:\Drive\Histology\for tutorial\SS096\processed\\';

probe_save_name_suffix = 'test'; % name the saved probe points, to avoid overwriting another set of probes going in the same folder





%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations (takes a while)
if ~exist('av')
    disp('loading reference...')
    av = readNPY('C:\Drive\Histology\for tutorial\annotation_volume_10um_by_index.npy');
    st = loadStructureTree('C:\Drive\Histology\for tutorial\structure_tree_safe_2017.csv');
end
if ~exist('tv')
    tv = readNPY('C:\Drive\Histology\for tutorial\template_volume_10um.npy');
end

% show histology in Slice Viewer
try; figure(slice_figure2); title('');
catch; slice_figure2 = figure('Name','Slice Viewer','Position', [121 542 822 542]); end
sliceBrowser(slice_figure2, processed_images_folder);

% use application in Atlas Viewer
f = AtlasTransformBrowser(tv,av,st, slice_figure2, processed_images_folder, probe_save_name_suffix); % use this function if you a processed_images_folder with appropriately processed .tif histology images
% save_location = processed_images_folder;
% f = allenAtlasBrowser(tv,av,st, save_location, probe_save_name_suffix); % use this function if you lack a processed_images_folder
   
