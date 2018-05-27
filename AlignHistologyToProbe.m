% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATION AND PROBE-SAVE-NAME


% directory of histology
processed_images_folder = 'C:\\Users\\Experiment\\Desktop\\brain volumes\\slices\\SS096\\processed\\';
probe_save_name_suffix = 'test';




%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations (takes a while)
if ~exist('av')
    disp('loading reference...')
    av = readNPY('\\ZSERVER.cortexlab.net\Lab\Atlas\allenCCF\annotation_volume_10um_by_index.npy');
    st = loadStructureTree('structure_tree_safe_2017.csv');
end
if ~exist('tv')
    tv = readNPY('\\ZSERVER.cortexlab.net\Lab\Atlas\allenCCF\template_volume_10um.npy');
end

% show histology in Slice Viewer
try; figure(slice_figure); title('');
catch; slice_figure = figure('Name','Slice Viewer'); end
sliceBrowser(slice_figure, processed_images_folder);

% use application in Atlas Viewer
% f = AtlasTransformBrowser(tv,av,st, slice_figure, processed_images_folder, probe_save_name_suffix); % use this function if you a processed_images_folder with appropriately processed .tif histology images
save_location = processed_images_folder;
f = allenAtlasBrowser(tv,av,st, save_location, probe_save_name_suffix); % use this function if you lack a processed_images_folder
   

