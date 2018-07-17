% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS

% directory of histology
image_folder = 'C:\Drive\Histology\for tutorial\SS096_raw\';
save_folder = 'C:\Drive\Histology\for tutorial\SS096_raw\';



% name of images, in order anterior to posterior or vice versa
image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'};
                        
% name to save cropped slices as
save_file_name = 'SS096_';

% parameters
microns_per_pixel = 3.233; %1.62; 
microns_per_pixel_after_downsampling = 10;



% additional parameters
gain = 1; % increase gain if for some reason the images are not bright enough
use_already_downsampled_image = true; 

adjust_slice_contrast = true; % adjust large histology image before cropping






% find or create folder location for processed images
reference_size = [800 1140]; % size in of reference atlas brain coronal slice
adjust_histology_contrast = true; % adjust image before moving on
folder_processed_images = [save_folder 'processed\\'];
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end



%% LOAD AND PROCESS SLICE PLATE IMAGES

close all


% open figure to view histology
try; figure(histology_figure);
catch; histology_figure = figure('Name','Histology Viewer'); end
warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');


    
% adjust contrast if desired
file_num = 1; % start with this file
if adjust_histology_contrast
    HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, file_num, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)
end
  




%% CROP AND SAVE SLICES -- only run once the above is done

close all


slice_figure = figure('Name','Slice Viewer'); 
histology_figure = figure('Name','Histology Viewer');




HistologyCropper(histology_figure, slice_figure, save_folder, image_file_names, reference_size, save_file_name)






%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER

close all


slice_figure = figure('Name','Slice Viewer');



SliceFlipper(slice_figure, folder_processed_images, reference_size)








