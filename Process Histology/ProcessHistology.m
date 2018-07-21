% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS

% directory of histology
image_folder = 'C:\Drive\Histology\for tutorial - sample data\SS096_raw';

% directory to save the processed images
save_folder = 'C:\Drive\Histology\for tutorial - sample data\SS096_raw';

% name of images, in order anterior to posterior or vice versa
image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'};
                        
% name to save cropped slices as
save_file_name = 'SS096_';

% pixel size parameters
microns_per_pixel = 3.233; %1.62; 
microns_per_pixel_after_downsampling = 10; %to match reference atlas



% additional parameters

% increase gain if for some reason the images are not bright enough
gain = 1; 

% use images that are already at reference atlas (here, 10um/pixel) resolution
use_already_downsampled_image = true; 

% adjust the contrast of the large histology images before cropping them
adjust_histology_contrast = false; 

% size in pixels of reference atlas brain coronal slice -- currently also
% hardcorded into certain parts of the code
reference_size = [800 1140]; 



% find or create folder location for processed images

folder_processed_images = fullfile(save_folder, 'processed');
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end





%% LOAD AND PROCESS SLICE PLATE IMAGES

close all
   

% if adjusting the contrast of each channel, or still need to downsample images
if adjust_histology_contrast || ~use_already_downsampled_image

    % Open figure to view histology
    try; figure(histology_figure);
    catch; histology_figure = figure('Name','Histology Viewer'); end
    warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');
    
    % Downsample and adjust histology image
    HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)
else
    disp('downsampled images already available')
end
  

%% CROP AND SAVE SLICES -- run once the above is done

close all

histology_figure = figure('Name','Histology Viewer');
HistologyCropper(histology_figure, slice_figure, save_folder, image_file_names, reference_size, save_file_name)



%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER

close all
            
slice_figure = figure('Name','Slice Viewer');
SliceFlipper(slice_figure, folder_processed_images, reference_size)

