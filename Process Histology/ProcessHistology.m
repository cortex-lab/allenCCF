% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS


% directory of histology (image with several slices)
% 
% if you already have individual, downsampled slice images, 
% put them in a folder called 'processed' inside of this image_folder
% and run this and the last cell. 
% 
% If you have high-res individual images, put
% them in this image_folder, and just skip the 'crop and save' cell below)
%
image_folder = 'C:\Drive\Histology\for tutorial - sample data\slices';

% directory to save the processed images -- can be the same as the above image_folder
save_folder = 'C:\Drive\Histology\for tutorial - sample data\slices';

% name of images, in order anterior to posterior or vice versa
% once these are downsampled, using the HistologyBrowser function, they
% will be saved in the save_folder, and named the ['original name' '_processed.tif']
image_file_names = dir([image_folder filesep '*.tif']);
image_file_names = natsortfiles({image_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'};

% if the images are individual slices as opposed to an image of multiple
% slices, which must each be cropped and saved
image_file_are_individual_slices = true;

% use images that are already at reference atlas (here, 10um/pixel) resolution
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images below is set to false);
% microns_per_pixel_after_downsampling should be set to 10, to match the atlas
microns_per_pixel = 3.233;
microns_per_pixel_after_downsampling = 10;


% ----------------------
% additional parameters
% ----------------------

% if the images are cropped (image_file_are_individual_slices = false),
% name to save cropped slices as; e.g. the third cropped slice from the 2nd
% image containing many slices will be saved as: save_folder/processed/save_file_name2_3.tif
save_file_name = 'SS096_';

% adjust the contrast of the histology images before cropping them, even if
% they are already downsampled
adjust_histology_contrast = true; 

% increase gain if for some reason the images are not bright enough
gain = 1; 

% size in pixels of reference atlas brain coronal slice
reference_size = [800 1140]; 







% finds or creates a folder location for processed images -- 
% a folder within save_folder called processed
folder_processed_images = fullfile(save_folder, 'processed');
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end


%% LOAD AND PROCESS SLICE PLATE IMAGES

% close all figures
close all
   

% if the images need to be downsampled to 10um pixels (use_already_downsampled_image = false), 
% this will downsample and allow you to adjust contrast of each channel of each image from image_file_names
%
% if the images are already downsampled (use_already_downsampled_image = true), this will allow
% you to adjust the contrast of each channel if adjust_histology_contrast = true
if  adjust_histology_contrast || ~use_already_downsampled_image
    % Open Histology Viewer figure
    try; figure(histology_figure);
    catch; histology_figure = figure('Name','Histology Viewer'); end
    warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');
    
    % Function to downsample and adjust histology image
    HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, folder_processed_images, image_file_are_individual_slices, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)
else
    % if use_already_downsampled_image = true and adjust_histology_contrast = false
    disp('downsampled images already available:')
    disp(image_file_names)
end
  

%% CROP AND SAVE SLICES -- run once the above is done, if image_file_are_individual_slices = false

% close all figures
close all

% run this function if the images from image_file_names have several
% slices, per image (e.g. an image of an entire histology slide), and there 
% are now versions of these images called 'original_file_name_processed.tif'
% in the save_folder -- this function allows you to crop all the slices you 
% would like to process im each image, by drawing rectangles around them in the figure. 
% These can then be further processed in the next cell
if ~image_file_are_individual_slices
    histology_figure = figure('Name','Histology Viewer');
    HistologyCropper(histology_figure, save_folder, image_file_names, reference_size, save_file_name, use_already_downsampled_image)
else
    disp('individually cropped slices already available')
end


%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER

% close all figures
close all
            
% this takes images from folder_processed_images ([save_folder/processed]),
% and allows you to rotate, flip, sharpen, switch order, and further crop them so they
% are in anterior->posterior or posterior->anterior order, and aesthetically pleasing
% note -- presssing left or right arrow saves the modified image, so be
% sure to do this even after modifying the last slice in the folder
slice_figure = figure('Name','Slice Viewer');
SliceFlipper(slice_figure, folder_processed_images, reference_size)


