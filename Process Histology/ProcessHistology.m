% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS

% directory of histology
image_folder = '\\zubjects\\Subjects\\Richards\\Histology\\';
save_folder = 'P:\brain volumes\slices\Richards\tutorial\';



% name of images, in order anterior to posterior or vice versa
image_file_names = {'Large Image 1.tif','Large Image 2.tif','Large Image 13.tif'}; %,'Large Image 12.tif',...
%                              'Large Image 11.tif', 'Large Image 9.tif','Large Image 8.tif','Large Image 7.tif','Large Image 6.tif',...
%                              'Large Image 10.tif','Large Image 5.tif','Large Image 4.tif'};
                            % note that image 14/15/16/17, which fit inbetween 2 and 12, have no dye visible
                        
% name to save cropped slices as
save_file_name = 'Tatum_';

% parameters
microns_per_pixel = 3.233; %1.62; 
microns_per_pixel_after_downsampling = 1;



% additional parameters
gain = 8; % increase gain if for some reason the images are not bright enough
use_already_downsampled_image = false; 

adjust_slice_contrast = true; % adjust large histology image before cropping




close all

% find or create folder location for processed images
reference_size = [800 1140]; % size in of reference atlas brain coronal slice
adjust_histology_contrast = true; % adjust image before moving on
folder_processed_images = [save_folder 'processed\\'];
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end; 



%% LOAD AND PROCESS SLICE PLATE IMAGES




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




try; figure(slice_figure);
catch; slice_figure = figure('Name','Slice Viewer'); end

try; figure(histology_figure);
catch; histology_figure = figure('Name','Histology Viewer'); end




HistologyCropper(histology_figure, slice_figure, save_folder, image_file_names, reference_size, save_file_name)






%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER




try; figure(slice_figure);
catch; slice_figure = figure('Name','Slice Viewer'); end



SliceFlipper(slice_figure, folder_processed_images, reference_size)








