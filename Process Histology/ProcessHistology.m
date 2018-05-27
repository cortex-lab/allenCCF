% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS

% directory of histology
image_folder = '\\zubjects\\Subjects\\Richards\\Histology\\';
save_folder = 'C:\\Users\\Experiment\\Desktop\\brain volumes\\slices\\Richards\\';
% name of images, in order anterior to posterior
image_file_names = {'Large Image 1.tif','Large Image 2.tif','Large Image 3.tif','Large Image 4.tif',...
                            'Large Image 5.tif','Large Image 6.tif','Large Image 7.tif','Large Image 8.tif',...
                            'Large Image 9.tif','Large Image 10.tif','Large Image 11.tif','Large Image 12.tif',...
                            'Large Image 13.tif','Large Image 14.tif','Large Image 15.tif','Large Image 16.tif','Large Image 17.tif'}

% name to save cropped slices as
save_file_name = 'Richards_';

% parameters
gain = 10;
use_already_downsampled_image = false;
microns_per_pixel = 3.233; %1.62; 
microns_per_pixel_after_downsampling = 10;
adjust_histology_contrast = true; 
adjust_slice_contrast = true;
reference_size = [800 1140];








%% LOAD IMAGE

% find or create folder location for processed images
folder_processed_images = [save_folder 'processed\\'];
if ~exist(folder_processed_images)
    mkdir(folder_processed_images)
end
warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');
slice_num = 1; close all

% loop across histology images
for file_num = 1:length(image_file_names)
    % load histology image
    disp('loading image...')
    try; figure(histology_figure);
    catch; histology_figure = figure('Name','Histology Viewer'); end
    
    % load already processed image
    if use_already_downsampled_image
        image = imread([image_folder image_file_names{file_num}(1:end-4) '_processed.tif']);
    else %process image now
        image = imread([image_folder image_file_names{file_num}]); 
        original_image_size = size(image);

        % resize (downsample) image to 25 micron pixels
        image = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
        original_image = image*gain;

    


%% ADJUST CONTRAST IF DESIRED

        % display image and set up user controls for contrast change        
        ud_histology.contrast = [0 1]; ud_histology.break = 0; ud_histology.key = 0; ud_histology.show_original = 0; ud_histology.contrast_type = 2;
        set(histology_figure, 'UserData', ud_histology);
        set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyHotkeyFcn(histology_figure, keydata));
        set(histology_figure, 'WindowScrollWheelFcn', @(src,evt)HistologyScrollFcn(histology_figure, evt))

        fprintf(1, '\n Controls: \n \n');
        fprintf(1, 'scroll: adjust contrast \n');
        fprintf(1, 'space: switch btwn adjusting upper and lower saturation points \n');
        fprintf(1, 'e: view original version \n');
        fprintf(1, 'any other key: return to modified version \n');
        fprintf(1, 'r: reset to original \n');
        fprintf(1, 'q: move to next channel \n');

        % adjust contrast in each channel
        if adjust_histology_contrast
            for channel = [3,2,1]
                image_channel = image(:,:,channel);
                imshow(original_image)
                while 1   
                    ud_histology = get(histology_figure, 'UserData');    
                    if ud_histology.key
                        if strcmp(ud_histology.key,'e'); imshow(original_image)
                        else
                            adjusted_image = imadjust(original_image(:,:,channel),ud_histology.contrast);
                            image(:,:,channel) = adjusted_image; imshow(image)
                        end
                    end
                    if ud_histology.break
                        ud_histology.break = 0; ud_histology.contrast = [0 1]; set(histology_figure, 'UserData', ud_histology); break
                    else
                        ud_histology.key = 0; set(histology_figure, 'UserData', ud_histology); pause(.02)
                    end
                end
            end
        end
        imwrite(image, [save_folder image_file_names{file_num}(1:end-4) '_processed.tif'])
    end
    imshow(image)

%% CROP, PROCESS, AND SAVE SLICES


    % display image and set up user controls for rotation
    try; figure(slice_figure);
    catch; slice_figure = figure('Name','Slice Viewer'); end
    ud_slice.rotate_angle = 0; ud_slice.break = 0; ud_slice.key = 0; ud_slice.grid = zeros(1); ud_slice.size = zeros(3,1);
    set(slice_figure, 'UserData', ud_slice);
    set(slice_figure, 'KeyPressFcn', @(slice_figure,keydata)SliceHotkeyFcn(slice_figure, keydata));
    set(slice_figure, 'WindowScrollWheelFcn', @(src,evt)SliceScrollFcn(slice_figure, evt))

    fprintf(1, '\n Controls: \n \n');
    fprintf(1, 'scroll: adjust angle \n');
    fprintf(1, 'f: flip horizontally \n');
    fprintf(1, 'c: crop image \n');
    fprintf(1, 'r: return to original \n');
    fprintf(1, 'g: add gridlines \n');
    fprintf(1, 's: sharpen image \n');
    fprintf(1, 'q: move to next slice \n \n');
    fprintf(1, 'space or close histology viewer: move to next histology image \n \n');
    
    % loop over slices
    while 1
        % crop particular slice (anterior to posterior)
        try; figure(histology_figure);
        catch; histology_figure = figure('Name','Histology Viewer'); end
        imshow(image)
        try
        disp('please select an ROI')
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        slice_image = image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
        catch; break; end
        
        original_slice_image = slice_image; crop_flip_slice_image = slice_image; 
        ud_slice.grid = zeros(size(slice_image),'uint8'); ud_slice.size = size(slice_image); set(slice_figure, 'UserData', ud_slice);
        

    
        % rotate/flip/crop image
        figure(slice_figure)
        imshow(slice_image)
        while 1
            ud_slice = get(slice_figure, 'UserData');   
            
            if ud_slice.key
                % flip image
                if strcmp(ud_slice.key,'f')
                    crop_flip_slice_image = flip(slice_image,2);
                    slice_image = crop_flip_slice_image;
                % rotate image
                elseif strcmp(ud_slice.key,'scroll')
                    slice_image = imrotate(crop_flip_slice_image,ud_slice.rotate_angle,'nearest','crop');
                % crop image
                elseif strcmp(ud_slice.key,'c') 
                    cropped_slice_rect = imrect;
                    slice_position = cropped_slice_rect.getPosition;
                    try
                        crop_flip_slice_image = slice_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);          
                        slice_image = crop_flip_slice_image; 
                        ud_slice.size = size(slice_image); 
                        ud_slice.grid = imresize(ud_slice.grid, ud_slice.size(1:2)); 
                    catch; disp('crop out of bounds'); end
                % enchance contrast
                elseif strcmp(ud_slice.key,'s')
                    slice_image = localcontrast(slice_image); crop_flip_slice_image = slice_image; 
                %reset image
                elseif strcmp(ud_slice.key,'r')
                    slice_image = original_slice_image; crop_flip_slice_image = original_slice_image;
                    ud_slice.size = size(slice_image); 
                    ud_slice.grid = imresize(ud_slice.grid, ud_slice.size(1:2)); 
                end
                % display image
                imshow(slice_image + uint16(ud_slice.grid))
            end
            
            % finished with this slice
            if ud_slice.break
                ud_slice.rotate_angle = 0; ud_slice.break = 0; set(slice_figure, 'UserData', ud_slice); break
            else
                ud_slice.key = 0; set(slice_figure, 'UserData', ud_slice); pause(.02)
            end
        end
    % pad and save slice image
    try; slice_image = padarray(slice_image, [floor((reference_size(1) - size(slice_image,1)) / 2) + mod(size(slice_image,1),2) ...
                                                      floor((reference_size(2) - size(slice_image,2)) / 2) + mod(size(slice_image,2),2)],0);
    imwrite(slice_image, [folder_processed_images save_file_name num2str(file_num) '_' num2str(slice_num) '.tif'])
    slice_num=slice_num+1;                                                  
    catch; disp('saving failed; image must be under reference brain image size'); end                                                 

    
    % finished with this entire histology image
    if strcmp(ud_slice.key,'space')
        break
    end
    
    end   
end


%% GO THROUGH AND FLIP HORIZONTAL SLICE ORIENTATION
processed_images = dir([folder_processed_images '*tif']);
total_num_files = size(processed_images,1); disp(['found ' num2str(total_num_files) ' processed slice images']);
slice_num = 1;

try; figure(flip_figure);
catch; flip_figure = figure('Name','Flip Viewer'); end
ud_flip.break = 0; ud_flip.flip = 0; ud_flip.key = 1; 
set(flip_figure, 'UserData', ud_flip);
set(flip_figure, 'KeyPressFcn', @(flip_figure, keydata)FlipHotkeyFcn(flip_figure, keydata));
fprintf(1, '\n Controls: \n \n');
fprintf(1, 'left/right arrow: go to next slice \n');
fprintf(1, 'f: flip horizontally \n');
fprintf(1, 's: switch order of current and following slice \n');

while 1
    ud_flip = get(flip_figure, 'UserData');   
    
    if ud_flip.key
        % load image
        if ud_flip.key==1
            processed_image_name = processed_images(slice_num).name;
            current_slice_image = imread([folder_processed_images processed_image_name]);
        elseif ud_flip.key=='f'
            current_slice_image = flip(current_slice_image,2);
        elseif ud_flip.key=='s'
            disp('switching order -- moving this image forward')
            if slice_num < total_num_files
                next_processed_image_name = processed_images(slice_num+1).name;
                next_slice_image = imread([folder_processed_images next_processed_image_name]);

                imwrite(next_slice_image, [folder_processed_images processed_image_name])            
                imwrite(current_slice_image, [folder_processed_images next_processed_image_name])
            end
        end
        % show image
        imshow(current_slice_image)
    end
    
    % save flipped image, and go to next image or break
    if ud_flip.break
        if ud_flip.flip
            imwrite(current_slice_image, [folder_processed_images processed_image_name])
        end
        if strcmp(ud_flip.key,'q')
            disp('done'), break
        else
            if strcmp(ud_flip.key,'leftarrow') && slice_num > 1
                direction = -1;
            elseif strcmp(ud_flip.key,'rightarrow') && slice_num < total_num_files
                direction = 1;
            else; direction = 0;
            end
            slice_num = slice_num + direction;
            ud_flip.flip = 0; ud_flip.break = 0; ud_flip.key = 1; set(flip_figure, 'UserData', ud_flip); 
         end
    else
        ud_flip.key = 0; set(flip_figure, 'UserData', ud_flip); pause(.02)
    end 
end




