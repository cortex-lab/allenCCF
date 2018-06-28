function HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, file_num, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)

% display image and set up user controls for contrast change        
ud_histology.contrast = [0 1]; 
ud_histology.break = 0; 
ud_histology.key = 0; 
ud_histology.show_original = 0; 
ud_histology.contrast_type = 2;
ud_histology.channel = 3;
ud_histology.file_num = file_num;
ud_histology.num_files = length(image_file_names);
ud_histology.save_folder = save_folder;
ud_histology.image_folder = image_folder;
ud_histology.microns_per_pixel = microns_per_pixel;
ud_histology.microns_per_pixel_after_downsampling = microns_per_pixel_after_downsampling;
ud_histology.gain = gain;

% load histology image
disp(['loading image ' num2str(file_num) '...'])


% load already processed image
if use_already_downsampled_image
    image = imread([save_folder image_file_names{file_num}(1:end-4) '_processed.tif']);
else %process image now
    image = imread([image_folder image_file_names{file_num}]);
    original_image_size = size(image);

    % resize (downsample) image to 25 micron pixels
    image = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
end
original_image = image*gain;
imshow(original_image);




ud_histology.original_image = original_image;
ud_histology.adjusted_image = original_image;



set(histology_figure, 'UserData', ud_histology);

set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyHotkeyFcn(histology_figure, keydata, image_file_names, use_already_downsampled_image));
set(histology_figure, 'WindowScrollWheelFcn', @(src,evt)HistologyScrollFcn(histology_figure, evt))

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'scroll: adjust contrast \n');
fprintf(1, 'space: switch btwn adjusting upper and lower saturation points \n');
fprintf(1, 'e: view original version \n');
fprintf(1, 'any other key: return to modified version \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'c: move to next channel \n');
fprintf(1, 's: save image \n');
fprintf(1, 'left/right arrow: save and move to next slide image \n');




function HistologyHotkeyFcn(fig, keydata, image_file_names, use_already_downsampled_image)

ud = get(fig, 'UserData');

switch lower(keydata.Key)    
    case 'e' % show original
        figure(fig);
        imshow(ud.original_image)
    case 'r' % return to original
        ud.contrast = [0 1];       
    case 'space'
        disp('switch contrast effect')
        if ud.contrast_type==2
            ud.contrast_type = 1;
        elseif ud.contrast_type==1
            ud.contrast_type = 2;
    end
    case 'c' % break
        disp('next channel')
        ud.channel = ud.channel + 1 - (ud.channel==3)*3;
    case 's'
        disp('saving downsampled and processed image');
        imwrite(ud.adjusted_image, [ud.save_folder image_file_names{ud.file_num}(1:end-4) '_processed.tif'])
        imshow(ud.adjusted_image)
    case 'leftarrow'
    disp('saving downsampled and processed image');
    imwrite(ud.adjusted_image, [ud.save_folder image_file_names{ud.file_num}(1:end-4) '_processed.tif'])
    imshow(ud.adjusted_image)           
        if ud.file_num > 1
            ud.file_num = ud.file_num - 1;
            move_on = true;
        else
            move_on = false;
        end
        case 'rightarrow'
    disp('saving downsampled and processed image');
    imwrite(ud.adjusted_image, [ud.save_folder image_file_names{ud.file_num}(1:end-4) '_processed.tif'])
    imshow(ud.adjusted_image)               
        if ud.file_num < ud.num_files;
            ud.file_num = ud.file_num + 1;
            move_on = true;
        else
            disp('that''s all, folks; save this image and continue to the next cell')
            move_on = false;
        end
end

if (strcmp(lower(keydata.Key),'leftarrow') || strcmp(lower(keydata.Key),'rightarrow')) && move_on
            
    % load histology image
    disp(['loading image ' num2str(ud.file_num) '...'])
    
    % load already processed image
    if use_already_downsampled_image
        image = imread([ud.save_folder image_file_names{ud.file_num}(1:end-4) '_processed.tif']);
    else %process image now
        image = imread([ud.image_folder image_file_names{ud.file_num}]);
        original_image_size = size(image);

        % resize (downsample) image to 25 micron pixels
        image = imresize(image, [round(original_image_size(1)*ud.microns_per_pixel/ud.microns_per_pixel_after_downsampling)  NaN]);
    end
    original_image = image*ud.gain;
    imshow(original_image);

    ud.original_image = original_image;
    ud.adjusted_image = original_image;


    
end


ud.key = keydata.Key;
set(fig, 'UserData', ud);





function HistologyScrollFcn(fig, evt)

ud = get(fig, 'UserData');
ud.key = 'scroll';


%modify based on scrolling
ud.contrast(ud.contrast_type) = ud.contrast(ud.contrast_type) + evt.VerticalScrollCount*.1;

% make sure within limit of 0 to 1
if ud.contrast(ud.contrast_type) < 0
    ud.contrast(ud.contrast_type) = 0;
elseif ud.contrast(ud.contrast_type) > 1
    ud.contrast(ud.contrast_type) = 1;
end

try
adjusted_image_slice = imadjust(ud.original_image(:,:,ud.channel),ud.contrast);
ud.adjusted_image(:,:,ud.channel) = adjusted_image_slice; 
imshow(ud.adjusted_image)
catch
    disp('parameter out of bounds')
end
    

set(fig, 'UserData', ud);



    
