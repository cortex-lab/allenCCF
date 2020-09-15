function HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, folder_processed_images, image_file_are_individual_slices, ...
                use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)

% display image and set up user controls for contrast change        
ud.show_original = 0; 
ud.adjusting_contrast = 0;
ud.file_num = 1;
ud.num_files = length(image_file_names);
if image_file_are_individual_slices
    ud.save_folder = folder_processed_images;
else
    ud.save_folder = save_folder;
end
ud.image_folder = image_folder;
ud.microns_per_pixel = microns_per_pixel;
ud.microns_per_pixel_after_downsampling = microns_per_pixel_after_downsampling;
ud.gain = gain;

% load histology image
disp(['loading image ' num2str(ud.file_num) '...'])

% load first image
image = imread(fullfile(ud.image_folder,image_file_names{ud.file_num}));


if ~use_already_downsampled_image
    % resize (downsample) image to reference atlas size
    disp('downsampling image...')
    original_image_size = size(image);
    image = imresize(image, [round(original_image_size(1)*microns_per_pixel/microns_per_pixel_after_downsampling)  NaN]);
end
ud.file_name_suffix = '_processed';
ud.channel = min( 3, size(image,3));
original_image = image(:,:,1:ud.channel)*gain;

imshow(original_image);
title(['Adjusting channel ' num2str(ud.channel) ' on image ' num2str(ud.file_num) ' / ' num2str(ud.num_files)],...
                    'color',[1==ud.channel 2==ud.channel 3==ud.channel])
                
ud.original_image = original_image;
ud.adjusted_image = original_image;

imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) ud.file_name_suffix '.tif']))

set(histology_figure, 'UserData', ud);

set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyHotkeyFcn(histology_figure, keydata, image_file_names, use_already_downsampled_image));

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'space: adjust contrast for current channel / return to image-viewing mode \n');
fprintf(1, 'e: view original version \n');
fprintf(1, 'any key: return to modified version \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'c: move to next channel \n');
fprintf(1, 's: save image \n');
fprintf(1, 'left/right arrow: save and move to next slide image \n');




% --------------------
% Respond to keypress
% --------------------
function HistologyHotkeyFcn(fig, keydata, image_file_names, use_already_downsampled_image)

ud = get(fig, 'UserData');

if strcmp(lower(keydata.Key), 'space') % adjust contrast
    ud.adjusting_contrast = ~ud.adjusting_contrast;

    if ud.adjusting_contrast
        disp(['adjust contrast on channel ' num2str(ud.channel)])
        imshow(ud.adjusted_image(:,:,ud.channel))
        imcontrast(fig)
    else
        adjusted_image_channel = fig.Children.Children.CData;
        ud.adjusted_image(:,:,ud.channel) = adjusted_image_channel;
    end   

% ignore commands while adjusting contrast    
elseif ~ud.adjusting_contrast     
    switch lower(keydata.Key)    
        case 'e' % show original
            ud.show_original = ~ud.show_original;
            if ud.show_original 
                disp('showing original image (press any key to return)')
                imshow(ud.original_image)
            end    
        case 'r' % return to original
            disp('revert to original image')
            ud.adjusted_image = ud.original_image;    
        case 'c' % break
            disp('next channel')
            ud.channel = ud.channel + 1 - (ud.channel==3)*3;

        case 's' % save image
            disp(['saving processed image ' num2str(ud.file_num)]);
            imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) ud.file_name_suffix '.tif']))
            imshow(ud.adjusted_image)
        case 'leftarrow' % save image and move to previous image
            disp(['saving processed image ' num2str(ud.file_num)]);
            imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) ud.file_name_suffix '.tif']))
        
            if ud.file_num > 1
                ud.file_num = ud.file_num - 1;
                move_on = true;
            else
                move_on = false;
            end
        case 'rightarrow' % save image and move to next image
            disp(['saving processed image ' num2str(ud.file_num)]);
            imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) ud.file_name_suffix '.tif']))
             
            if ud.file_num < ud.num_files;
                ud.file_num = ud.file_num + 1;
                move_on = true;        
            else
                disp('that''s all, folks; continue to the next cell')
                move_on = false;
            end
    end
    if (strcmp(lower(keydata.Key),'leftarrow') || strcmp(lower(keydata.Key),'rightarrow')) && move_on

        % load image
        image = imread(fullfile(ud.image_folder, image_file_names{ud.file_num}) );
        disp(['image ' num2str(ud.file_num) ' loaded'])
        
        if ~use_already_downsampled_image
            % resize (downsample) image to reference size
            disp('downsampling image...')
            original_image_size = size(image);
            image = imresize(image, [round(original_image_size(1)*ud.microns_per_pixel/ud.microns_per_pixel_after_downsampling)  NaN]);
        end
        original_image = image*ud.gain;

        ud.original_image = original_image;
        ud.adjusted_image = original_image;

        % save immediately
        imwrite(ud.adjusted_image, fullfile(ud.save_folder, [image_file_names{ud.file_num}(1:end-4) ud.file_name_suffix '.tif']))

    end
else % if pressing commands while adjusting contrast
    disp(' ')
    disp('Please press space to exit contrast adjustment before issuing other commands')
    disp('If you are dissatisfied with your changes, you can then press ''r'' to revert to the original image')
end



% show the image, unless in other viewing modes
figure(fig)
if ~(ud.adjusting_contrast || (strcmp(lower(keydata.Key),'e')&&ud.show_original) )
    imshow(ud.adjusted_image)
end
title(['Adjusting channel ' num2str(ud.channel) ' on image ' num2str(ud.file_num) ' / ' num2str(ud.num_files)],...
            'color',[1==ud.channel 2==ud.channel 3==ud.channel])

set(fig, 'UserData', ud);


