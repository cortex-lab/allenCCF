function HistologyBrowser(histology_figure, original_image, save_folder, image_file_names, file_num)

% display image and set up user controls for contrast change        
ud_histology.contrast = [0 1]; 
ud_histology.break = 0; 
ud_histology.key = 0; 
ud_histology.show_original = 0; 
ud_histology.contrast_type = 2;
ud_histology.channel = 3;

ud_histology.original_image = original_image;
ud_histology.adjusted_image = original_image;

figure(histology_figure); imshow(original_image);

set(histology_figure, 'UserData', ud_histology);

set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyHotkeyFcn(histology_figure, keydata, save_folder, image_file_names, file_num));
set(histology_figure, 'WindowScrollWheelFcn', @(src,evt)HistologyScrollFcn(histology_figure, evt))

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'scroll: adjust contrast \n');
fprintf(1, 'space: switch btwn adjusting upper and lower saturation points \n');
fprintf(1, 'e: view original version \n');
fprintf(1, 'any other key: return to modified version \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'c: move to next channel \n');
fprintf(1, 's: save image \n');




function HistologyHotkeyFcn(fig, keydata, save_folder, image_file_names, file_num)

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
        imwrite(ud.adjusted_image, [save_folder image_file_names{file_num}(1:end-4) '_processed.tif'])
        imshow(ud.adjusted_image)
end

ud.key = keydata.Key;
set(fig, 'UserData', ud);





function HistologyScrollFcn(fig, evt)

ud = get(fig, 'UserData');
ud.key = 'scroll';


%modify based on scrolling
ud.contrast(ud.contrast_type) = ud.contrast(ud.contrast_type) + evt.VerticalScrollCount*.05;

% make sure within limit of 0 to 1
if ud.contrast(ud.contrast_type) < 0
    ud.contrast(ud.contrast_type) = 0;
elseif ud.contrast(ud.contrast_type) > 1
    ud.contrast(ud.contrast_type) = 1;
end


adjusted_image_slice = imadjust(ud.original_image(:,:,ud.channel),ud.contrast);
ud.adjusted_image(:,:,ud.channel) = adjusted_image_slice; 
imshow(ud.adjusted_image)


set(fig, 'UserData', ud);



    
