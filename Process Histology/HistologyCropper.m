function HistologyCropper(histology_figure, slice_figure, save_folder, image_file_names, reference_size, save_file_name)


% set up histology figure
ud_histology.file_num = 1;
ud_histology.slice_num = ones(length(image_file_names),1);

ud_histology.hist_image = imread([save_folder image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']);
figure(histology_figure); imshow(ud_histology.hist_image);

try % get first slice ROI
disp('please select an ROI')
cropped_slice_rect = imrect;
slice_position = cropped_slice_rect.getPosition;
ud.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
ud.slice_image = localcontrast(ud.slice_image);
catch; disp('cropping failed'); end

%set up slice figure
figure(slice_figure);
imshow(ud.slice_image);
ud.original_slice_image = ud.slice_image;

ud.size = size(ud.slice_image);
ud.save_file_name = save_file_name;


set(histology_figure, 'UserData', ud_histology);
set(slice_figure, 'UserData', ud);


% crop function for slice
set(slice_figure, 'KeyPressFcn', @(slice_figure,keydata)SliceCropHotkeyFcn(slice_figure, histology_figure, keydata, save_folder, reference_size));

% crop and switch image function for histology
set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyCropHotkeyFcn(slice_figure, histology_figure, keydata, save_folder, image_file_names));



fprintf(1, '\n Controls: \n \n');
fprintf(1, 'space: next image \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'c: crop slice further \n');






function SliceCropHotkeyFcn(slice_figure, histology_figure, keydata, save_folder, reference_size);

ud = get(slice_figure, 'UserData');
ud_histology = get(histology_figure, 'UserData');

switch lower(keydata.Key)    
    case 'c' % crop
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        try
            crop_flip_slice_image = ud.slice_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);          
            ud.slice_image = crop_flip_slice_image; 
            ud.size = size(ud.slice_image); 
            imshow(ud.slice_image)
        catch; disp('crop out of bounds');
        end
        
    case 'r' % return to original
        ud.slice_image = ud.original_slice_image;
        imshow(ud.slice_image)
        
    case 'space' % move onto next image
    % pad and save slice image
        try; ud.slice_image = padarray(ud.slice_image, [floor((reference_size(1) - size(ud.slice_image,1)) / 2) + mod(size(ud.slice_image,1),2) ...
                                                          floor((reference_size(2) - size(ud.slice_image,2)) / 2) + mod(size(ud.slice_image,2),2)],0);
        catch; disp('image must be under reference brain image size -- make sure to crop in next stage of preprocessing!');
        end              
        imwrite(ud.slice_image, [save_folder 'processed\\' ...
                    ud.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) '.tif'])
        
        ud_histology.slice_num(ud_histology.file_num) = ud_histology.slice_num(ud_histology.file_num) + 1;        
        disp([ud.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) ' saved!'])
        
        figure(histology_figure);
        
      
        try
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        ud.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
        ud.slice_image = localcontrast(ud.slice_image);
        ud.original_slice_image = ud.slice_image;
        figure(slice_figure);
        imshow(ud.slice_image);
        ud.original_slice_image = ud.slice_image;

        ud.size = size(ud.slice_image);          
        catch; disp('cropping failed'); end 
        

end


set(histology_figure, 'UserData', ud_histology);
set(slice_figure, 'UserData', ud);





function HistologyCropHotkeyFcn(slice_figure, histology_figure, keydata, save_folder, image_file_names)

ud = get(slice_figure, 'UserData');
ud_histology = get(histology_figure, 'UserData');

switch lower(keydata.Key)    
        
    case 'space' % move onto next image
        ud_histology.file_num = ud_histology.file_num + 1;
        
        ud_histology.hist_image = imread([save_folder image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']);
        figure(histology_figure); imshow(ud_histology.hist_image);
        
        try % get first slice ROI
        disp('please select an ROI')
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        ud.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
        ud.slice_image = localcontrast(ud.slice_image);
        ud.original_slice_image = ud.slice_image;
        catch; disp('cropping failed'); end
end


set(slice_figure, 'UserData', ud_slice);
set(histology_figure, 'UserData', ud_histology);



