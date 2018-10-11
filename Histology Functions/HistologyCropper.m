function HistologyCropper(histology_figure, save_folder, image_file_names, reference_size, save_file_name, use_already_downsampled_image)

% set up histology figure
ud_histology.file_num = 1;
ud_histology.num_files = length(image_file_names);
ud_histology.slice_num = ones(length(image_file_names),1);
ud_histology.save_file_name =save_file_name;

if use_already_downsampled_image
    ud_histology.file_name_suffix = '';
else
    ud_histology.file_name_suffix = '_processed';
end
ud_histology.hist_image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) ud_histology.file_name_suffix '.tif']));
figure(histology_figure); 
ud_histology.im = imshow(ud_histology.hist_image);
warning('off', 'MATLAB:colon:nonIntegerIndex');

ud_histology.cropped_slice_rect = {};

set(histology_figure, 'UserData', ud_histology);

% crop and switch image function for histology
set(histology_figure, 'KeyPressFcn', @(histology_figure,keydata)HistologyCropHotkeyFcn(histology_figure, keydata, save_folder, image_file_names, reference_size));

fprintf(1, '\n Controls: \n \n');
fprintf(1, 'click and drag: crop image \n');
fprintf(1, 'space: next image \n \n');

disp('select ROIs')
ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);

end

% --------------------------------------------------------------------------
% use imrect to crop, and then save the cropped image in 'processed' folder
% --------------------------------------------------------------------------
function ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size)

    try % get first slice ROI
        ud_histology.cropped_slice_rect{end+1} = imrect;
        slice_position = ud_histology.cropped_slice_rect{end}.getPosition;
        ud_histology.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);

        % pad if possible (if small enough)
        try; ud_histology.slice_image = padarray(ud_histology.slice_image, [floor((reference_size(1) - size(ud_histology.slice_image,1)) / 2) + mod(size(ud_histology.slice_image,1),2) ...
                                                      floor((reference_size(2) - size(ud_histology.slice_image,2)) / 2) + mod(size(ud_histology.slice_image,2),2)],0);
        end         

        % save cropped slice
        imwrite(ud_histology.slice_image, fullfile(save_folder, 'processed', ...
                [ud_histology.save_file_name num2str(ud_histology.file_num,'%.2d') '.' num2str(ud_histology.slice_num(ud_histology.file_num),'%.3d') '.tif']))
        disp([ud_histology.save_file_name num2str(ud_histology.file_num,'%.2d') '.' num2str(ud_histology.slice_num(ud_histology.file_num),'%.3d') ' saved!'])

        ud_histology.slice_num(ud_histology.file_num) = ud_histology.slice_num(ud_histology.file_num) + 1;

    catch; ud_histology.slice_image = zeros(100,100,'int16');
    end
    
    set(histology_figure, 'UserData', ud_histology);
    
    try
    ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);
    catch; disp('')
    end
end



% -------------------------------
% respond to keypress (space bar)
% -------------------------------
function HistologyCropHotkeyFcn(histology_figure, keydata, save_folder, image_file_names, reference_size)

ud_histology = get(histology_figure, 'UserData');

    switch lower(keydata.Key)    
        
    case 'space' % move onto next image
        if ud_histology.file_num + 1 <= ud_histology.num_files
            ud_histology.file_num = ud_histology.file_num + 1;

            ud_histology.hist_image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) ud_histology.file_name_suffix '.tif']));
            figure(histology_figure); 
            set(ud_histology.im, 'CData', ud_histology.hist_image); 

            for i = 1:length(ud_histology.cropped_slice_rect)
            delete(ud_histology.cropped_slice_rect{i})
            end
            ud_histology.cropped_slice_rect = {};
            
            ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);

        else
           disp('That was the last file -- close and move on to the next cell') 
        end
    end

try   
set(histology_figure, 'UserData', ud_histology);
catch; disp('')
end


end
