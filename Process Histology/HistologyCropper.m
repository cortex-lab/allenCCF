function HistologyCropper(histology_figure, slice_figure, save_folder, image_file_names, reference_size, save_file_name)

% set up histology figure
ud_histology.file_num = 1;
ud_histology.num_files = length(image_file_names);
ud_histology.slice_num = ones(length(image_file_names),1);
ud_histology.save_file_name =save_file_name;

ud_histology.hist_image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']));
figure(histology_figure); 
ud_histology.im = imshow(ud_histology.hist_image);

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

        % save cropped slice
        try; ud_histology.slice_image = padarray(ud_histology.slice_image, [floor((reference_size(1) - size(ud_histology.slice_image,1)) / 2) + mod(size(ud_histology.slice_image,1),2) ...
                                                      floor((reference_size(2) - size(ud_histology.slice_image,2)) / 2) + mod(size(ud_histology.slice_image,2),2)],0);
        catch; disp(''); disp('image must be under reference brain image size -- make sure to crop in next stage of preprocessing!');
        end         

        imwrite(ud_histology.slice_image, fullfile(save_folder, 'processed', ...
                [ud_histology.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) '.tif']))
        disp([ud_histology.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) ' saved!'])

        ud_histology.slice_num(ud_histology.file_num) = ud_histology.slice_num(ud_histology.file_num) + 1;

    catch; disp('cropping failed'); ud_histology.slice_image = zeros(100,100,'int16');
    end
    
    set(histology_figure, 'UserData', ud_histology);
    
    ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);
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

            ud_histology.hist_image = imread(fullfile(save_folder, [image_file_names{ud_histology.file_num}(1:end-4) '_processed.tif']));
            figure(histology_figure); 
            set(ud_histology.im, 'CData', ud_histology.hist_image); 

            for i = 1:length(ud_histology.cropped_slice_rect)
            delete(ud_histology.cropped_slice_rect{i})
            end
            ud_histology.cropped_slice_rect = {};
            
            ud_histology = crop_and_save_image(ud_histology, histology_figure, save_folder, reference_size);
%             
%             try % get first slice ROI
%             disp('please select an ROI')
%             ud_histology.cropped_slice_rect{end+1} = imrect;
%             slice_position = ud_histology.cropped_slice_rect{end}.getPosition;
%             ud_histology.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
%             
%             % save cropped slice
%             try; ud_histology.slice_image = padarray(ud_histology.slice_image, [floor((reference_size(1) - size(ud_histology.slice_image,1)) / 2) + mod(size(ud_histology.slice_image,1),2) ...
%                                                           floor((reference_size(2) - size(ud_histology.slice_image,2)) / 2) + mod(size(ud_histology.slice_image,2),2)],0);
%             catch; disp(''); disp('image must be under reference brain image size -- make sure to crop in next stage of preprocessing!');
%             end         
%             
%             imwrite(ud_histology.slice_image, fullfile(save_folder, 'processed', ...
%                     [ud_histology.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) '.tif']))
%             disp([ud_histology.save_file_name num2str(ud_histology.file_num) '.' num2str(ud_histology.slice_num(ud_histology.file_num)) ' saved!'])
%             
%             ud_histology.slice_num(ud_histology.file_num) = ud_histology.slice_num(ud_histology.file_num) + 1;    
%         
%             figure(histology_figure);
%             try
%             ud_histology.cropped_slice_rect{end+1} = imrect;
%             slice_position = ud_histology.cropped_slice_rect{end}.getPosition;
%             ud_histology.slice_image = ud_histology.hist_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);
%             ud_histology.slice_image = localcontrast(ud_histology.slice_image);
%             ud_histology.original_slice_image = ud_histology.slice_image;
%             figure(slice_figure);
%             imshow(ud_histology.slice_image);
% 
%             ud_histology.size = size(ud_histology.slice_image);    
% 
%             set(histology_figure, 'UserData', ud_histology);       
%             catch; 
%                 disp('');
%             end 
%         
%             catch; 
%                 disp('cropping failed'); 
%             end
        else
           disp('That was the last file -- close and move on to the next cell') 
        end
    end

set(histology_figure, 'UserData', ud_histology);

end
