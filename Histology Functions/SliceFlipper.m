function SliceFlipper(slice_figure, folder_processed_images, reference_size)
% crop, sharpen, and flip slice images

processed_images = dir([folder_processed_images filesep '*tif']);
ud.processed_image_names = natsortfiles({processed_images.name});
ud.total_num_files = size(processed_images,1); disp(['found ' num2str(ud.total_num_files) ' processed slice images']);

ud.slice_num = 1;
ud.flip = 0; 
ud.rotate_angle = 0;
ud.reference_size = reference_size;

ud.grid = 0;
ud = load_next_slice(ud,folder_processed_images);

ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image)); 
ud.grid(1:50:end,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16')); 
ud.grid(:,1:50:end,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));    

imshow(ud.current_slice_image + ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])
set(slice_figure, 'UserData', ud);

     

% key function for slice
set(slice_figure, 'KeyPressFcn', @(slice_figure,keydata)SliceCropHotkeyFcn(keydata, slice_figure, folder_processed_images));
% scroll function for slice
set(slice_figure, 'WindowScrollWheelFcn', @(src,evt)SliceScrollFcn(slice_figure, evt))


fprintf(1, '\n Controls: \n \n');
fprintf(1, 'right: save and see next image \n');
fprintf(1, 'left: save and see previous image \n');
fprintf(1, 'scroll: rotate slice \n');
fprintf(1, 's: sharpen \n');
fprintf(1, 'g: toggle grid \n');
fprintf(1, 'c: crop slice further \n');
fprintf(1, 'f: flip horizontally \n');
fprintf(1, 'w: switch order (move image forward) \n');
fprintf(1, 'r: reset to original \n');
fprintf(1, 'delete: delete current image \n');

% --------------------
% respond to keypress
% --------------------
function SliceCropHotkeyFcn(keydata, slice_figure, folder_processed_images)

ud = get(slice_figure, 'UserData');



switch lower(keydata.Key)   
    case 'leftarrow' % save and previous slice
        imwrite(ud.current_slice_image, fullfile(folder_processed_images, ud.processed_image_name))        
        ud.slice_num = ud.slice_num - 1*(ud.slice_num>1);
        ud = load_next_slice(ud,folder_processed_images);

    case 'rightarrow' % save and next slice      
        imwrite(ud.current_slice_image, fullfile(folder_processed_images, ud.processed_image_name))
        ud.slice_num = ud.slice_num + 1*(ud.slice_num < length(ud.processed_image_names));
        ud = load_next_slice(ud,folder_processed_images);
        
    case 'delete' % delete and previous slice
        delete(fullfile(folder_processed_images, ud.processed_image_name));
        ud.slice_num = ud.slice_num - 1*(ud.slice_num>1);  
        
        processed_images = dir([folder_processed_images filesep '*tif']);
        ud.processed_image_names = natsortfiles({processed_images.name});
        ud.total_num_files = size(processed_images,1); disp(['found ' num2str(ud.total_num_files) ' processed slice images']);
        ud = load_next_slice(ud,folder_processed_images);

    case 'g' % grid
        if sum(ud.grid(:)) == 0
            ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image)); 
            ud.grid(1:50:end,:,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16')); 
            ud.grid(:,1:50:end,:) = 150 + 20000*(isa(ud.original_slice_image,'uint16'));             
        else
            ud.grid = zeros(size(ud.current_slice_image),class(ud.original_slice_image)); 
        end
    case 'c' % crop
        cropped_slice_rect = imrect;
        slice_position = cropped_slice_rect.getPosition;
        try
            ud.current_slice_image = ud.current_slice_image(slice_position(2):slice_position(2)+slice_position(4),slice_position(1):slice_position(1)+slice_position(3),:);          
        catch; disp('crop out of bounds'); 
        end        
        
        try; ud.current_slice_image = padarray(ud.current_slice_image, [floor((ud.reference_size(1) - size(ud.current_slice_image,1)) / 2) + ...
                                mod(size(ud.current_slice_image,1),2) floor((ud.reference_size(2) - size(ud.current_slice_image,2)) / 2) + ...
                                mod(size(ud.current_slice_image,2),2)],0);
            ud.original_ish_slice_image = ud.current_slice_image;                            
        catch; disp('cropping failed');
        end              
        
        ud.size = size(ud.current_slice_image); 
        ud.grid = imresize(ud.grid, ud.size(1:2));         
        
	case 's' % sharpen
        ud.current_slice_image = localcontrast(ud.current_slice_image);
        ud.original_ish_slice_image = localcontrast(ud.original_ish_slice_image);
        
    case 'w' % switch order
        if ud.slice_num < length(ud.processed_image_names)
            disp('switching order -- moving this image forward')
            next_processed_image_name = ud.processed_image_names{ud.slice_num+1};
            next_slice_image = imread(fullfile(folder_processed_images, next_processed_image_name));

            imwrite(next_slice_image, fullfile(folder_processed_images, ud.processed_image_name))            
            imwrite(ud.current_slice_image, fullfile(folder_processed_images, next_processed_image_name))
            
            ud.current_slice_image = next_slice_image; 
            ud.size = size(ud.current_slice_image); 
            ud.grid = imresize(ud.grid, ud.size(1:2));             
        end
        
    case 'f' % flip horizontally
        ud.current_slice_image = flip(ud.current_slice_image,2);
        ud.original_ish_slice_image = flip(ud.original_ish_slice_image,2);
        
    case 'r' % return to original image
        ud.current_slice_image = ud.original_slice_image;
        ud.original_ish_slice_image = ud.original_slice_image;
        ud.size = size(ud.current_slice_image); 
        ud.grid = imresize(ud.grid, ud.size(1:2)); 
        ud.rotate_angle = 0;

end



% in all cases, update image and title
imshow(ud.current_slice_image+ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])


set(slice_figure, 'UserData', ud);

function ud = load_next_slice(ud,folder_processed_images)
    ud.processed_image_name = ud.processed_image_names{ud.slice_num};
    ud.current_slice_image = imread(fullfile(folder_processed_images, ud.processed_image_name));

    % pad if possible (if small enough)
    try; ud.current_slice_image = padarray(ud.current_slice_image, [floor((ud.reference_size(1) - size(ud.current_slice_image,1)) / 2) + mod(size(ud.current_slice_image,1),2) ...
                                                  floor((ud.reference_size(2) - size(ud.current_slice_image,2)) / 2) + mod(size(ud.current_slice_image,2),2)],0);
    end          

    ud.original_slice_image = ud.current_slice_image;             
    ud.original_ish_slice_image = ud.current_slice_image;        

    ud.size = size(ud.current_slice_image); 
    if ud.size(1) > ud.reference_size(1)+1 || ud.size(2) > ud.reference_size(2)+2
        disp(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files) ' is ' num2str(ud.size(1)) 'x' num2str(ud.size(2)) ' pixels:']);
        disp(['I suggest you crop this image down to under ' num2str(ud.reference_size(1)) ' x ' num2str(ud.reference_size(2)) ' pxl'])
    end        
    ud.grid = imresize(ud.grid, ud.size(1:2)); 
    ud.rotate_angle = 0;
    
    imwrite(ud.current_slice_image, fullfile(folder_processed_images, ud.processed_image_name)) 


% function to rotate slice by scrolling
function SliceScrollFcn(fig, evt)

ud = get(fig, 'UserData');

%modify based on scrolling
ud.rotate_angle = ud.rotate_angle + evt.VerticalScrollCount*.75;

ud.current_slice_image = imrotate(ud.original_ish_slice_image,ud.rotate_angle,'nearest','crop');
imshow(ud.current_slice_image+ud.grid)
title(['Slice ' num2str(ud.slice_num) ' / ' num2str(ud.total_num_files)])

set(fig, 'UserData', ud);



