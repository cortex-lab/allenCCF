function sliceBrowser(slice_figure, processed_images_folder, f, reference_size)
% -----------------------------------------------------------------
% browser to go through histology along with the reference browser
%-----------------------------------------------------------------

% initialize user data variables held by the figure
processed_images = dir([processed_images_folder filesep '*.tif']);
ud_slice.processed_image_names = natsortfiles({processed_images.name});
total_num_files = size(processed_images,1); disp(['found ' num2str(total_num_files) ' processed slice images']);
ud_slice.total_num_files = total_num_files;
ud_slice.break = 0; 
ud_slice.slice_num = 1;
ud_slice.key = 1; 
ud_slice.pointList = []; 
ud_slice.pointHands = [];
ud_slice.getPoint = 0;
ud_slice.ref_size = reference_size(2:3);

ud_slice.processed_images = processed_images;
ud_slice.processed_images_folder = processed_images_folder;
ud_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
hold(ud_slice.sliceAx, 'on');
set(ud_slice.sliceAx, 'HitTest', 'off');
ud_slice.im = plotTVslice(zeros(ud_slice.ref_size, 'uint8'));

% create functions needed to interact with the figure
set(ud_slice.im, 'ButtonDownFcn', @(slice_figure,k)sliceClickCallback(slice_figure, k));
set(slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(slice_figure, keydata, f));
set(slice_figure, 'UserData', ud_slice)

% adjust figure to user's screen size
try; screen_size = get(0,'ScreenSize'); screen_size = [max(screen_size(3:4)) min(screen_size(3:4))]./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end
set(slice_figure,'Position', [150*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
movegui(slice_figure,'onscreen')

% set up first slice image
ud_slice = updateSliceImage(ud_slice);

% ------------------------------------------------    
% Clicking function to register transform points  
% ------------------------------------------------
function sliceClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');


if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));

    ud.pointList(end+1, :) = [clickX, ud.ref_size(1) - clickY];
    ud.pointHands(end+1) = plot(ud.sliceAx, clickX, clickY, 'ro', 'color', [0 .5 0],'linewidth',2,'markers',4);    
    
     if clickX < 100 && (ud.ref_size(1) - clickY) < 100 % if click in corner, break
        ud.pointList = []; 
        set(ud.pointHands(:), 'Visible', 'off');     
     end
    
end
set(f, 'UserData', ud);

% ------------------------
% react to keyboard press
% ------------------------
function SliceAtlasHotkeyFcn(fig, keydata, f)

ud = get(fig, 'UserData');

% left arrow -- go to previous slice
if strcmp(keydata.Key,'leftarrow')    
    if ud.slice_num > 1
        ud.slice_num = ud.slice_num - 1;
        ud = updateSliceImage(ud);
    end
    
% right arrow -- go to next slice    
elseif strcmp(keydata.Key,'rightarrow') 
    if ud.slice_num < ud.total_num_files
        ud.slice_num = ud.slice_num + 1;
        ud = updateSliceImage(ud);
    end
% d -- delete current transform points
elseif strcmp(keydata.Key,'d') 
%     disp('current transform points deleted')
%     set(ud.pointHands(:), 'Visible', 'off'); 
%     ud.pointList = [];    
    
    % Try to delete only most recent point
    set(ud.pointHands(end), 'Visible', 'off'); 
    ud.pointHands = ud.pointHands(1:end-1); 
    ud.pointList = ud.pointList(1:end-1,:); 
    disp('transform point deleted')
% t -- transform point mode
elseif strcmp(keydata.Key,'t')
    ud.getPoint = ~ud.getPoint;
        if ud.getPoint; disp('transform point mode!'); end
else
% otherwise -- call function to atlas browser       
    figure(f);
    fcn = get(f, 'KeyPressFcn'); 
    fcn(f, keydata);
end


set(fig, 'UserData', ud);


function ud = updateSliceImage(ud)

    title_ending = '';
    
    processed_image_name = ud.processed_image_names{ud.slice_num};
    current_slice_image = flip(imread(fullfile(ud.processed_images_folder, processed_image_name)),1);
    % reduce to 3 channels at most
    color_channels = min( 3, size(current_slice_image,3));
    current_slice_image = current_slice_image(:,:,1:color_channels);
    % reduce to reference atlas size
    if size(current_slice_image,1) > ud.ref_size(1)+2 || size(current_slice_image,2) > ud.ref_size(2)+2
        disp(['shrinking image to reference size ' num2str(ud.ref_size(1)) ' x ' num2str(ud.ref_size(2)) ' pxl'])
        current_slice_image = imresize(current_slice_image, ud.ref_size);
    end          
    set(ud.im, 'CData', current_slice_image); 


    file_transformations = fullfile(ud.processed_images_folder, 'transformations\\' ,...
                            [processed_image_name(1:end-4) '_transform_data.mat']);

    set(ud.pointHands(:), 'Visible', 'off'); 
    ud.pointList = [];
    
    if exist(file_transformations,'file')
        % load transform data
        transform_data = load(file_transformations);
        transform_data = transform_data.save_transform;
        if ~isempty(transform_data.transform_points{2})
            ud.pointList = transform_data.transform_points{2};
            title_ending = ' (transform points loaded)';
        end       
    end
    title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files) title_ending])    
