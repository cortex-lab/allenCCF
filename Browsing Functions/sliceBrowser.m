function sliceBrowser(slice_figure, processed_images_folder)
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

ud_slice.processed_images = processed_images;
ud_slice.processed_images_folder = processed_images_folder;
ud_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
hold(ud_slice.sliceAx, 'on');
set(ud_slice.sliceAx, 'HitTest', 'off');
ud_slice.im = plotTVslice(zeros(800,1140, 'uint8'));

% create functions needed to interact with the figure
set(ud_slice.im, 'ButtonDownFcn', @(slice_figure,k)sliceClickCallback(slice_figure, k));
set(slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(slice_figure, keydata));
set(slice_figure, 'UserData', ud_slice)

% adjust figure to user's screen size
try; screen_size = get(0,'ScreenSize'); screen_size = screen_size(3:4)./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end
set(slice_figure,'Position', [121*screen_size(1) 542*screen_size(2) 822*screen_size(1) 542*screen_size(2)])
movegui(slice_figure,'onscreen')

% set up first slice image
processed_image_name = ud_slice.processed_image_names{ud_slice.slice_num};
current_slice_image = flip(imread(fullfile(ud_slice.processed_images_folder, processed_image_name)),1);
if size(current_slice_image,1) > 802 || size(current_slice_image,1) > 1142
    disp('shrinking image to reference 800 x 1140 pxl')
    current_slice_image = imresize(current_slice_image, [800 1140]);
end  
set(ud_slice.im, 'CData', current_slice_image);
title('Slice Viewer');


% ------------------------------------------------    
% Clicking function to register transform points  
% ------------------------------------------------
function sliceClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');


if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));

    ud.pointList(end+1, :) = [clickX, 800 - clickY];
    ud.pointHands(end+1) = plot(ud.sliceAx, clickX, clickY, 'ro', 'color', [0 .5 0],'linewidth',2,'markers',4);    
    
     if clickX < 100 && (800 - clickY) < 100 % if click in corner, break
        ud.pointList = []; 
        set(ud.pointHands(:), 'Visible', 'off');     
     end
    
end
set(f, 'UserData', ud);

% ------------------------
% react to keyboard press
% ------------------------
function SliceAtlasHotkeyFcn(fig, keydata)

ud = get(fig, 'UserData');

% left arrow -- go to previous slice
if strcmp(keydata.Key,'leftarrow')    
    if ud.slice_num > 1
        ud.slice_num = ud.slice_num - 1;
        
        processed_image_name = ud.processed_image_names{ud.slice_num};
        current_slice_image = flip(imread(fullfile(ud.processed_images_folder, processed_image_name)),1);
        if size(current_slice_image,1) > 802 || size(current_slice_image,1) > 1142
            disp('shrinking image to reference size 800 x 1140 pxl')
            current_slice_image = imresize(current_slice_image, [800 1140]);
        end        
        set(ud.im, 'CData', current_slice_image); 
        
        title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files)])
        file_transformations = fullfile(ud.processed_images_folder, 'transformations\\' ,...
                                [processed_image_name(1:end-4) '_transform_data.mat']);
                            
        set(ud.pointHands(:), 'Visible', 'off'); 
        
        if exist(file_transformations,'file')
            % load transform data
            transform_data = load(file_transformations);
            transform_data = transform_data.save_transform;
            if ~isempty(transform_data.transform_points{2})
                ud.pointList = transform_data.transform_points{2};
            end
        else
            ud.pointList = [];
        end
    end
    
% left arrow -- go to next slice    
elseif strcmp(keydata.Key,'rightarrow') 
    if ud.slice_num < ud.total_num_files
        ud.slice_num = ud.slice_num + 1;
        
        processed_image_name = ud.processed_image_names{ud.slice_num};
        current_slice_image = flip(imread(fullfile(ud.processed_images_folder, processed_image_name)),1);
        if size(current_slice_image,1) > 802 || size(current_slice_image,1) > 1142
            disp('shrinking image to reference size 800 x 1140 pxl')
            current_slice_image = imresize(current_slice_image, [800 1140]);
        end          
        set(ud.im, 'CData', current_slice_image); 
        
        title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files)])
        
        file_transformations = fullfile(ud.processed_images_folder, 'transformations\\' ,...
                                [processed_image_name(1:end-4) '_transform_data.mat']);
                            
        set(ud.pointHands(:), 'Visible', 'off'); 
        
        if exist(file_transformations,'file')
            % load transform data
            transform_data = load(file_transformations);
            transform_data = transform_data.save_transform;
            if ~isempty(transform_data.transform_points{2})
                ud.pointList = transform_data.transform_points{2};
            end
        else
            ud.pointList = [];
        end
    end
% d -- delete current transform points
elseif strcmp(keydata.Key,'d') 
    disp('current transform points deleted')
    set(ud.pointHands(:), 'Visible', 'off'); 
    ud.pointList = [];    
elseif strcmp(keydata.Key,'t')
    ud.getPoint = ~ud.getPoint;
        if ud.getPoint; disp('transform point mode!'); end
end


set(fig, 'UserData', ud);
