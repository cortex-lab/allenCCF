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
ud_slice.pointList = [];  % point coordinates
ud_slice.pointHands = []; % graphic objects
ud_slice.getPoint = 0;
ud_slice.ref_size = reference_size(2:3);

ud_slice.processed_images = processed_images;
ud_slice.processed_images_folder = processed_images_folder;
ud_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
hold(ud_slice.sliceAx, 'on');
set(ud_slice.sliceAx, 'HitTest', 'off');
ud_slice.im = plotTVslice(zeros(ud_slice.ref_size, 'uint8'));

ud_slice.pointsText = annotation('textbox', [0.88 0.03 0.1 0.05], ...
    'String', '0 point', 'EdgeColor', 'none', 'Color', 'k', 'HorizontalAlignment', 'right');
ud_slice.pointsText.Visible = 'off';

check_file_status(ud_slice)

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

set(slice_figure, 'UserData', ud_slice);


% ------------------------------------------------    
% Clicking function to register transform points  
% ------------------------------------------------
function sliceClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');

if isempty(ud.pointHands)
    ud.pointHands = gobjects(0);
end

if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));

    % ud.pointList(end+1, :) = [clickX, ud.ref_size(1) - clickY]; % often pointList is shorter than pointHands
    set(ud.pointHands,'Color',[.7 .3 .3])
    ud.pointHands(end+1) = plot(ud.sliceAx, clickX, clickY, 'o', 'color', [0 .9 0],'linewidth',2,'markers',4);

    % make sure what you see are what you have
    phx = zeros(length(ud.pointHands), 1);
    phy = zeros(length(ud.pointHands), 1);
    for i = 1:length(ud.pointHands)
        phx(i) = ud.pointHands(i).XData;
        phy(i) = ud.pointHands(i).YData;
    end
    ud.pointList = [phx, ud.ref_size(1) - phy]; %TODO

    if clickX < 100 && (ud.ref_size(1) - clickY) < 100 % if click in corner, break
        ud.pointList = [];
        set(ud.pointHands(:), 'Visible', 'off');
    end

    ud.pointsText.String = sprintf('%d point(s)', length(ud.pointHands));
    
end
set(f, 'UserData', ud);
end

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
    if isempty(ud.pointHands)
        disp('There is no transform point to delete.')
    else
        set(ud.pointHands(end), 'Visible', 'off'); %TODO isn't this better to be deleted? But if you do, you have a problem for callbacks        
        ud.pointHands = ud.pointHands(1:end-1);
        
        % make sure the length of ud.pointHands and the rows of ud.pointList match
        % make sure what you see are what you have
        phx = zeros(length(ud.pointHands), 1);
        phy = zeros(length(ud.pointHands), 1);
        for i = 1:length(ud.pointHands)
            phx(i) = ud.pointHands(i).XData;
            phy(i) = ud.pointHands(i).YData;
        end
        ud.pointList = [phx, ud.ref_size(1) - phy];
        if ~isempty(ud.pointHands)
            set(ud.pointHands(end),'color', [0 .9 0]);
        end

        disp('transform point deleted')
        ud.pointsText.String = sprintf('%d point(s)', length(ud.pointHands));

    end
% t -- transform point mode
elseif strcmp(keydata.Key,'t')
    ud.getPoint = ~ud.getPoint;
    if ud.getPoint
        disp('transform point mode!'); 
        set(ud.pointHands(:), 'Visible', 'on')
        ud.pointsText.Visible = 'on';
        ud.pointsText.String = sprintf('%d point(s)', length(ud.pointHands));    
    else
        ud.pointsText.Visible = 'off';
        set(ud.pointHands(:), 'Visible', 'off')
    end
else
% otherwise -- call function to atlas browser       
    figure(f);
    fcn = get(f, 'KeyPressFcn'); 
    fcn(f, keydata);
end


set(fig, 'UserData', ud);
end

function ud = updateSliceImage(ud)
    %TODO occasionally pointList has more items than pointHands

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

    try
        delete(findobj(ud.im.Parent, 'Type','line','Marker','o')); % delete all the circles
    end

    ud.pointHands = gobjects(0);
    ud.pointList = [];
    
    if exist(file_transformations,'file')
        % load transform data
        transform_data = load(file_transformations);
        transform_data = transform_data.save_transform;
        if ~isempty(transform_data.transform_points{2})
            ud.pointList = transform_data.transform_points{2};
            for i = 1:size(ud.pointList,1)
                ud.pointHands(end+1) = plot(ud.sliceAx, ud.pointList(i,1), ud.ref_size(1) - ud.pointList(i,2), 'ro', 'color', [0 .5 0],'linewidth',2,'markers',4);
            end
            set(ud.pointHands(end),'color', [0 .9 0]);

            title_ending = ' (transform points loaded)';
        end       
    end
    if ud.getPoint
        set(ud.pointHands(:), 'Visible', 'on');
        ud.pointsText.String = sprintf('%d point(s)', length(ud.pointHands));        
    else
        set(ud.pointHands(:), 'Visible', 'off'); %Hide because not in 't' mode
        ud.pointsText.String = sprintf('%d point(s)', length(ud.pointHands));        
    end
    title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files) title_ending])

end

end


function check_file_status(ud_slice)

proc_file_names = [string({ud_slice.processed_images.name})]';
processed_images_folder = ud_slice.processed_images_folder;


T_files = table(proc_file_names, zeros(length(proc_file_names),1),'VariableNames',{'file_names','status'});

if isfile(fullfile(processed_images_folder, 'file_status.mat'))

    S = load(fullfile(processed_images_folder, 'file_status.mat'));
    for i = 1:height(T_files)
        % apply saved values for the existing images
        T_files.status(i) = S.T_files{S.T_files.file_names == T_files.file_names(i) , 2};
    end

else
    error('You have to run HistologyBrowser first for the image %s', this_image_name)    

end

T_files.Properties.VariableDescriptions = ["",...
    "0, raw; 1, HistologyBrowser; 2, SliceFlipper"];

for i = 1:height(T_files)
    switch T_files.status(i) 
        case 0
            error('You have to run HistologyBrowser and SliceFlipper for the image %s', this_image_name)
        case 1
            error('You have to run SliceFlipper for the image %s', this_image_name)      
        case 2
            % good to go!
    end
end

end
