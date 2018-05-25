function sliceBrowser(slice_figure, processed_images_folder)

% go through histology at the same time
processed_images = dir(processed_images_folder);
total_num_files = size(processed_images,1) - 2; disp(['found ' num2str(total_num_files) ' processed slice images']);
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

set(ud_slice.im, 'ButtonDownFcn', @(slice_figure,k)sliceClickCallback(slice_figure, k));
set(slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(slice_figure, keydata));
set(slice_figure, 'UserData', ud_slice)


processed_image_name = ud_slice.processed_images(ud_slice.slice_num + 2).name;
current_slice_image = flip(imread([ud_slice.processed_images_folder processed_image_name]),1);
set(ud_slice.im, 'CData', current_slice_image);


    
    
function sliceClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');


if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    
    ud.pointList(end+1, :) = [clickX, 800 - clickY];
    ud.pointHands(end+1) = plot(ud.sliceAx, clickX, clickY, 'ro', 'color', 'green','linewidth',2);
end
set(f, 'UserData', ud);




function SliceAtlasHotkeyFcn(fig, keydata)

ud = get(fig, 'UserData');

%display relevant points
% if ~isempty(ud.pointHands)
%        arrayfun(@(x)pointVis(ud.pointHands(x), ud.currentSlice, ud.pointList(x,:)), 1:numel(ud.pointHands));
% end

if strcmp(keydata.Key,'leftarrow')    
    if ud.slice_num > 1
        ud.slice_num = ud.slice_num - 1;
        
        processed_image_name = ud.processed_images(ud.slice_num + 2).name;
        current_slice_image = flip(imread([ud.processed_images_folder processed_image_name]),1);
        set(ud.im, 'CData', current_slice_image);
        ud.pointList = [];


    end
elseif strcmp(keydata.Key,'rightarrow') % break
    if ud.slice_num+1 < ud.total_num_files
        ud.slice_num = ud.slice_num + 1;
        
        processed_image_name = ud.processed_images(ud.slice_num + 2).name;
        current_slice_image = flip(imread([ud.processed_images_folder processed_image_name]),1);
        set(ud.im, 'CData', current_slice_image);
        set(ud.pointHands(:), 'Visible', 'off'); 
        ud.pointList = [];


    end
elseif strcmp(keydata.Key,'t')
    ud.getPoint = ~ud.getPoint;
        if ud.getPoint; disp('transform point mode!'); end
end







set(fig, 'UserData', ud);


