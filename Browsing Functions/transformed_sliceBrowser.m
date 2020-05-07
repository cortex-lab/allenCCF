function transformed_slice_figure = transformed_sliceBrowser(transformed_slice_figure, save_location, f, highlight_point, relevant_slice, min_dist, ...
                                                clickX, clickY, point_ind, add)

% go through histology at the same time

transformed_images_folder = fullfile(save_location, ['transformations' filesep]);
transformed_images = dir([transformed_images_folder filesep '*.tif']);

try
    ud_transformed_slice.transformed_image_names = natsortfiles({transformed_images.name});
    browse_transformed_slices = true;
catch
    disp('no transformed images found.'); ud_transformed_slice.transformed_image_names = {};
    browse_transformed_slices = false;
end

if browse_transformed_slices

    try
        figure(transformed_slice_figure);
    catch
        transformed_slice_figure = figure('Name','Transformed Slice & Probe Point Viewer');
    end
    clf(transformed_slice_figure)


    ud_transformed_slice.current_plot_handles = []; 
    ud_transformed_slice.save_location = save_location;

    ud = get(f, 'UserData');

    ud_transformed_slice.transformed_images = transformed_images;
    ud_transformed_slice.transformed_images_folder = transformed_images_folder;
    ud_transformed_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
    hold(ud_transformed_slice.sliceAx, 'on');
    set(ud_transformed_slice.sliceAx, 'HitTest', 'off');
    title('Probe on Slice Viewer');

    ud_transformed_slice.im = plotTVslice(zeros(ud.ref_size, 'uint8'));

    % for probe view mode
    if highlight_point
        ud_transformed_slice.all_slices_slice_num = relevant_slice;
        if min_dist < 10
            x_plot = ud.pointList{ud.currentProbe,1}(point_ind,1);
            y_plot = ud.ref_size(1)-ud.pointList{ud.currentProbe,1}(point_ind,2);
        else
            x_plot = clickX;
            y_plot = ud.ref_size(1) - clickY;        
        end
            figure(transformed_slice_figure)
            ud_transformed_slice.quiver_plot{1} = scatter(ud_transformed_slice.sliceAx, x_plot, y_plot, 200, 'ro', ...
            'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [.3 .3 .3], ...
            'MarkerFaceAlpha',.5, 'MarkerFaceAlpha',.1);

        ud_transformed_slice.quiver_plot{2} = quiver( x_plot - 30 - 10, y_plot + 30 + 10, 30, -30, 1, 'color', 'white', 'linewidth',1,'MaxHeadSize',10);
    else       
        ud_transformed_slice.all_slices_slice_num = ud.slice_at_shift_start+ud.slice_shift+add;
        ud_transformed_slice.quiver_plot = [];
    end

    processed_images = dir([ud_transformed_slice.save_location filesep '*.tif']);
    processed_image_names = natsortfiles({processed_images.name});
    ud_transformed_slice.total_num_files = length(processed_image_names); 

    ud_transformed_slice.slice_num =  find(strcmp([processed_image_names{ud_transformed_slice.all_slices_slice_num}(1:end-4)  '_transformed.tif'],ud_transformed_slice.transformed_image_names));

    % show image
    if size(ud_transformed_slice.slice_num,2)
        processed_image_name = ud_transformed_slice.transformed_image_names{ud_transformed_slice.slice_num};
        current_slice_image = flip(imread([transformed_images_folder processed_image_name]),1);
        ud_transformed_slice.extra_text = ' (transformed)';
    else
        processed_image_name = processed_image_names{ud_transformed_slice.all_slices_slice_num};
        current_slice_image = flip(imread(fullfile(save_location, processed_image_name)),1);
        ud_transformed_slice.extra_text = ' (not transformed)';
    end


    set(ud_transformed_slice.im, 'CData', current_slice_image);
    title(['Probe on Slice Viewer ' num2str(ud_transformed_slice.all_slices_slice_num) '/' num2str(ud_transformed_slice.total_num_files) ud_transformed_slice.extra_text])

    set(transformed_slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(transformed_slice_figure, keydata, f, transformed_images_folder));
    set(transformed_slice_figure, 'UserData', ud_transformed_slice)
end

% ------------------------
% react to keyboard press
% ------------------------
function SliceAtlasHotkeyFcn(fig, keydata, f, transformed_images_folder)

ud = get(fig, 'UserData');
ud_atlas_viewer = get(f, 'UserData');

processed_images = dir([ud.save_location filesep '*.tif']);
processed_image_names = natsortfiles({processed_images.name});

transformed_images = dir([transformed_images_folder filesep '*.tif']);
try
    ud.transformed_image_names = natsortfiles({transformed_images.name});
catch
    ud.transformed_image_names = {};
end

ud.transformed_images = transformed_images;


if strcmp(keydata.Key,'leftarrow')    
    if ud.all_slices_slice_num > 1
        ud.all_slices_slice_num = ud.all_slices_slice_num - 1;
        
        ud.slice_num =  find(strcmp([processed_image_names{ud.all_slices_slice_num}(1:end-4)  '_transformed.tif'],ud.transformed_image_names));
        
        % show next image
        if size(ud.slice_num,2)
            processed_image_name = ud.transformed_image_names{ud.slice_num};
            current_slice_image = flip(imread(fullfile(ud.transformed_images_folder, processed_image_name)),1);
            ud.extra_text = ' (transformed)';
        else
            processed_image_name = processed_image_names{ud.all_slices_slice_num};
            current_slice_image = flip(imread(fullfile(ud.save_location, processed_image_name)),1);
            ud.extra_text = ' (not transformed)';
        end
        % reduce to 3 channels at most
        color_channels = min( 3, size(image,3));
        current_slice_image = current_slice_image(:,:,1:color_channels);
        set(ud.im, 'CData', current_slice_image);


    end
elseif strcmp(keydata.Key,'rightarrow') % break
    if ud.all_slices_slice_num < ud.total_num_files
       ud.all_slices_slice_num = ud.all_slices_slice_num + 1;
        
        ud.slice_num =  find(strcmp([processed_image_names{ud.all_slices_slice_num}(1:end-4)  '_transformed.tif'],ud.transformed_image_names));
        
        % show next image
        if size(ud.slice_num,2)
            processed_image_name = ud.transformed_image_names{ud.slice_num};
            current_slice_image = flip(imread(fullfile(ud.transformed_images_folder, processed_image_name)),1);
            ud.extra_text = ' (transformed)';
        else
            processed_image_name = processed_image_names{ud.all_slices_slice_num};
            current_slice_image = flip(imread(fullfile(ud.save_location, processed_image_name)),1);
            ud.extra_text = ' (not transformed)';
        end
        % reduce to 3 channels at most
        color_channels = min( 3, size(image,3));
        current_slice_image = current_slice_image(:,:,1:color_channels);
        set(ud.im, 'CData', current_slice_image);
    end
else
% otherwise -- call function to atlas browser    
    fcn = get(f, 'KeyPressFcn'); 
    fcn(f, keydata);    
end
        title(['Probe on Slice Viewer -- Slice ' num2str(ud.all_slices_slice_num) '/' num2str(ud.total_num_files) ud.extra_text ])        
        % plot probe points for that slice
        set(ud.current_plot_handles(:), 'Visible', 'off'); ud.current_plot_handles = [];
        
        for probe = 1:size(ud_atlas_viewer.pointList,1)
            for point = 1:size(ud_atlas_viewer.pointList{probe,1},1)
                if ud.all_slices_slice_num == ud_atlas_viewer.pointList{probe,2}(point)
                    ud.current_plot_handles(end+1) = scatter(ud.sliceAx, ud_atlas_viewer.pointList{probe,1}(point,1), ...
                        ud_atlas_viewer.ref_size(1)-ud_atlas_viewer.pointList{probe,1}(point,2), 20, 'ro', 'MarkerFaceColor',[ .1 .1 .1],...
                    'MarkerEdgeColor', ud_atlas_viewer.ProbeColors(probe, :), 'MarkerFaceAlpha',.4,'LineWidth',1.5);
                end
            end
        end

set(fig, 'UserData', ud);
