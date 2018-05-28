function transformed_sliceBrowser(transformed_slice_figure, transformed_images_folder, f, highlight_point, relevant_slice, min_dist, ...
                                                clickX, clickY, point_ind)


% go through histology at the same time
figure(transformed_slice_figure)
clf(transformed_slice_figure)
transformed_images = dir([transformed_images_folder '*.tif']);
total_num_files = size(transformed_images,1); 
ud_transformed_slice.total_num_files = total_num_files;
ud_transformed_slice.current_plot_handles = []; 

ud = get(f, 'UserData');

ud_transformed_slice.transformed_images = transformed_images;
ud_transformed_slice.transformed_images_folder = transformed_images_folder;
ud_transformed_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
hold(ud_transformed_slice.sliceAx, 'on');
set(ud_transformed_slice.sliceAx, 'HitTest', 'off');
title('Probe on Slice Viewer');

ud_transformed_slice.im = plotTVslice(zeros(800,1140, 'uint8'));


if highlight_point
    ud_transformed_slice.slice_num = relevant_slice;
    if min_dist < 40
        x_plot = ud.pointList{ud.currentProbe,1}(point_ind,1);
        y_plot = 800-ud.pointList{ud.currentProbe,1}(point_ind,2);
    else
        x_plot = clickX;
        y_plot = 800 - clickY;        
    end
        figure(transformed_slice_figure)
        ud_transformed_slice.quiver_plot{1} = scatter(ud_transformed_slice.sliceAx, x_plot, y_plot, 200, 'ro', ...
        'MarkerFaceColor',[1 1 1],'MarkerEdgeColor', [.3 .3 .3], ...
        'MarkerFaceAlpha',.5, 'MarkerFaceAlpha',.1);
    
    ud_transformed_slice.quiver_plot{2} = quiver( x_plot - 30 - 10, y_plot + 30 + 10, 30, -30, 1, 'color', 'white', 'linewidth',1,'MaxHeadSize',10);
else       
    ud_transformed_slice.slice_num = ud.curr_slice_num;
    ud_transformed_slice.quiver_plot = [];
end

processed_image_name = ud_transformed_slice.transformed_images(ud_transformed_slice.slice_num).name;
current_slice_image = flip(imread([ud_transformed_slice.transformed_images_folder processed_image_name]),1);
set(ud_transformed_slice.im, 'CData', current_slice_image);

set(transformed_slice_figure, 'KeyPressFcn', @(slice_figure, keydata)SliceAtlasHotkeyFcn(transformed_slice_figure, keydata, f));
set(transformed_slice_figure, 'UserData', ud_transformed_slice)



function SliceAtlasHotkeyFcn(fig, keydata, f)

ud = get(fig, 'UserData');
ud_atlas_viewer = get(f, 'UserData');

% set(ud.quiver_plot{1},'Visible','off');
% set(ud.quiver_plot{2},'Visible','off');

if strcmp(keydata.Key,'leftarrow')    
    if ud.slice_num > 1
        ud.slice_num = ud.slice_num - 1;
        
        % show next image
        processed_image_name = ud.transformed_images(ud.slice_num).name;
        current_slice_image = flip(imread([ud.transformed_images_folder processed_image_name]),1);
        set(ud.im, 'CData', current_slice_image);
        
        title(['Probe on Slice Viewer ' num2str(ud.slice_num) '/' num2str(ud.total_num_files)])
        
        % plot probe points for that slice
        set(ud.current_plot_handles(:), 'Visible', 'off'); ud.current_plot_handles = [];
        
        for probe = 1:size(ud_atlas_viewer.pointList,1)
            for point = 1:size(ud_atlas_viewer.pointList{probe,1},1)
                if ud.slice_num == ud_atlas_viewer.pointList{probe,2}(point)
                    ud.current_plot_handles(end+1) = scatter(ud.sliceAx, ud_atlas_viewer.pointList{probe,1}(point,1), ...
                        800-ud_atlas_viewer.pointList{probe,1}(point,2), 20, 'ro', 'MarkerFaceColor',[ .1 .1 .1],...
                    'MarkerEdgeColor', ud_atlas_viewer.ProbeColors(probe, :), 'MarkerFaceAlpha',.4,'LineWidth',1.5);
                end
            end
        end
             
            
            

    end
elseif strcmp(keydata.Key,'rightarrow') % break
    if ud.slice_num < ud.total_num_files
        ud.slice_num = ud.slice_num + 1;
        
        processed_image_name = ud.transformed_images(ud.slice_num).name;
        current_slice_image = flip(imread([ud.transformed_images_folder processed_image_name]),1);
        set(ud.im, 'CData', current_slice_image);
        
        title(['Probe on Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files)])        
        
        % plot probe points for that slice
        set(ud.current_plot_handles(:), 'Visible', 'off'); ud.current_plot_handles = [];
        
        for probe = 1:size(ud_atlas_viewer.pointList,1)
            for point = 1:size(ud_atlas_viewer.pointList{probe,1},1)
                if ud.slice_num == ud_atlas_viewer.pointList{probe,2}(point)
                    ud.current_plot_handles(end+1) = scatter(ud.sliceAx, ud_atlas_viewer.pointList{probe,1}(point,1), ...
                        800-ud_atlas_viewer.pointList{probe,1}(point,2), 20, 'ro', 'MarkerFaceColor',[ .1 .1 .1],...
                    'MarkerEdgeColor', ud_atlas_viewer.ProbeColors(probe, :), 'MarkerFaceAlpha',.4,'LineWidth',1.5);
                end
            end
        end


    end
elseif strcmp(keydata.Key,'uparrow') || strcmp(keydata.Key,'downarrow')
    % update probe points    
        % plot probe points for that slice
        set(ud.current_plot_handles(:), 'Visible', 'off'); ud.current_plot_handles = [];
        
        for probe = 1:size(ud_atlas_viewer.pointList,1)
            for point = 1:size(ud_atlas_viewer.pointList{probe,1},1)
                if ud.slice_num == ud_atlas_viewer.pointList{probe,2}(point)
                    ud.current_plot_handles(end+1) = scatter(ud.sliceAx, ud_atlas_viewer.pointList{probe,1}(point,1), ...
                        800-ud_atlas_viewer.pointList{probe,1}(point,2), 20, 'ro', 'MarkerFaceColor', [ .1 .1 .1], ... %ud_atlas_viewer.ProbeColors(probe, :),...
                    'MarkerEdgeColor', ud_atlas_viewer.ProbeColors(probe, :), 'MarkerFaceAlpha',.4,'LineWidth',1.5);
                end
            end
        end
    
end







set(fig, 'UserData', ud);


