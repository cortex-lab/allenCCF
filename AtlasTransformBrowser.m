

function f = allenAtlasBrowser(templateVolume, annotationVolume, structureTree, slice_figure, save_location, save_suffix)
% Browser for the allen atlas ccf data in matlab.
%
% Inputs templateVolume, annotationVolume, and structureTree are the data describing the atlas.
% The annotation volume should be the "by_index" version
%


fprintf(1, 'Controls: \n');
fprintf(1, '--------- \n');
fprintf(1, 'scroll: move between slices \n');
fprintf(1, 'g: add/remove gridlines \n');
fprintf(1, 'o: add/remove overlay of current region extent \n');
fprintf(1, 'h: add/remove overlay of current histology slice \n');
fprintf(1, 'a: switch to viewing annotations (or switch back) \n');
fprintf(1, 'p: enable/disable mode where clicks are logged for probe or switch probes \n');
fprintf(1, 't: enable/disable mode where clicks are logged for transform \n');
fprintf(1, 'x: save transform \n');
fprintf(1, 'l: load transform for current slice \n');
fprintf(1, 'n: add a new probe \n');
fprintf(1, 's: save current probe \n');
fprintf(1, 'w: enable/disable probe viewer mode for current probe  \n');
fprintf(1, 'd: delete most recent probe point or all transform points \n');

fprintf(1, 'up: scroll through A/P angles \n');
fprintf(1, 'right: scroll through M/L angles \n');
fprintf(1, 'down: scroll through slices \n');

f = figure('Name','Atlas Viewer','Position',[956 401 883 657]); 

ud.bregma = allenCCFbregma; 

ud.currentSlice = ud.bregma(1); 
ud.currentAngle = zeros(2,1);
ud.scrollMode = 0;
ud.offset_map = zeros(800,1140); % size of reference image
ud.oldContour = [];
ud.showContour = false;
ud.showOverlay = false; ud.overlayAx = [];
ud.pointList = cell(1,3); ud.pointList{1} = zeros(0,3); 
ud.pointHands = cell(1,3);
ud.probe_view_mode = false;
ud.currentProbe = 0; ud.ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; 1 .6 0; .7 0 0; .5 0 .6]; 
ud.ProbeColor =  {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','orange','red','purple'};
ud.getPoint_for_transform = false; ud.pointList_for_transform = zeros(0,2); ud.pointHands_for_transform = [];
ud.current_pointList_for_transform = zeros(0,2); ud.curr_slice_num = 1;
ud.showAtlas = false;
ud.histology_overlay = 0; 
ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);
ud.transform = [];
ud.transformed_slice_figure = [];
ud.slice_shift = 0;
ud.loaded_slice = 0;
ud.slice_at_shift_start = 1;
ud.text = [];

ud.im = plotTVslice(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.ref = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_im = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_slice_trans = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.im_annotation = zeros(800,1140,'uint16');
ud.atlas_boundaries = zeros(800,1140,'uint16');;
ud.loaded = 10;

set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));

ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;

hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'off');

set(f, 'UserData', ud);

set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, slice_figure, k, allData, save_location, save_suffix));
set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData, slice_figure, save_location))
set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData, slice_figure, save_location)); % Set the motion detector.


function hotkeyFcn(f, slice_figure, keydata, allData, save_location, save_suffix)

ud = get(f, 'UserData');
ud_slice = get(slice_figure, 'UserData');

key_letter = lower(keydata.Key);
switch key_letter  
    case 'g' % toggle showing Gridlines
        if ~isfield(ud, 'gridlines') || isempty(ud.gridlines)
            axes(ud.atlasAx); hold on;
            gridY = 100:100:700; % assuming the size of the atlas for this for now
            gridX = 70:100:1140; 
            xl = xlim(); yl = ylim();
            gx = arrayfun(@(x)plot(x*[1 1], yl, 'w'), gridX, 'uni', false);
            gy = arrayfun(@(y)plot(xl, y*[1 1], 'w'), gridY, 'uni', false);
            ud.gridlines = [gx gy];
        elseif strcmp(get(ud.gridlines{1}, 'Visible'), 'on');
            cellfun(@(x)set(x, 'Visible', 'off'), ud.gridlines);
        elseif strcmp(get(ud.gridlines{1}, 'Visible'), 'off');
            cellfun(@(x)set(x, 'Visible', 'on'), ud.gridlines);
        end
    case 'o' % toggle Overlay
        ud.showOverlay = ~ud.showOverlay;
        if ~ud.showOverlay
            delete(ud.overlayAx); ud.overlayAx = [];
            disp('Overlay OFF');
        else; disp('Overlay on!');
        end
    case 'p' % toggle mode to register clicks as Probe points
        ud.probe_view_mode = 0;
        if ud.currentProbe < size(ud.pointList,1)
            ud.currentProbe = ud.currentProbe + 1;
            
            disp(['probe point mode -- selecting probe ' num2str(ud.currentProbe) ' (' ud.ProbeColor{ud.currentProbe} ')']); 
            ud.getPoint_for_transform = false; 
            
            % show Transformed Slice & Probage Viewer, if not already showing
            if ~ud.slice_at_shift_start; add = 1; else; add = 0; end
        slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift+add}(1:end-4);
        folder_transformations = [save_location 'transformations\\'];
            try; load([folder_transformations slice_name '_transform_data.mat']);
                try; figure(ud.transformed_slice_figure); 
                catch
                    ud.transformed_slice_figure = figure('Name','Transformed Slice & Probe Point Viewer','Position', [265 37 560 420]);
                    highlight_point = false;
                    transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, highlight_point, [], [], [], [], [], add)
                end; figure(f);
            end
        else
            ud.currentProbe = 0;
            disp(['probe point mode OFF']);
        end
    case 'w'
        for probe_plotted = 1:size(ud.pointHands,1)
            set(ud.pointHands{probe_plotted,1}(:), 'Visible', 'off'); 
            set(ud.pointHands{probe_plotted,3}(:), 'Visible', 'off'); 
        end
        if ~ud.currentProbe
            ud.currentProbe = 1;
        end
            ud.probe_view_mode = ~ud.probe_view_mode; 
            
            if ud.probe_view_mode
                 % load probe points if none are already up
                if ~size(ud.pointList{1,1},1)
                        probe_points = load([save_location 'probe_points' save_suffix]);  disp('loading probe points')
                        ud.pointList = probe_points.pointList.pointList;
                        ud.pointHands = probe_points.pointList.pointHands;
                end; disp(['activate probe view mode for probe ' num2str(ud.currentProbe)])
                
                for probe_point = 1:size(ud.pointHands{ud.currentProbe, 1}(:),1)
                    ud.pointHands{ud.currentProbe, 1}(probe_point) = scatter(ud.atlasAx, ...
                        ud.pointList{ud.currentProbe,1}(probe_point,1), ud.pointList{ud.currentProbe,1}(probe_point,2), 20, 'ro', ...
                    'MarkerFaceColor', ud.ProbeColors(ud.currentProbe, :),'MarkerEdgeColor', ud.ProbeColors(ud.currentProbe, :), ...
                    'LineWidth',3);
                end
                set(ud.pointHands{ud.currentProbe,1}, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));
                
                curr_probePoints = ud.pointList{ud.currentProbe,1}(:, [3 2 1]);

                % get line of best fit through points
                % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
                [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
                
                min_y = min(ud.pointList{ud.currentProbe,1}(:,2));
                max_y = max(ud.pointList{ud.currentProbe,1}(:,2));
                min_x = m(3) + (min_y - m(2))  * p(3) / p(2);
                max_x = m(3) + (max_y - m(2))  * p(3) / p(2);
                ud.pointHands{ud.currentProbe, 3} = plot([min_x max_x],[min_y max_y],'color',ud.ProbeColors(ud.currentProbe,:),'linestyle',':');
                set(ud.pointHands{ud.currentProbe,3}, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));
                
                % ensure proper orientation: want 0 at the top of the brain 
                % and positive distance goes down into the brain
                if p(2)<0
                    p = -p;
                end

                % use the z value of the mean as the selected slice
                ud.currentSlice = round(m(1)); 

                % calculate slice angle along probe track
                ud.currentAngle(1) = round (400*p(1) / p(2) ) ;
                
                z_shift_above_probe_center = m(2)*p(1) / p(2);
                ud.currentAngle(2) = round( 570 * (ud.currentAngle(1) - z_shift_above_probe_center) / (m(3) - 570) );
                p_abs = abs(p);
                ud.currentAngle(2) = round(ud.currentAngle(2) * ( p_abs(3) / (p_abs(1) + p_abs(3))));
                  

                % update slice
                update.VerticalScrollCount = 0; ud.scrollMode = 0; ud.histology_overlay = 0; set(f, 'UserData', ud);
                updateSlice(f, update, allData, slice_figure, save_location); ud = get(f, 'UserData');   
                fill([5 5 200 200],[5 50 50 5],[0 0 0]);

                
                % show Transformed Slice & Probage Viewer, if not already showing
                try; figure(ud.transformed_slice_figure); 
                catch
                    ud.transformed_slice_figure = figure('Name','Transformed Slice & Probe Point Viewer');
                    highlight_point = false;
                    transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, highlight_point, [], [], [], [], [],0)
                end; figure(f);                
            
                set(ud.im, 'CData', ud.ref);
                ud.curr_im = ud.ref; set(f, 'UserData', ud);
            else; disp('probe view mode OFF'); end

            
        
    case 't' % toggle mode to register clicks as Points -- 0, 1, or 2
        
        if ~size(ud.current_pointList_for_transform,1)
            ud.getPoint_for_transform = ud.getPoint_for_transform + 1 - 3*(ud.getPoint_for_transform==2);
        else
            ud.getPoint_for_transform = ~ud.getPoint_for_transform;
        end
            
        ud.loaded = false;
        
        if ud.getPoint_for_transform
            if ud.getPoint_for_transform
                disp('transform point mode on -- press t again before clicking for automatic transform point mode'); 
            else
                disp('automatic transform point mode activated (press top left of slice to abandon)!')
            end
                ud.currentProbe = 0;

            % launch transform point mode
            if ud_slice.slice_num ~= (ud.slice_at_shift_start+ud.slice_shift) || ~size(ud.current_pointList_for_transform,1) % ud_slice.slice_num ~= ud.curr_slice_num
                ud.curr_slice_num = ud_slice.slice_num;
                ud.current_pointList_for_transform = zeros(0,2);
                set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
                num_hist_points = size(ud_slice.pointList,1);
                template_point = 1; template_points_shown = 0;
                updateBoundaries(f,ud, allData); ud = get(f, 'UserData');

                while true && ud.getPoint_for_transform == 2
                    % select template point automatically
                    if template_point == 1
                        findX = 1140 / 2;
                        findY = min(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 4
                        findX = 1140 / 2;
                        findY = max(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 3
                        findX = 1140 / 2;
                        findY = round(prctile(find(ud.atlas_boundaries(:,findX)) , 66));   %median(find(ud.atlas_boundaries(:,findX)));            
                    elseif template_point == 2
                        findX = 1140 / 2;
                        findY = round(prctile(find(ud.atlas_boundaries(:,findX)) , 25));                       
                    elseif template_point == 12
                        findY = 800 / 2;
                        findX = max(find(ud.atlas_boundaries(findY,:)));
                    elseif template_point == 7
                        findY = 800 / 2;
                        findX = min(find(ud.atlas_boundaries(findY,:)));         
                    elseif template_point == 10
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 66));
                        findY = min(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 14
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 66));
                        findY = max(find(ud.atlas_boundaries(:,findX)));     
                    elseif template_point == 5
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 33));
                        findY = min(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 9
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 33));
                        findY = max(find(ud.atlas_boundaries(:,findX)));   
                    elseif template_point == 6
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 15));
                        findY = min(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 8
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 15));
                        findY = max(find(ud.atlas_boundaries(:,findX)));                           
                    elseif template_point == 11
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 85));
                        findY = min(find(ud.atlas_boundaries(:,findX)));
                    elseif template_point == 13
                        findX = round(prctile(find(ud.atlas_boundaries(800/2,:)) , 85));
                        findY = max(find(ud.atlas_boundaries(:,findX)));                           
                    end

                    % apply to template
                    if template_points_shown < template_point
                        ud.pointList_for_transform(end+1, :) = [findX, findY];
                        ud.current_pointList_for_transform(end+1, :) = [findX, findY];
                        set(ud.pointHands_for_transform(:), 'color', [0 0 0]);
                        ud.pointHands_for_transform(end+1) = plot(ud.atlasAx, findX, findY, 'ro', 'color', 'green','LineWidth',3);  
                        template_points_shown = template_points_shown + 1;
                    end
                    % move on if point was selected on histology
                    ud_slice = get(slice_figure, 'UserData');
                    if size(ud_slice.pointList,1) > num_hist_points
                         if ud_slice.pointList(end,1) < 100 && ud_slice.pointList(end,2) < 100% if click in corner, break
                            ud.current_pointList_for_transform = zeros(0,2); 
                            set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
                            disp('exit transform mode'); break; end
                        num_hist_points = size(ud_slice.pointList,1);
                        template_point = template_point + 1;
                    end

                    % stop after last point
                    if template_point > 14
                        disp('press h to transform histology (whereupon more points can be added if necessary)')
                        break
                    end

                    pause(.1);
                end
            end
        else; disp('transform point mode OFF');    
        end        
    case 'a' % toggle View between template/annotation
        
        ud.showAtlas = ~ud.showAtlas;
        
        if ud.showAtlas % superimpose boundaries           
            updateBoundaries(f,ud, allData);
        else % return to image
            set(ud.im, 'CData', ud.curr_im)
        end

    case 'uparrow' % scroll angles along M/L axis
        ud.scrollMode = 1;
        disp('switch scroll mode -- tilt D/V')
    case 'rightarrow' % scroll angles along A/P axis
        ud.scrollMode = 2;
        disp('switch scroll mode -- tilt M/L')
    case 'downarrow' % scroll along A/P axis
        ud.scrollMode = 0;
        disp('switch scroll mode -- scroll along A/P axis')
    case 'leftarrow' % scroll along A/P axis
        ud.scrollMode = 3;
        if ~ud.slice_at_shift_start
            ud.slice_at_shift_start = ud_slice.slice_num;
        end
        disp('switch scroll mode -- scroll along slice images')        
    case 'h'
        ud.histology_overlay = ud.histology_overlay + 1 - 3*(ud.histology_overlay==2);
        slice_points = ud_slice.pointList;
        
        slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
        folder_transformations = [save_location 'transformations\\'];
        if size(ud.current_pointList_for_transform,1)  && size(slice_points,1) && ud.slice_at_shift_start+ud.slice_shift == ud_slice.slice_num
            key_letter = 'x'; % save transform automatically
        end
        
        if (ud.histology_overlay == 1 || ud.histology_overlay == 2) && ...
                ( (size(ud.current_pointList_for_transform,1) && size(slice_points,1)) || ud.loaded)

        if ud.slice_shift > 0
            
            ud.curr_slice_trans = imread([folder_transformations slice_name '_transformed.tif']);
            
        else
                set(ud.text,'Visible','off');
                fill([5 5 200 200],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');            

            reference_points = ud.current_pointList_for_transform;
            slice_points = ud_slice.pointList;
            
            current_slice_image = flip(get(ud_slice.im, 'CData'),1);
            if ~ud.loaded  % use loaded version if 'l' was just pressed 
                ud.transform = fitgeotrans(slice_points,reference_points,'projective'); %can use 'affine', 'projective', or 'pwl'
            end
            R = imref2d(size(ud.ref));
            ud.curr_slice_trans = imwarp(current_slice_image, ud.transform, 'OutputView',R);
        end
            
            image_blend =  imfuse(uint8(ud.ref*.6), ud.curr_slice_trans(:,:,:),'blend','Scaling','none');
            if ud.histology_overlay == 2 % 2 ~ blend
                disp('Slice + Reference mode!');
                set(ud.im, 'CData', image_blend);
                ud.curr_im = image_blend;
            else % 1 ~ just see slice
                disp('Slice mode!');
                set(ud.im, 'CData', ud.curr_slice_trans); 
                ud.curr_im = ud.curr_slice_trans;
            end

        else % ud.histology_overlay == 0
            ud.histology_overlay = 0;
            disp('Reference mode!');
            set(ud.im, 'CData', ud.ref);
            ud.curr_im = ud.ref; set(f, 'UserData', ud);
        end
        if ud.showAtlas
            updateBoundaries(f,ud, allData);
        end
    case 'n' % new probe
        new_num_probes = size(ud.pointList,1) + 1; disp(['probe ' num2str(new_num_probes) ' added! (' ud.ProbeColor{new_num_probes} ')']);
        probe_point_list = cell(new_num_probes,1); probe_hands_list = cell(new_num_probes,3); 
        for prev_probe = 1:new_num_probes-1
            probe_point_list{prev_probe,1} = ud.pointList{prev_probe,1};
            probe_point_list{prev_probe,2} = ud.pointList{prev_probe,2};
            probe_point_list{prev_probe,3} = ud.pointList{prev_probe,3};
            probe_hands_list{prev_probe, 1} = ud.pointHands{prev_probe, 1};
            probe_hands_list{prev_probe, 2} = ud.pointHands{prev_probe, 2};
            probe_hands_list{prev_probe, 3} = ud.pointHands{prev_probe, 3};
        end; probe_point_list{new_num_probes,1} = zeros(0,3);
        ud.pointList = probe_point_list; ud.pointHands = probe_hands_list;
        ud.currentProbe = new_num_probes;
    case 's' % save probe trajectory and points of each probe per histology image (and associated histology name/number)
        pointList.pointList = ud.pointList;
        pointList.pointHands = ud.pointHands;
        save([save_location  'probe_points' save_suffix], 'pointList'); disp('probe points saved');        
    case 'l' % load transform and current slice position and angle
        slice_name = ud_slice.processed_image_names{ud_slice.slice_num}(1:end-4);
        folder_transformations = [save_location 'transformations\\'];
        
        try
            
        if ud.loaded_slice+ud.slice_shift ~= ud_slice.slice_num
            
            ud.curr_slice_num = ud_slice.slice_num;
            
            % load transform data
            transform_data = load([folder_transformations slice_name '_transform_data.mat']);  
            transform_data = transform_data.save_transform;

            % load new transform
            ud.transform = transform_data.transform;
           
            if ~isempty(transform_data.transform_points{1}) && ~isempty(transform_data.transform_points{2})
                ud.current_pointList_for_transform = transform_data.transform_points{1};
                ud_slice.pointList = transform_data.transform_points{2};
                set(slice_figure, 'UserData', ud_slice);
            end            
            
            % load allen ref location
            ud.currentSlice = transform_data.allen_location{1}; ud.currentAngle = transform_data.allen_location{2};

            % create transformed histology image
            current_slice_image = flip(get(ud_slice.im, 'CData'),1); R = imref2d(size(ud.ref));
            ud.curr_slice_trans = imwarp(current_slice_image, ud.transform, 'OutputView',R);

            % update figure
            update.VerticalScrollCount = 0; temp_scroll_mode = ud.scrollMode; ud.scrollMode = 4; set(f, 'UserData', ud);
            updateSlice(f, update, allData, slice_figure, save_location); ud = get(f, 'UserData');
            ud.scrollMode = temp_scroll_mode;
            ud.loaded = true;

            ud.slice_at_shift_start = ud_slice.slice_num; ud.slice_shift = 0;
            ud.loaded_slice = ud_slice.slice_num;
            
            if ~isempty(ud.text)
                set(ud.text,'Visible','off');
                fill([5 5 200 200],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');            
            end
            
            disp('transform loaded -- press ''l'' again now to load probe points');
        else % load probe points
            if ~size(ud.pointList{1,1},1)
                probe_points = load([save_location 'probe_points' save_suffix]);  disp('probe points loaded')
                ud.pointList = probe_points.pointList.pointList;
                ud.pointHands = probe_points.pointList.pointHands;
            end
        end
        
        for probe = 1:size(ud.pointList,1)
             set(ud.pointHands{probe, 1}(:),'Visible','off')
            for probe_point = 1:size(ud.pointHands{probe, 1}(:),1)
                ud.pointHands{probe, 1}(probe_point) = scatter(ud.atlasAx, ...
                    ud.pointList{probe,1}(probe_point,1), ud.pointList{probe,1}(probe_point,2), 30, 'ro', ...
                'MarkerFaceColor', [0 0 0],'MarkerEdgeColor', ud.ProbeColors(probe, :), ...
               'LineWidth',3);
            end
             
            set(ud.pointHands{probe, 3}(:),'Visible','off'); ud.pointHands{probe, 3} = [];
            for probe_point = 1:size(ud.pointList{probe,1},1)
                slice_point_belongs_to = ud.pointHands{probe, 2}(probe_point);
                if slice_point_belongs_to == ud_slice.slice_num
                    set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'on');
                else
                    set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'off'); 
                end

            end
    end
         
            
        
        ud.slice_shift = 0;
        catch; 
            disp(['loading failed']); end
            
    case 'd' % delete current transform or most recent probe point
        if ud.getPoint_for_transform
            ud.current_pointList_for_transform = zeros(0,2); set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
            ud.pointHands_for_transform = []; ud_slice.pointList = []; set(slice_figure, 'UserData', ud_slice);
            disp('current transform erased');
        elseif ud.currentProbe
            ud.pointList{ud.currentProbe,1} = ud.pointList{ud.currentProbe,1}(1:end-1,:);
            ud.pointList{ud.currentProbe,2} = ud.pointList{ud.currentProbe,2}(1:end-1,:);
            ud.pointList{ud.currentProbe,3} = ud.pointList{ud.currentProbe,3}(1:end-1,:);
            set(ud.pointHands{ud.currentProbe, 1}(end), 'Visible', 'off'); 
            ud.pointHands{ud.currentProbe, 1} = ud.pointHands{ud.currentProbe, 1}(1:end-1);
            ud.pointHands{ud.currentProbe, 2} = ud.pointHands{ud.currentProbe, 2}(1:end-1);
            disp('probe point deleted')
        end
end


if strcmp(key_letter,'x') % save transform and current slice position and angle
        
        % find or create folder location for transformations
        try
        folder_transformations = [save_location 'transformations\\'];
        if ~exist(folder_transformations)
            mkdir(folder_transformations)
        end

        slice_name = ud_slice.processed_image_names{ud_slice.slice_num}(1:end-4);
        
        % store transform
        save_transform.transform = ud.transform;
        % store transform points
        transform_points = cell(2,1); transform_points{1} = ud.current_pointList_for_transform;
        if ~isempty(ud_slice.pointList)
            transform_points{2} = ud_slice.pointList;
        end
        
        save_transform.transform_points = transform_points;
        % store reference location
        allen_location = cell(2,1); allen_location{1} = ud.currentSlice; allen_location{2} = ud.currentAngle; 
        save_transform.allen_location = allen_location;
        % save all this
        save([folder_transformations slice_name '_transform_data.mat'], 'save_transform');
        
        % save transformed histology image
        current_slice_image = flip(get(ud_slice.im, 'CData'),1); R = imref2d(size(ud.ref));
        curr_slice_trans = imwarp(current_slice_image, ud.transform, 'OutputView',R);
        imwrite(curr_slice_trans, [folder_transformations slice_name '_transformed.tif'])
        
        disp('transform and atlas location saved.')
        catch
            disp('transform not saved')
        end
end
        
set(f, 'UserData', ud);


function updateSlice(f, evt, allData, slice_figure, save_location)

ud = get(f, 'UserData');

% scroll through slices
if ud.scrollMode==0
    ud.currentSlice = ud.currentSlice+evt.VerticalScrollCount*3;

    if ud.currentSlice>size(allData.tv,1); ud.currentSlice = 1; end %wrap around
    if ud.currentSlice<1; ud.currentSlice = size(allData.tv,1); end %wrap around
    
    
% scroll through A/P angles        
elseif ud.scrollMode==1 %&&  abs(ud.currentAngle(1)) < 130
  ud.currentAngle(1) = ud.currentAngle(1)+evt.VerticalScrollCount*3;

% scroll through M/L angles
elseif ud.scrollMode==2 %&&  abs(ud.currentAngle(2)) < 130
  ud.currentAngle(2) = ud.currentAngle(2)+evt.VerticalScrollCount*3; 
  
% scroll through slices
elseif ud.scrollMode == 3
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
  
  ud_slice = get(slice_figure, 'UserData');
  
  try
  ud.slice_shift = ud.slice_shift-evt.VerticalScrollCount;
  slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
  catch
  ud.slice_shift = ud.slice_shift+evt.VerticalScrollCount;    
  slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
  end
  folder_transformations = [save_location 'transformations\\'];
  
%   ud.slice_at_shift_start = ud.loaded_slice; %ud_slice.slice_num;
  
    % set probe points from other slices invisible and from this slice visible
    for probe = 1:size(ud.pointList,1)
        set(ud.pointHands{probe, 3}(:),'Visible','off'); ud.pointHands{probe, 3} = [];     
        for probe_point = 1:size(ud.pointList{probe,1},1)
            slice_point_belongs_to = ud.pointList{probe, 2}(probe_point);
            if slice_point_belongs_to == ud.slice_at_shift_start+ud.slice_shift
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'on');
            else
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'off'); 
            end

        end
    end  
  
    try; load([folder_transformations slice_name '_transform_data.mat']);
       
        % load transform data
        transform_data = load([folder_transformations slice_name '_transform_data.mat']);  
        transform_data = transform_data.save_transform;
        
        % load new transform
        ud.transform = transform_data.transform;
        
        if ~isempty(transform_data.transform_points{1}) && ~isempty(transform_data.transform_points{2})
            ud.current_pointList_for_transform = transform_data.transform_points{1};
            ud_slice.pointList = transform_data.transform_points{2};
        else
            ud_slice.pointList = [];
        end
        set(slice_figure, 'UserData', ud_slice);
        
        % load allen ref location
        ud.currentSlice = transform_data.allen_location{1}; ud.currentAngle = transform_data.allen_location{2};

        % create transformed histology image
        ud.curr_slice_trans = imread([folder_transformations slice_name '_transformed.tif']);
        
        % update figure
        update.VerticalScrollCount = 0; set(f, 'UserData', ud);
%         updateSlice(f, update, allData, slice_figure, save_location); ud = get(f, 'UserData');
        ud.loaded = true;
        
        ud.histology_overlay = 1;
        
        set(ud.text,'Visible','off');
        fill([5 5 200 200],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');
    catch;
        % if no transform, just show reference
        ud.histology_overlay = 0;
        set(ud.im, 'CData', ud.ref);
        ud.curr_im = ud.ref; set(f, 'UserData', ud);   
        set(ud.text,'Visible','off');
        fill([5 5 200 200],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift) ' - no transform found'],'color','white');        
    end  
        
end  


% ---------------------------------
% if no angle, just do normal thing
% ---------------------------------
if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
    
    set(ud.im, 'CData', squeeze(allData.tv(ud.currentSlice,:,:)));
 
   % update title/overlay with brain region
    pixel = getPixel(ud.atlasAx);
    updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText);
    [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);
    updateTitle(ud.atlasAx, name, acr)
    if ud.showOverlay    
        updateOverlay(f, allData, ann, slice_figure, save_location);
    end  
    ud.ref = uint8(get(ud.im, 'CData'));
    set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
    set(f, 'UserData', ud);
% ---------------------------
% if angle, do angle thing    
% ---------------------------
else 
  image_size = size(squeeze(allData.av(ud.currentSlice,:,:)));
  angle_slice = zeros(image_size);
  
  if ud.currentAngle(1)==0; offset_AP = 0;
  else; offset_AP = -ud.currentAngle(1):sign(ud.currentAngle(1)):ud.currentAngle(1);
  end; start_index_AP = 1; 
 
  
  % loop through AP offsets
  for cur_offset_AP = offset_AP
      if cur_offset_AP == ud.currentAngle(1); end_index_AP = image_size(1);
      else; end_index_AP = start_index_AP + floor( image_size(1) / length(offset_AP)) - 1;
      end
      
       if ud.currentAngle(2)==0;  offset_ML = 0;
       else; offset_ML = -ud.currentAngle(2):sign(ud.currentAngle(2)):ud.currentAngle(2);
       end; start_index_ML = 1;
    % nested: loop through ML offsets
    for cur_offset_ML = offset_ML
      if cur_offset_ML == ud.currentAngle(2)
         end_index_ML = image_size(2);
      else
         end_index_ML = start_index_ML + floor( image_size(2) / length(offset_ML)) - 1;
      end
          
      % update current slice
     angle_slice(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = ...
         squeeze(allData.tv(ud.currentSlice + cur_offset_AP + cur_offset_ML,start_index_AP:end_index_AP,start_index_ML:end_index_ML));

    
      ud.im_annotation(start_index_AP:end_index_AP,start_index_ML:end_index_ML) = squeeze(allData.av(ud.currentSlice + cur_offset_AP + cur_offset_ML,...
                                                            start_index_AP:end_index_AP,start_index_ML:end_index_ML));
      ud.offset_map(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = cur_offset_AP + cur_offset_ML;
      
      start_index_ML = end_index_ML + 1;
    end
      start_index_AP = end_index_AP + 1;
  end     
  if ~ud.showAtlas
    set(ud.im, 'CData', angle_slice);
  end

  ud.ref = uint8(angle_slice);
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
end


% update histology overlay
if ud.histology_overlay == 1 || ud.histology_overlay == 2
    updateHistology(f,ud); ud = get(f, 'UserData');
else
    ud.curr_im = ud.ref;
end


% then update boundary overlay
if ud.showAtlas
    updateBoundaries(f,ud, allData)
end

% in normal mode, show only the probe points for that slice
% in probe view mode, color depends on distance from point
if ud.probe_view_mode && ud.currentProbe
    probe = ud.currentProbe; mean_distance = 0;
    for probe_point = 1:size(ud.pointList{probe,1},1)
        distance = abs(ud.pointList{probe,1}(probe_point,3) - ...
            (ud.currentSlice + ud.offset_map(ud.pointList{probe,1}(probe_point,2),ud.pointList{probe,1}(probe_point,1))));
        mean_distance = mean_distance + distance / size(ud.pointList{probe,1},1);
        if distance < 50
            color = abs((ud.ProbeColors(probe,:) * (50 - distance) + [0 0 0] * distance) / 50);
        else
            color = [.33 .33 .33];
        end
            set(ud.pointHands{probe, 1}(probe_point), 'MarkerFaceColor',color,'MarkerEdgeColor', color, 'SizeData', 20)
    end
disp(['mean distance from this slice to probe points is ' num2str(round(mean_distance*10)) ' microns'])
    if mean_distance < 50
        color = abs((ud.ProbeColors(probe,:) * (50 - mean_distance) + [0 0 0] * mean_distance) / 50);
    else
        color = [.33 .33 .33];
    end    
    set(ud.pointHands{probe, 3}, 'Color',color)
    
end
  set(f, 'UserData', ud);


% update the image shown if histology is currently being overlaid
function updateHistology(f, ud)
    if ud.histology_overlay == 2
        image_blend =  imfuse(uint8(ud.ref*.6), ud.curr_slice_trans(:,:,:),'blend','Scaling','none');
        set(ud.im, 'CData', image_blend);
        ud.curr_im = image_blend;
    elseif ud.histology_overlay == 1
        set(ud.im, 'CData', ud.curr_slice_trans);
        ud.curr_im = ud.curr_slice_trans;
    end
    set(f, 'UserData', ud);

% update the position of the region boundary image
function updateBoundaries(f, ud, allData)
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        curr_annotation = squeeze(allData.av(ud.currentSlice,:,:));
    else
        curr_annotation = ud.im_annotation;
    end
    
    atlas_vert_1 = double(curr_annotation(1:end-2,:));
    atlas_vert_2 = double(curr_annotation(3:end,:));
    atlas_vert_offset = abs( atlas_vert_1 - atlas_vert_2 ) > 0;
    shifted_atlas_vert1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert1(3:end,:) = atlas_vert_offset;
    shifted_atlas_vert2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert2(1:end-2,:) = atlas_vert_offset;

    atlas_horz_1 = double(curr_annotation(:,1:end-2));
    atlas_horz_2 = double(curr_annotation(:,3:end));
    atlas_horz_offset = abs( atlas_horz_1 - atlas_horz_2 )>0;
    shifted_atlas_horz1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz1(:,3:end) = atlas_horz_offset;
    shifted_atlas_horz2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz2(:,1:end-2) = atlas_horz_offset;

    shifted_atlas = shifted_atlas_horz1 + shifted_atlas_horz2 + shifted_atlas_vert1 + shifted_atlas_vert2;

    atlas_boundaries = (shifted_atlas>0); ud.atlas_boundaries = atlas_boundaries;

    if ud.showAtlas
        image_blend =  uint8( imfuse(ud.curr_im, atlas_boundaries/3.5*(1+.35*isa(ud.curr_im,'uint16')),'blend','Scaling','none') )* 2;
        set(ud.im, 'CData', image_blend); 
    end
    
    set(f, 'UserData', ud);
    

function atlasClickCallback(im, keydata, slice_figure, save_location)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');
ud_slice = get(slice_figure, 'UserData');

if ud.probe_view_mode && ud.currentProbe
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
    
    if ud.showOverlay 
        clickY = size(ud.ref,1) - clickY;
    end
    % find the probe point closest to this clicked point
    [min_dist, point_ind] = min( sqrt(sum(([clickX clickY clickZ] - ud.pointList{ud.currentProbe}).^2,2)));
    
    % find the slice corresponding to that point
    relevant_slice = ud.pointList{ud.currentProbe,2}(point_ind);
    highlight_point = true;
    transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, ...
        highlight_point, relevant_slice, min_dist, clickX, clickY, point_ind,0)
    figure(f);
    
elseif ud.currentProbe > 0
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
       
    ud.pointList{ud.currentProbe,1}(end+1, :) = [clickX, clickY, clickZ];
    ud.pointList{ud.currentProbe,2}(end+1, :) = ud.slice_at_shift_start + ud.slice_shift;
    slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start + ud.slice_shift}(1:end-4);
    ud.pointList{ud.currentProbe,3}{end+1} = slice_name;
   
    
    ud.pointHands{ud.currentProbe, 1}(end+1) = scatter(ud.atlasAx, clickX, clickY, 30, 'ro', ...
                'MarkerFaceColor', [0 0 0],'MarkerEdgeColor', ud.ProbeColors(ud.currentProbe, :), ...
                'LineWidth',3);
        
                
    ud.pointHands{ud.currentProbe, 2}(end+1) = ud_slice.slice_num + ud.slice_shift;
    
elseif ud.getPoint_for_transform
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    
    if ud.curr_slice_num ~= ud_slice.slice_num;
        ud.current_pointList_for_transform = zeros(0,2);
        ud_slice.pointList = []; set(slice_figure, 'UserData', ud_slice);
        ud.curr_slice_num = ud_slice.slice_num;
        disp('transforming new slice');
    end
    
    ud.pointList_for_transform(end+1, :) = [clickX, clickY];
    ud.current_pointList_for_transform(end+1, :) = [clickX, clickY];
    set(ud.pointHands_for_transform(:), 'color', [0 0 0]);
    ud.pointHands_for_transform(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro', 'color', [0 .9 0],'LineWidth',2,'markers',4);    
        
    ud.slice_at_shift_start = ud_slice.slice_num;
    ud.slice_shift = 0;
    ud.loaded = 0;
end
set(f, 'UserData', ud);


function fh_wbmfcn(f, allData, slice_figure, save_location)
% WindowButtonMotionFcn for the figure.

ud = get(f, 'UserData');

ax = ud.atlasAx;

pixel = getPixel(ax);

%get offset due to angling
if 0<pixel(1) && pixel(1)<=800 && 0<pixel(2) && pixel(2)<=1140
    offset = ud.offset_map(pixel(1),pixel(2));
else; offset = 0;
end

% show bregma coords
updateStereotaxCoords(ud.currentSlice + offset, pixel, ud.bregma, ud.bregmaText);

% get annotation for this pixel
[name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice+offset);

updateTitle(ax, name, acr);

if ~isempty(name)
    if ud.showContour
        if ~isempty(ud.oldContour)
            delete(ud.oldContour);
        end
        [~,ch] = contour(squeeze(allData.av(ud.currentSlice,:,:)==ann), 1, 'r');
        ud.oldContour = ch;
        set(f, 'UserData', ud);
    end
    
    if ud.showOverlay
        updateOverlay(f, allData, ann, slice_figure, save_location)
    end    
end


function updateStereotaxCoords(currentSlice, pixel, bregma, bregmaText)
atlasRes = 0.010; % mm
ap = -(currentSlice-bregma(1))*atlasRes;
dv = (pixel(1)-bregma(2))*atlasRes;
ml = (pixel(2)-bregma(3))*atlasRes;
set(bregmaText, 'String', sprintf('%.2f AP, %.2f DV, %.2f ML', ap, dv, ml));


function pixel = getPixel(ax)

currPoint = get(ax,'currentpoint');  % The current point w.r.t the axis.

Cx = currPoint(1,1); Cy = currPoint(1,2);
pixel = round([Cy Cx]);


function updateOverlay(f, allData, ann, slice_figure, save_location)
ud = get(f, 'UserData');
if isempty(ud.overlayAx) % first time
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        avo = plotAVoverlay(fliplr(squeeze(allData.av(ud.currentSlice,:,:))'), ann, ud.atlasAx);
    else
        avo = plotAVoverlay(fliplr(ud.im_annotation)', ann, ud.atlasAx);
    end
    ud.overlayAx = avo;
    set(ud.overlayAx, 'HitTest', 'off');
    set(f, 'UserData', ud);
else
    ovIm = get(ud.overlayAx, 'Children');
%     set(ovIm, 'HitTest', 'off');
    
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        thisSlice = squeeze(allData.av(ud.currentSlice,:,:));
    else
        thisSlice = ud.im_annotation;
    end

    set(ovIm, 'CData', flipud(thisSlice));    
    plotAVoverlay(fliplr(thisSlice'), ann, ud.atlasAx, ud.overlayAx);
    set(ovIm, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));

end


function [name, acr, ann] = getPixelAnnotation(allData, pixel, currentSlice)
if pixel(1)>0&&pixel(1)<size(allData.av,2) && pixel(2)>0&&pixel(2)<=size(allData.av,3)
    ann = allData.av(currentSlice,pixel(1),pixel(2));
    name = allData.st.safe_name(ann);
    acr = allData.st.acronym(ann);
else
    ann = []; name = []; acr = [];
end

function updateTitle(ax, name, acr)
if ~isempty(name)
    title(ax, [name{1} ' (' acr{1} ')']);
else
    title(ax, 'not found');
end
