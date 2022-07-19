function f = allenAtlasBrowser(f, templateVolume, annotationVolume, structureTree, slice_figure, save_location, save_suffix, plane)
% ------------------------------------------------
% Browser for the allen atlas ccf data in matlab.
% ------------------------------------------------
%
% Inputs templateVolume, annotationVolume, and structureTree are the data describing the atlas.
% The annotation volume should be the "by_index" version
%

% print instructions
display_controls(plane)

% create figure and adjust to user's screen size
% f = figure('Name','Atlas Viewer'); 
figure(f);
try; screen_size = get(0,'ScreenSize'); screen_size = [max(screen_size(3:4)) min(screen_size(3:4))]./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end 
set(f,'Position', [1050*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
movegui(f,'onscreen')

% initialize user data variables held by the figure
ud.plane = plane;
ud.bregma = allenCCFbregma; 
ud.currentSlice = ud.bregma(1); 
ud.currentAngle = zeros(2,1);
ud.scrollMode = 0;

ud.transform_type = 'projective'; %can change to 'affine' or 'pwl'

ud.oldContour = [];
ud.showContour = false;
ud.showOverlay = false; ud.overlayAx = [];
ud.pointList = cell(1,3); ud.pointList{1} = zeros(0,3); 
ud.pointHands = cell(1,3);
ud.probe_view_mode = false;
ud.currentProbe = 0; ud.ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .5 0 .6; 1 .6 0]; 
ud.ProbeColor =  {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','red','purple','orange'};
ud.getPoint_for_transform =false; ud.pointList_for_transform = zeros(0,2); ud.pointHands_for_transform = [];
ud.current_pointList_for_transform = zeros(0,2); ud.curr_slice_num = 1;
ud.clicked = false;
ud.showAtlas = false;
ud.viewColorAtlas = false;
ud.histology_overlay = 0; 
ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);
ud.transform = [];
ud.transformed_slice_figure = [];
ud.slice_shift = 0;
ud.loaded_slice = 0;
ud.slice_at_shift_start = 1;
ud.text = [];

reference_image = squeeze(templateVolume(ud.currentSlice,:,:));
ud.im = plotTVslice(reference_image);
ud.ref_size = size(reference_image);
ud.ref = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_im = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_slice_trans = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.im_annotation = squeeze(annotationVolume(ud.currentSlice,:,:));
ud.atlas_boundaries = zeros(ud.ref_size,'uint16');;
ud.offset_map = zeros(ud.ref_size);
ud.loaded = 0;

% create functions needed to interact with the figure
set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));
ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

ud.angleText = annotation('textbox', [.7 0.95 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;
hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'off');
set(f, 'UserData', ud);
set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, slice_figure, k, allData, save_location, save_suffix));
set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData, slice_figure, save_location))
set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData, slice_figure, save_location)); % Set the motion detector.

% display user controls in the console
function display_controls(plane)
fprintf(1, '\n Controls: \n');
fprintf(1, '--------- \n');
fprintf(1, 'Navigation: \n');
switch plane
    case 'coronal'
        fprintf(1, 'up: scroll through D/V angles (for coronal sections)\n');
        fprintf(1, 'right: scroll through M/L angles  (for coronal sections)\n');
        fprintf(1, 'down: scroll through A/P atlas slices \n');
        fprintf(1, 'left: scroll for Slice Viewer \n');
        
        fprintf(1, 'scroll: move between slices or angles \n');
    case 'sagittal'
        fprintf(1, 'up: scroll through D/V angles (for sagittal sections)\n');
        fprintf(1, 'right: scroll through A/P angles  (for sagittal sections)\n');
        fprintf(1, 'down: scroll through M/L atlas slices \n');
        fprintf(1, 'left: scroll for Slice Viewer \n');

        fprintf(1, 'scroll: move between slices or angles \n');

    case 'transverse'
        fprintf(1, 'up: scroll through M/L angles (for transverse sections)\n');
        fprintf(1, 'right: scroll through A/P angles  (for transverse sections)\n');
        fprintf(1, 'down: scroll through D/V atlas slices \n');
        fprintf(1, 'left: scroll for Slice Viewer \n');

        fprintf(1, 'scroll: move between slices or angles \n');

end

fprintf(1, '\n Registration: \n');
fprintf(1, 't: toggle mode where clicks are logged for transform \n');
fprintf(1, 'h: toggle overlay of current histology slice \n');
fprintf(1, 'p: toggle mode where clicks are logged for probe or switch probes \n');
fprintf(1, 'n: add a new probe \n');
fprintf(1, 'x: save transform and current atlas location \n');
fprintf(1, 'l: load transform for current slice; press again to load probe points \n');
fprintf(1, 's: save current probe \n');
fprintf(1, 'd: delete most recent probe point or transform point \n');
fprintf(1, 'w: enable/disable probe viewer mode for current probe  \n');

fprintf(1, '\n Viewing modes: \n');
fprintf(1, 'o: toggle overlay of current region extent \n');
fprintf(1, 'a: toggle to viewing boundaries \n');
fprintf(1, 'v: toggle to color atlas mode \n');
fprintf(1, 'g: toggle gridlines \n');

fprintf(1, '\n space: display controls \n \n');


% ------------------------
% react to keyboard press
% ------------------------
function hotkeyFcn(f, slice_figure, keydata, allData, save_location, save_suffix)

% retrieve user data from figure
ud = get(f, 'UserData');
ud_slice = get(slice_figure, 'UserData');


key_letter = lower(keydata.Key);
switch key_letter  
% space -- display controls    
    case 'space'
        display_controls(ud.plane)
% o -- toggle showing brain region overlay
    case 'o'
        ud.showOverlay = ~ud.showOverlay;
        if ~ud.showOverlay
            delete(ud.overlayAx); ud.overlayAx = [];
            disp('Overlay OFF');
        elseif ~ud.viewColorAtlas; disp('Overlay on!');
        end
% g -- toggle showing Gridlines    
    case 'g' 
        if ~isfield(ud, 'gridlines') || isempty(ud.gridlines)
            axes(ud.atlasAx); hold on;
            gridY = 100:100:ud.ref_size(1); % assuming the size of the atlas for this for now
            gridX = 70:100:ud.ref_size(2); 
            xl = xlim(); yl = ylim();
            gx = arrayfun(@(x)plot(x*[1 1], yl, 'w'), gridX, 'uni', false);
            gy = arrayfun(@(y)plot(xl, y*[1 1], 'w'), gridY, 'uni', false);
            ud.gridlines = [gx gy];
        elseif strcmp(get(ud.gridlines{1}, 'Visible'), 'on');
            cellfun(@(x)set(x, 'Visible', 'off'), ud.gridlines);
        elseif strcmp(get(ud.gridlines{1}, 'Visible'), 'off');
            cellfun(@(x)set(x, 'Visible', 'on'), ud.gridlines);
        end  
% p -- toggle mode to register clicks as probe points        
    case 'p' 
        ud.probe_view_mode = 0;
        if ud.currentProbe < size(ud.pointList,1)
            ud.currentProbe = ud.currentProbe + 1;
                        
            disp(['probe point mode -- selecting probe ' num2str(ud.currentProbe) ' (' ud.ProbeColor{ud.currentProbe} ')']); 
            ud.getPoint_for_transform = false; 
            
            % show Transformed Slice & Probage Viewer, if not already showing
            if ~ud.slice_at_shift_start; add = 1; else; add = 0; end
                slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift+add}(1:end-4);
                folder_transformations = fullfile(save_location, ['transformations' filesep]);
                try; load([folder_transformations slice_name '_transform_data.mat']);
                    try; figure(ud.transformed_slice_figure); 
                    catch

                        ud.transformed_slice_figure = figure('Name','Transformed Slice & Probe Point Viewer');
                        try; screen_size = get(0,'ScreenSize'); screen_size = [max(screen_size(3:4)) min(screen_size(3:4))]./[2560 1440];
                        catch; screen_size = [1900 1080]./[2560 1440];
                        end

                        set(ud.transformed_slice_figure,'Position', [256*screen_size(1) 12*screen_size(2) 560*screen_size(1) 420*screen_size(2)]);
                        movegui(ud.transformed_slice_figure,'onscreen')        

                        highlight_point = false;
                        transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, highlight_point, [], [], [], [], [], add);
                    end; figure(f);
                end
        else
            ud.currentProbe = 0;
            disp(['probe point mode OFF']);
        end
% w -- toggle mode to visualize probe trajectory           
    case 'w'
        % remove color atlas
        ud.viewColorAtlas = false;
        set(ud.im, 'CData', ud.ref)
        colormap(ud.atlasAx, 'gray'); caxis(ud.atlasAx, [0 400]);
        
        % remove overlay
        ud.showOverlay = 0;
        delete(ud.overlayAx); ud.overlayAx = [];
            
        for probe_plotted = 1:size(ud.pointHands,1)
            set(ud.pointHands{probe_plotted,1}(:), 'Visible', 'off'); 
            set(ud.pointHands{probe_plotted,3}(:), 'Visible', 'off'); 
        end
        if ~ud.currentProbe
            ud.currentProbe = 1;
        end
            ud.probe_view_mode = ~ud.probe_view_mode; 

            if ud.probe_view_mode && ~isempty(ud.pointList{ud.currentProbe,1})
                 % load probe points if none are already up
                if ~size(ud.pointList{1,1},1)
                    probe_points = load(fullfile(save_location, ['probe_points' save_suffix]));  disp('loading probe points')
                    ud.pointList = probe_points.pointList.pointList;
                    ud.pointHands = probe_points.pointList.pointHands;
                for probe = 1:size(ud.pointList,1)
                    for probe_point = 1:size(ud.pointList{probe, 1},1)
                        ud.pointHands{probe, 1}(probe_point) = scatter(ud.atlasAx, ...
                            ud.pointList{probe,1}(probe_point,1), ud.pointList{probe,1}(probe_point,2), 20, 'ro', ...
                        'MarkerFaceColor', ud.ProbeColors(probe, :),'MarkerEdgeColor', ud.ProbeColors(probe, :), ...
                        'LineWidth',3);
                    end
                    if probe ~= ud.currentProbe
                        set(ud.pointHands{probe, 1}(:), 'Visible', 'off'); 
                    end
                end
                    
                end; disp(['activate probe view mode for probe ' num2str(ud.currentProbe)])
            
                % show probe points
                set(ud.pointHands{ud.currentProbe, 1}(:), 'Visible', 'on'); 
                
                curr_probe_handle = ud.pointHands{ud.currentProbe,1};
                set(curr_probe_handle, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));
                
                % switch the order to z y x
                curr_probePoints = ud.pointList{ud.currentProbe,1}(:, [3 2 1]);

                % get line of best fit through points
                % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
                [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
                atlas_res = 0.010;
                disp(['root mean square error to line of ' num2str(round(sqrt(s / size(curr_probePoints,1))*(atlas_res*1000),1)) ' micron'])
                
                min_y = min(ud.pointList{ud.currentProbe,1}(:,2));
                max_y = max(ud.pointList{ud.currentProbe,1}(:,2));
                min_x = m(3) + (min_y - m(2))  * p(3) / p(2);
                max_x = m(3) + (max_y - m(2))  * p(3) / p(2);
                max_z = m(1) + (max_y - m(2))  * p(1) / p(2);
                ud.pointHands{ud.currentProbe, 3} = plot([min_x max_x],[min_y max_y],'color',ud.ProbeColors(ud.currentProbe,:),'linestyle',':','linewidth',1.5);
                set(ud.pointHands{ud.currentProbe,3}, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));
  
                % ensure proper orientation: want 0 at the top of the brain 
                % and positive distance goes down into the brain
                if p(2)<0
                    p = -p;
                end

                                
                % calculate slice angle along probe track -- do so by 
                % constraining either ML angle or DV angle to zero
                position_at_x_center_point = m + (ud.ref_size(2)/2 - m(3)) * p / p(3);
                position_at_y_center_point = m + (ud.ref_size(1)/2 - m(2)) * p / p(2);
                
                angle_DV_if_constraining_ML = round(ud.ref_size(1)/2 / (ud.ref_size(1)/2 - m(2)) * (position_at_y_center_point(1)-m(1)) );
                angle_ML_if_constraining_DV = round(ud.ref_size(2)/2 / (ud.ref_size(2)/2 - m(3)) * (position_at_x_center_point(1)-m(1)) );
                
                if abs(angle_DV_if_constraining_ML) < abs(angle_ML_if_constraining_DV) || abs(angle_DV_if_constraining_ML) < 150
                    ud.currentSlice = round(position_at_y_center_point(1));
                    ud.currentAngle(1) = angle_DV_if_constraining_ML;
                    ud.currentAngle(2) = 0;
                else
                    ud.currentSlice = round(position_at_x_center_point(1));
                    ud.currentAngle(2) = angle_ML_if_constraining_DV;
                    ud.currentAngle(1) = 0;
                end
                    
                % if the angles are too extreme, just show mean slice without angles
                if abs(ud.currentAngle(1))+abs(ud.currentAngle(2)) >= ud.currentSlice || ...
                        abs(ud.currentAngle(1))+abs(ud.currentAngle(2)) >= ud.ref_size(2)-ud.currentSlice
                    ud.currentAngle(1) = 0;
                    ud.currentAngle(2) = 0;
                    ud.currentSlice = round(m(1));
                    disp('probe angle not ideal for viewing in coronal slice -- angles set to 0')
                end
               
                % report estimated probe angle
                AP_angle = round(atand(angle_DV_if_constraining_ML/(ud.ref_size(1)/2)),1);
                ML_angle = round(atand((max_x - min_x)/(max_y - min_y)),1);
                position_at_bregma_depth = [ (m(3) + (0 - m(2)) * p(3) / p(2)) ((m(1) + (0 - m(2)) * p(1) / p(2)))];
                ML_position = round((position_at_bregma_depth(1)-ud.bregma(3))*atlas_res,2);
                AP_position = round((ud.bregma(1) - position_at_bregma_depth(2))*atlas_res,2);
                disp(' ');
                disp('---estimated probe insertion---')
                disp(['entry position at DV = 0: AP = ' num2str(AP_position) ' mm, ' ...
                                                'ML = ' num2str(ML_position) ' mm']) 
                insertion_dist = sqrt((position_at_bregma_depth(1) - max_x)^2+(position_at_bregma_depth(2)-max_z)^2+(0-max_y)^2);
                disp(['insertion distance from the above position: ' num2str(round(insertion_dist*atlas_res,3)) ' mm'])
                direction_AP = {'posterior','anterior'};
                disp([num2str(abs(AP_angle)) ' degrees in the ' direction_AP{(AP_angle<0)+1} ' direction'])
                direction_ML = {'medial','lateral'};
                disp([num2str(abs(ML_angle)) ' degrees in the ' direction_ML{(ML_angle<0&ML_position<0)+1} ' direction']);       
                disp(' ');   
                
                % update slice
                update.VerticalScrollCount = 0; ud.scrollMode = 0; ud.histology_overlay = 0; set(f, 'UserData', ud);
                updateSlice(f, update, allData, slice_figure, save_location); ud = get(f, 'UserData');   
                fill([5 5 250 250],[5 50 50 5],[0 0 0]);

                
                % show Transformed Slice & Probage Viewer, if not already showing
                try; figure(ud.transformed_slice_figure); 
                catch
                    ud.transformed_slice_figure = figure('Name','Transformed Slice & Probe Point Viewer');
                    try; screen_size = get(0,'ScreenSize'); screen_size = screen_size(3:4)./[2560 1440];
                    catch; screen_size = [1900 1080]./[2560 1440];
                    end

                    set(ud.transformed_slice_figure,'Position', [256*screen_size(1) 37*screen_size(2) 560*screen_size(1) 420*screen_size(2)]);
                    movegui(ud.transformed_slice_figure,'onscreen');                  
                    highlight_point = false;
                    transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, highlight_point, [], [], [], [], [],0);
                end; figure(f);                
            
                set(ud.im, 'CData', ud.ref);
                ud.curr_im = ud.ref; set(f, 'UserData', ud);
            else; disp('probe view mode OFF'); end
            
% t -- toggle mode to register clicks as Points   
    case 't' 
        ud.getPoint_for_transform = ~ud.getPoint_for_transform;
        ud.loaded = false;
        
        if ud.getPoint_for_transform
            disp('transform point mode on'); 

            ud.currentProbe = 0;

            % launch transform point mode
            if ~size(ud.current_pointList_for_transform,1)% ) (ud.curr_slice_num ~= (ud.slice_at_shift_start+ud.slice_shift) ||  && ~ud.loaded 
                ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift; %ud_slice.slice_num;
                ud.current_pointList_for_transform = zeros(0,2);
                set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
                num_hist_points = size(ud_slice.pointList,1);
                template_point = 1; template_points_shown = 0;
                updateBoundaries(f,ud, allData); ud = get(f, 'UserData');
            end
        else; disp('transform point mode OFF');    
        end     
% a -- toggle viewing of annotation boundaries  
    case 'a' 
        
        ud.showAtlas = ~ud.showAtlas;
        
        if ud.showAtlas % superimpose boundaries           
            updateBoundaries(f,ud, allData);
        else % return to image
            set(ud.im, 'CData', ud.curr_im)
        end
        
% v -- toggle viewing of Color Atlas        
    case 'v' 
        ud.viewColorAtlas = ~ud.viewColorAtlas;
        ud.histology_overlay = 0;
        
        if ud.viewColorAtlas  
            % remove overlay
            ud.showOverlay = 0;
            ref_mode = false;
            delete(ud.overlayAx); ud.overlayAx = [];            
            set(ud.im, 'CData', ud.im_annotation)
            ud.curr_im = ud.im_annotation;
            cmap = allen_ccf_colormap('2017');
            colormap(ud.atlasAx, cmap); caxis(ud.atlasAx, [1 size(cmap,1)]);   
            fill([5 5 250 250],[5 50 50 5],[1 1 1],'edgecolor',[1 1 1]);
        else           
            set(ud.im, 'CData', ud.ref)
            ud.curr_im = ud.ref;
            colormap(ud.atlasAx, 'gray'); caxis(ud.atlasAx, [0 400]);
            fill([5 5 250 250],[5 50 50 5],[0 0 0]);
        end       
        
        if ud.showAtlas % superimpose boundaries           
            updateBoundaries(f,ud, allData);
        end

    case 'uparrow' % scroll angles along M/L axis
        ud.scrollMode = 1;
        if strcmp(ud.plane,'coronal')
            disp('switch scroll mode -- tilt D/V')
        elseif strcmp(ud.plane,'sagittal')
            disp('switch scroll mode -- tilt D/V')
        elseif strcmp(ud.plane,'transverse')
            disp('switch scroll mode -- tilt M/L')
        end
    case 'rightarrow' % scroll angles along A/P axis
        ud.scrollMode = 2;
        if strcmp(ud.plane,'coronal')
            disp('switch scroll mode -- tilt M/L')
        elseif strcmp(ud.plane,'sagittal')
            disp('switch scroll mode -- tilt A/P')
        elseif strcmp(ud.plane,'transverse')
            disp('switch scroll mode -- tilt A/P')
        end
        
    case 'downarrow' % scroll along A/P axis
        ud.scrollMode = 0;
        if strcmp(ud.plane,'coronal')
            disp('switch scroll mode -- scroll along A/P axis')
        elseif strcmp(ud.plane,'sagittal')
            disp('switch scroll mode -- scroll along M/L axis')
        elseif strcmp(ud.plane,'transverse')
            disp('switch scroll mode -- scroll along D/V axis')
        end
    case 'leftarrow' % scroll along A/P axis
        ud.scrollMode = 3;
        if ~ud.slice_at_shift_start
            ud.slice_at_shift_start = ud_slice.slice_num;
            ud.slice_shift = 0;
        end
        disp('switch scroll mode -- scroll along slice images')     
% h -- toggle viewing of histology / histology overlay
    case 'h'
        disp('  ');
        % remove color atlas
        ud.viewColorAtlas = false;
        set(ud.im, 'CData', ud.ref)
        colormap(ud.atlasAx, 'gray'); caxis(ud.atlasAx, [0 400]);
        % remove overlay
        ud.showOverlay = 0;
        ref_mode = false;
        delete(ud.overlayAx); ud.overlayAx = [];
        % toggle which mode is active
        if ~ud.clicked || ~ud.histology_overlay;
            ud.histology_overlay = ud.histology_overlay + 1 - 3*(ud.histology_overlay==2);
        end
        ud.clicked = false;
        % get clicked points and slice info from the slice figure
        slice_points = ud_slice.pointList;
        slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
        folder_transformations = [save_location filesep 'transformations' filesep];
        if size(ud.current_pointList_for_transform,1)  && size(slice_points,1) && ud.slice_at_shift_start+ud.slice_shift == ud_slice.slice_num && ~ud.probe_view_mode
            key_letter = 'x'; % save transform automatically
        end
        % if in one of the transformation modes, perform the transform
        if (ud.histology_overlay == 1 || ud.histology_overlay == 2) && ...
                ( (size(ud.current_pointList_for_transform,1) && size(slice_points,1)) || ud.loaded)
            try
                if ud.slice_shift > 0
                    ud.curr_slice_trans = imread([folder_transformations slice_name '_transformed.tif']);
                else
                    set(ud.text,'Visible','off');
                    fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');            

                    reference_points = ud.current_pointList_for_transform;
                    slice_points = ud_slice.pointList;

                    current_slice_image = flip(get(ud_slice.im, 'CData'),1);
                    
                    % ** this is where the transform happens, using the
                    % clicked points from the reference and histology images
                    if ~ud.loaded  % use loaded version if 'l' was just pressed 
                        ud.transform = fitgeotrans(slice_points,reference_points,ud.transform_type); %can use 'affine', 'projective', 'polynomial', or 'pwl'
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
            % if wrong number of points clicked
            catch
                ref_mode = true;
                disp(['Unable to transform -- ' num2str(size(ud_slice.pointList,1)) ...
                     ' slice points and ' num2str(size(ud.current_pointList_for_transform,1)) ' reference points']);
                key_letter = 'h'; 
            end
        end
        % if not doing transform, just show reference atlas
        if ud.histology_overlay == 0 || ref_mode
            ud.histology_overlay = 0;
            disp('Reference mode!');
            set(ud.im, 'CData', ud.ref);
            ud.curr_im = ud.ref; set(f, 'UserData', ud);
        end
        if ud.showAtlas
            updateBoundaries(f,ud, allData);
        end
% n -- start marking a new probe        
    case 'n' 
        new_num_probes = size(ud.pointList,1) + 1; 
        ud.getPoint_for_transform = false;
        if new_num_probes <= size(ud.ProbeColors,1)
            disp(['probe ' num2str(new_num_probes) ' added! (' ud.ProbeColor{new_num_probes} ')']);
%         ud.probe_view_mode = 0;
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
        end
 % s -- save probe trajectory and points of each probe per histology image (and associated histology name/number)     
    case 's'
        pointList.pointList = ud.pointList;
        for probe = 1:size(ud.pointList,1)
            set(ud.pointHands{probe,3}, 'ButtonDownFcn', '');  
            set(ud.pointHands{probe,1}, 'ButtonDownFcn', '');        
        end
        pointList.pointHands = ud.pointHands;
        warning('off', 'MATLAB:Figure:FigureSavedToMATFile');
        save(fullfile(save_location, ['probe_points' save_suffix]), 'pointList'); disp('probe points saved');        
% l -- load transform and current slice position and angle        
    case 'l' 
        slice_name = ud_slice.processed_image_names{ud_slice.slice_num}(1:end-4);
        folder_transformations = fullfile(save_location, ['transformations' filesep]);
        ud.clicked = false;
        try
        if ud.loaded_slice+ud.slice_shift ~= ud_slice.slice_num
            
            ud.curr_slice_num = ud_slice.slice_num;
            
            % remove overlay
            ud.showOverlay = 0;
            delete(ud.overlayAx); ud.overlayAx = [];
            
            % load transform data
            transform_data = load(fullfile(folder_transformations, [slice_name '_transform_data.mat']));  
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
                fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');            
            end
            
            disp('transform loaded -- press ''l'' again now to load probe points');
        else % load probe points
            if ~size(ud.pointList{1,1},1)
                probe_points = load(fullfile(save_location, ['probe_points' save_suffix]));  disp('probe points loaded')
                ud.pointList = probe_points.pointList.pointList;
                ud.pointHands = probe_points.pointList.pointHands;
            else
                disp('probe points not loaded -- there are already some current probe points')
            end
        end

        for probe = 1:size(ud.pointList,1)
            % create point plot handles anew
            try; set(ud.pointHands{probe, 1}(:),'Visible','off'); end
            for probe_point = 1:size(ud.pointList{probe, 1},1)
                ud.pointHands{probe, 1}(probe_point) = scatter(ud.atlasAx, ...
                    ud.pointList{probe,1}(probe_point,1), ud.pointList{probe,1}(probe_point,2), 20, 'ro', ...
                'MarkerFaceColor', [.1 .1 .1],'MarkerEdgeColor', ud.ProbeColors(probe, :), ...
                'LineWidth',2);
            
            % set the point plot from the current slice visible
            slice_point_belongs_to = ud.pointHands{probe, 2}(probe_point);
            if slice_point_belongs_to == ud_slice.slice_num
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'on');
            else
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'off'); 
            end  
            
            end
            % turn off best fit line if applicable
            set(ud.pointHands{probe, 3}(:),'Visible','off'); ud.pointHands{probe, 3} = [];
        end
        
        ud.slice_shift = 0;
        catch; 
            disp(['loading failed']); end
% d -- delete current transform or most recent probe point            
    case 'd' 
        if ud.getPoint_for_transform
%             ud.current_pointList_for_transform = zeros(0,2); set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
%             ud.pointHands_for_transform = []; ud_slice.pointList = []; set(slice_figure, 'UserData', ud_slice);
%             disp('current transform erased');
            
            % Try to delete only the most recent point
            ud.current_pointList_for_transform = ud.current_pointList_for_transform(1:end-1,:);
            ud_slice.pointList = ud_slice.pointList(1:end-1,:); 
            set(slice_figure, 'UserData', ud_slice);
            if ud.pointHands_for_transform
                % remove circle for most revent point
                set(ud.pointHands_for_transform(end), 'Visible', 'off');
                ud.pointHands_for_transform = ud.pointHands_for_transform(1:end-1); 
                if ud.pointHands_for_transform
                    % recolor points
                    set(ud.pointHands_for_transform(end), 'color', [0 .9 0]);
                end
            end
            
            disp('transform point deleted');
            
        elseif ud.currentProbe
            ud.pointList{ud.currentProbe,1} = ud.pointList{ud.currentProbe,1}(1:end-1,:);
            ud.pointList{ud.currentProbe,2} = ud.pointList{ud.currentProbe,2}(1:end-1,:);
            ud.pointList{ud.currentProbe,3} = ud.pointList{ud.currentProbe,3}(:,1:end-1);
            set(ud.pointHands{ud.currentProbe, 1}(end), 'Visible', 'off'); 
            ud.pointHands{ud.currentProbe, 1} = ud.pointHands{ud.currentProbe, 1}(1:end-1);
            ud.pointHands{ud.currentProbe, 2} = ud.pointHands{ud.currentProbe, 2}(1:end-1);
            disp('probe point deleted')
        end
end
% x -- save transform and current slice position and angle
if strcmp(key_letter,'x') 
        
        % find or create folder location for transformations
        try
        folder_transformations = fullfile(save_location, ['transformations' filesep]);
        if ~exist(folder_transformations)
            mkdir(folder_transformations)
        end
    
        % find name of slice
        if ud.slice_shift || ud.slice_at_shift_start ~= ud_slice.slice_num
            slice_name = ud_slice.processed_image_names{(ud.slice_at_shift_start+ud.slice_shift)}(1:end-4);
        else
            slice_name = ud_slice.processed_image_names{ud_slice.slice_num}(1:end-4);
        end
        

        if isempty(ud.current_pointList_for_transform)
            ud.transform = [];
        end
        % store transform, if applicable
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
        save(fullfile(folder_transformations, [slice_name '_transform_data.mat']), 'save_transform');
        disp('atlas location saved')
        
        % save transformed histology image
        current_slice_image = imread(fullfile(save_location, [slice_name '.tif']));
%         current_slice_image = flip(get(ud_slice.im, 'CData'),1);
        R = imref2d(size(ud.ref));
        curr_slice_trans = imwarp(current_slice_image, ud.transform, 'OutputView',R);
        imwrite(curr_slice_trans, fullfile(folder_transformations, [slice_name '_transformed.tif']))
        
        disp('transform saved')
        catch
            disp('transform not saved')
        end
end
        
set(f, 'UserData', ud);

% -----------------------------------------
% Update slice (from scrolling or loading)
% -----------------------------------------
function updateSlice(f, evt, allData, slice_figure, save_location)

ud = get(f, 'UserData');


% scroll through slices
if ud.scrollMode==0
    ud.currentSlice = ud.currentSlice+evt.VerticalScrollCount*3;

    if ud.currentSlice>size(allData.tv,1); ud.currentSlice = 1; end %wrap around
    if ud.currentSlice<1; ud.currentSlice = size(allData.tv,1); end %wrap around
    
% scroll through D/V angles        
elseif ud.scrollMode==1 %&&  abs(ud.currentAngle(1)) < 130
  ud.currentAngle(1) = ud.currentAngle(1)+evt.VerticalScrollCount*3;

% scroll through M/L angles
elseif ud.scrollMode==2 %&&  abs(ud.currentAngle(2)) < 130
  ud.currentAngle(2) = ud.currentAngle(2)+evt.VerticalScrollCount*3; 
  
% scroll through slices (left arrow pressed)
elseif ud.scrollMode == 3
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
  ud.showOverlay = 0;
  delete(ud.overlayAx); ud.overlayAx = [];  
  ud_slice = get(slice_figure, 'UserData');
  
  try
    ud.slice_shift = ud.slice_shift-evt.VerticalScrollCount;
    slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
  catch
    ud.slice_shift = ud.slice_shift+evt.VerticalScrollCount;    
    slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
  end
  folder_transformations = fullfile(save_location, ['transformations' filesep]);
    
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
  
    ud.current_pointList_for_transform = zeros(0,2);
    try; load([folder_transformations slice_name '_transform_data.mat']);
       
        ud.clicked = false;
        
        % load transform data
        transform_data = load(fullfile(folder_transformations, [slice_name '_transform_data.mat']));  
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
        ud.loaded = true;
        
        ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift;
        
        ud.histology_overlay = 1;
        
        set(ud.text,'Visible','off');
        fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');
    catch;
        % if no transform, just show reference
        ud.histology_overlay = 0;
        ud.current_pointList_for_transform = zeros(0,2);
        set(ud.im, 'CData', ud.ref);
        ud.curr_im = ud.ref; set(f, 'UserData', ud);   
        set(ud.text,'Visible','off');
        fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift) ' - no transform'],'color','white');        
    end  
        
end  

% update coordinates at the top
pixel = getPixel(ud.atlasAx);
updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText, ud.angleText, ud.currentSlice, ud.currentAngle(1), ud.currentAngle(2), ud.ref_size, ud.plane);
    
% ----------------------------------------
% if no angle, just change reference slice
% ----------------------------------------
if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
    
    reference_slice = squeeze(allData.tv(ud.currentSlice,:,:));
    ud.im_annotation = squeeze(allData.av(ud.currentSlice,:,:));    
 
    if ud.viewColorAtlas
        set(ud.im, 'CData', ud.im_annotation);
    else
        set(ud.im, 'CData', reference_slice);
    end 
    
    
   % update title/overlay with brain region
    [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);
    updateTitle(ud.atlasAx, name, acr)
    if ud.showOverlay    
        updateOverlay(f, allData, ann, slice_figure, save_location);
    end  
    ud.ref = uint8(reference_slice);
    set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
    ud.offset_map = zeros(ud.ref_size); 
    set(f, 'UserData', ud);
    
% ---------------------------
% if angle, angle the atlas
% ---------------------------
else 
  
  image_size = size(squeeze(allData.av(ud.currentSlice,:,:)));
  angle_slice = zeros(image_size);
  
  if ud.currentAngle(1)==0; offset_DV = 0;
  else; offset_DV = -ud.currentAngle(1):sign(ud.currentAngle(1)):ud.currentAngle(1);
  end; start_index_DV = 1; 
 
  
  % loop through AP offsets
  num_DV_iters_add_ind = floor( (image_size(1) - floor( image_size(1) / length(offset_DV))*length(offset_DV)) / 2);
  for curr_DV_iter = 1:length(offset_DV)
      cur_offset_DV = offset_DV(curr_DV_iter);
      if cur_offset_DV == ud.currentAngle(1)
          end_index_DV = image_size(1);
      elseif curr_DV_iter <= num_DV_iters_add_ind  || length(offset_DV - curr_DV_iter) < num_DV_iters_add_ind
          end_index_DV = start_index_DV + floor( image_size(1) / length(offset_DV));
      else
          end_index_DV = start_index_DV + floor( image_size(1) / length(offset_DV)) - 1;
      end
      
       if ud.currentAngle(2)==0;  offset_ML = 0;
       else; offset_ML = -ud.currentAngle(2):sign(ud.currentAngle(2)):ud.currentAngle(2);
       end; start_index_ML = 1;
    % nested: loop through ML offsets
  num_ML_iters_add_ind = floor( (image_size(2) - floor( image_size(2) / length(offset_ML))*length(offset_ML)) / 2);
  for curr_ML_iter = 1:length(offset_ML)
      cur_offset_ML = offset_ML(curr_ML_iter);
      if cur_offset_ML == ud.currentAngle(2)
          end_index_ML = image_size(2);
      elseif curr_ML_iter <= num_ML_iters_add_ind  || length(offset_ML - curr_ML_iter) < num_ML_iters_add_ind
          end_index_ML = start_index_ML + floor( image_size(2) / length(offset_ML));
      else
          end_index_ML = start_index_ML + floor( image_size(2) / length(offset_ML)) - 1;
      end
          
      % update current slice
      try
     angle_slice(start_index_DV:end_index_DV, start_index_ML:end_index_ML) = ...
         squeeze(allData.tv(ud.currentSlice + cur_offset_DV + cur_offset_ML,start_index_DV:end_index_DV,start_index_ML:end_index_ML));
      catch
          disp('')
      end
    
      ud.im_annotation(start_index_DV:end_index_DV,start_index_ML:end_index_ML) = squeeze(allData.av(ud.currentSlice + cur_offset_DV + cur_offset_ML,...
                                                            start_index_DV:end_index_DV,start_index_ML:end_index_ML));
      ud.offset_map(start_index_DV:end_index_DV, start_index_ML:end_index_ML) = cur_offset_DV + cur_offset_ML;
      
      start_index_ML = end_index_ML + 1;
    end
      start_index_DV = end_index_DV + 1;
  end     
  if ud.viewColorAtlas
      set(ud.im, 'CData', ud.im_annotation);
  elseif ~ud.showAtlas  
      set(ud.im, 'CData', angle_slice);
  end  

  ud.ref = uint8(angle_slice);
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
end


% in all cases. update histology overlay
if ud.histology_overlay == 1 || ud.histology_overlay == 2
    updateHistology(f,ud); ud = get(f, 'UserData');
else
    if ud.viewColorAtlas
        ud.curr_im = ud.im_annotation;
    else
        ud.curr_im = ud.ref;
    end    
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
            set(ud.pointHands{probe, 1}(probe_point),'MarkerEdgeColor', color, 'SizeData', 20)
    end
% disp(['mean distance from this slice to probe points is ' num2str(round(mean_distance*10)) ' microns'])
    if mean_distance < 50
        color = abs((ud.ProbeColors(probe,:) * (50 - mean_distance) + [0 0 0] * mean_distance) / 50);
    else
        color = [.33 .33 .33];
    end    
    set(ud.pointHands{probe, 3}, 'Color',color)
    
end
  set(f, 'UserData', ud);

% ---------------------------------------------------------------
% update the image shown if histology is currently being overlaid
% ---------------------------------------------------------------
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

    
% -------------------------------------------------    
% update the position of the region boundary image
% -------------------------------------------------
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
    
% ----------------
% react to clicks
% ----------------
function atlasClickCallback(im, keydata, slice_figure, save_location)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');
ud_slice = get(slice_figure, 'UserData');

% probe view mode
if ud.probe_view_mode && ud.currentProbe
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));   
    if ud.showOverlay 
        clickY = size(ud.ref,1) - clickY;
    end
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
    
    
    % find the probe point closest to this clicked point
    [min_dist, point_ind] = min( sqrt(sum(([clickX clickY clickZ] - ud.pointList{ud.currentProbe}).^2,2)));
    
    % find the slice corresponding to that point
    relevant_slice = ud.pointList{ud.currentProbe,2}(point_ind);
    highlight_point = true;
    ud.transformed_slice_figure = transformed_sliceBrowser(ud.transformed_slice_figure, save_location, f, ...
        highlight_point, relevant_slice, min_dist, clickX, clickY, point_ind,0);
    figure(f);
    
% selecting probe points mode    
elseif ud.currentProbe > 0
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));   
    if ud.showOverlay; clickY = size(ud.ref,1) - clickY; end % overlay inverts Y
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
       
    ud.pointList{ud.currentProbe,1}(end+1, :) = [clickX, clickY, clickZ];
    ud.pointList{ud.currentProbe,2}(end+1, :) = ud.slice_at_shift_start + ud.slice_shift;
    slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start + ud.slice_shift}(1:end-4);
    ud.pointList{ud.currentProbe,3}{end+1} = slice_name;
   
    ud.pointHands{ud.currentProbe, 1}(end+1) = scatter(ud.atlasAx, clickX, clickY, 20, 'ro', ...
                'MarkerFaceColor', [.1 .1 .1],'MarkerEdgeColor', ud.ProbeColors(ud.currentProbe, :), ...
                'LineWidth',2);
                
    ud.pointHands{ud.currentProbe, 2}(end+1) = ud_slice.slice_num + ud.slice_shift;
    
% transform mode    
elseif ud.getPoint_for_transform
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    if ud.showOverlay; clickY = size(ud.ref,1) - clickY; end
    
    if ud.curr_slice_num ~= ud.slice_at_shift_start+ud.slice_shift 
        if ~ud.loaded
            ud.current_pointList_for_transform = zeros(0,2);
            disp('transforming new slice');
        end
    end
    ud.pointList_for_transform(end+1, :) = [clickX, clickY];
    ud.current_pointList_for_transform(end+1, :) = [clickX, clickY];
    set(ud.pointHands_for_transform(:), 'color', [.7 .3 .3]);
    ud.pointHands_for_transform(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro', 'color', [0 .9 0],'LineWidth',2,'markers',4);    
        
    ud.slice_at_shift_start = ud_slice.slice_num;
    ud.slice_shift = 0;
    ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift;
    ud.loaded = 0;
    ud.clicked = true;
end
set(f, 'UserData', ud);

% ------------------------
% react to mouse hovering
% ------------------------
function fh_wbmfcn(f, allData, slice_figure, save_location)
% WindowButtonMotionFcn for the figure.

ud = get(f, 'UserData');
ax = ud.atlasAx;
pixel = getPixel(ax);

%get offset due to angling
if 0<pixel(1) && pixel(1)<=ud.ref_size(1) && 0<pixel(2) && pixel(2)<=ud.ref_size(2)
    offset = ud.offset_map(pixel(1),pixel(2));
else; offset = 0;
end

% show bregma coords
updateStereotaxCoords(ud.currentSlice + offset, pixel, ud.bregma, ud.bregmaText, ud.angleText, ud.currentSlice, ud.currentAngle(1), ud.currentAngle(2), ud.ref_size, ud.plane);

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

% ---------------------------------------------
% update the coordinates shown in the top left
% ---------------------------------------------
function updateStereotaxCoords(currentSlice, pixel, bregma, bregmaText, angleText, slice_num, ap_angle, ml_angle, ref_size, plane)
atlasRes = 0.010; % mm
if strcmp(plane,'coronal')
    ap = -(currentSlice-bregma(1))*atlasRes;
    dv = (pixel(1)-bregma(2))*atlasRes;
    ml = (pixel(2)-bregma(3))*atlasRes;
    set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, ML angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
elseif strcmp(plane,'sagittal')
    ap = -(pixel(2)-bregma(1))*atlasRes;
    dv = (pixel(1)-bregma(2))*atlasRes;
    ml = -(currentSlice-bregma(3))*atlasRes;
    set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, AP angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
elseif strcmp(plane,'transverse')
    ap = -(pixel(2)-bregma(1))*atlasRes;
    dv = (currentSlice-bregma(2))*atlasRes;
    ml = -(pixel(1)-bregma(3))*atlasRes;
    set(angleText, 'String', ['Slice ' num2str(bregma(1) - slice_num) ', ML angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, AP angle ' num2str(round(atand(ml_angle/(ref_size(2)/2)),1)) '^{\circ}']);
end
set(bregmaText, 'String', sprintf('%.2f AP, %.2f DV, %.2f ML', ap, dv, ml));

% ---------------------------------
% update the current mouse location
% ---------------------------------
function pixel = getPixel(ax)

currPoint = get(ax,'currentpoint');  % The current point w.r.t the axis.

Cx = currPoint(1,1); Cy = currPoint(1,2);
pixel = round([Cy Cx]);

% ---------------------------------
% update the overlaid brain region
% ---------------------------------
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

% ---------------------------------
% find the region being hovered on
% ---------------------------------
function [name, acr, ann] = getPixelAnnotation(allData, pixel, currentSlice)
if pixel(1)>0&&pixel(1)<size(allData.av,2) && pixel(2)>0&&pixel(2)<=size(allData.av,3)
    ann = allData.av(currentSlice,pixel(1),pixel(2));
    name = allData.st.safe_name(ann);
    acr = allData.st.acronym(ann);
else
    ann = []; name = []; acr = [];
end

% ---------------------------------
% update the title, showing region
% ---------------------------------
function updateTitle(ax, name, acr)
if ~isempty(name)
    title(ax, [name{1} ' (' acr{1} ')']);
else
    title(ax, 'not found');
end
