

function f = allenAtlasBrowser(templateVolume, annotationVolume, structureTree, save_location, save_suffix)
% Browser for the allen atlas ccf data in matlab.
%
% Inputs templateVolume, annotationVolume, and structureTree are the data describing the atlas.
% The annotation volume should be the "by_index" version
%

if nargin<4
    save_location = '';
    save_suffix = '';
end

fprintf(1, 'Controls: \n');
fprintf(1, '--------- \n');
fprintf(1, 'scroll: move between slices \n');
fprintf(1, 'g: add/remove gridlines \n');
fprintf(1, 'o: add/remove overlay of current region extent \n');
fprintf(1, 'a: switch to viewing boundaries (or switch back) \n');
fprintf(1, 'p: enable/disable mode where clicks are logged for probe or switch probes \n');
fprintf(1, 'n: trace a new probe \n');
fprintf(1, 'b: trace a previous probe \n');
fprintf(1, 's: save current probe \n');
fprintf(1, 'w: enable/disable probe viewer mode for current probe  \n');
fprintf(1, 'd: delete most recent probe point \n');

fprintf(1, 'up: scroll through A/P angles \n');
fprintf(1, 'right: scroll through M/L angles \n');
fprintf(1, 'down: scroll through slices \n');

f = figure('Name','Atlas Viewer'); 

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
ud.currentProbe = 0; ud.ProbeColors = [1 1 1; 1 .75 0; .3 1 1; .7 0 .8; 1 0 0; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .8 .95 .5]; 
ud.showAtlas = false;
ud.histology_overlay = 0; 
ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);


ud.im = plotTVslice(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.ref = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_im = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_slice_trans = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.im_annotation = zeros(800,1140,'uint16');
ud.atlas_boundaries = zeros(800,1140,'uint16');;
ud.loaded = 10;

set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, save_location));


ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;

hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'on');

set(f, 'UserData', ud);

set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, k, allData, save_location, save_suffix));
set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData))
set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData)); % Set the motion detector.


function hotkeyFcn(f, keydata, allData, save_location, save_suffix)

ud = get(f, 'UserData');

switch lower(keydata.Key)    
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
            ProbeColor =  {'white','gold','turquoise','purple','red','fern','bubble gum','overcast sky', 'rawhide', 'green apple'};
            disp(['probe point mode -- selecting probe ' num2str(ud.currentProbe) ' (' ProbeColor{ud.currentProbe} ')']); 
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
                    'MarkerFaceColor', ud.ProbeColors(ud.currentProbe, :),'MarkerEdgeColor', [0 0 0], ...
                    'MarkerFaceAlpha',.4,'LineWidth',1.5);
                end
                
                curr_probePoints = ud.pointList{ud.currentProbe,1}(:, [3 2 1]);

                if curr_probePoints(1,3)>570 % analyze all on same side (necessary?)
                        curr_probePoints(:,3) = 1140-curr_probePoints(:,3);
                end

                % get line of best fit through points
                % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
                [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
                
                min_y = min(ud.pointList{ud.currentProbe,1}(:,2));
                max_y = max(ud.pointList{ud.currentProbe,1}(:,2));
                min_x = m(3) + (min_y - m(2))  * p(3) / p(2);
                max_x = m(3) + (max_y - m(2))  * p(3) / p(2);
                ud.pointHands{ud.currentProbe, 3} = plot([min_x max_x],[min_y max_y],'color',ud.ProbeColors(ud.currentProbe,:),'linestyle',':');
                
                % ensure proper orientation: want 0 at the top of the brain 
                % and positive distance goes down into the brain
                if p(2)<0
                    p = -p;
                end

                % use the z value of the mean as the selected slice
                ud.currentSlice = round(m(1)); 

                % calculate slice angle along probe track
                ud.currentAngle(1) = round (400*p(1) / p(2) );
                z_shift_above_probe_center = m(2)*p(1) / p(2);
                ud.currentAngle(2) = round( 570 * (ud.currentAngle(1) - z_shift_above_probe_center) / (m(3) - 570) );
 
                % update slice
                update.VerticalScrollCount = 0; set(f, 'UserData', ud);
                updateSlice(f, update, allData); ud = get(f, 'UserData');                              
            
            else; disp('probe view mode OFF'); end

    case 'a' % toggle View boundaries
        
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
 
    case 'n' % new probe
        new_num_probes = size(ud.pointList,1) + 1; disp(['probe ' num2str(new_num_probes) ' added!']);
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
        
    case 'd' % delete current transform or most recent probe point
        if ud.currentProbe
            ud.pointList{ud.currentProbe,1} = ud.pointList{ud.currentProbe,1}(1:end-1,:);
            ud.pointList{ud.currentProbe,2} = ud.pointList{ud.currentProbe,2}(1:end-1,:);
            ud.pointList{ud.currentProbe,3} = ud.pointList{ud.currentProbe,3}(1:end-1,:);
            set(ud.pointHands{ud.currentProbe, 1}(end), 'Visible', 'off'); 
            ud.pointHands{ud.currentProbe, 1} = ud.pointHands{ud.currentProbe, 1}(1:end-1);
            ud.pointHands{ud.currentProbe, 2} = ud.pointHands{ud.currentProbe, 2}(1:end-1);
            disp('probe point deleted')
        end
end
set(f, 'UserData', ud);


function updateSlice(f, evt, allData)

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
        updateOverlay(f, allData, ann);
    end  
    ud.ref = uint8(get(ud.im, 'CData'));
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

end
ud.curr_im = ud.ref;


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
disp(['mean distance from probe track to points is ' num2str(round(mean_distance*10)) ' microns'])
    if mean_distance < 50
        color = abs((ud.ProbeColors(probe,:) * (50 - mean_distance) + [0 0 0] * mean_distance) / 50);
    else
        color = [.33 .33 .33];
    end    
    set(ud.pointHands{probe, 3}, 'Color',color)
    
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
        image_blend =  uint8( imfuse(ud.curr_im, atlas_boundaries/5,'blend','Scaling','none') )* 2;
        set(ud.im, 'CData', image_blend); 
    end
    
    set(f, 'UserData', ud);
    

function atlasClickCallback(im, keydata, slice_figure, save_location)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');

    
if ud.currentProbe > 0
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
    
    ud.pointList{ud.currentProbe,1}(end+1, :) = [clickX, clickY, clickZ];
    ud.pointList{ud.currentProbe,2}(end+1, :) = 0;
    slice_name = 'no histology image';
    ud.pointList{ud.currentProbe,3}{end+1} = slice_name;
   
    
    ud.pointHands{ud.currentProbe, 1}(end+1) = scatter(ud.atlasAx, clickX, clickY, 20, 'ro', 'MarkerFaceColor', [ .1 .1 .1], ... %ud_atlas_viewer.ProbeColors(probe, :),...
                    'MarkerEdgeColor', ud.ProbeColors(ud.currentProbe, :), 'MarkerFaceAlpha',.4,'LineWidth',1.5);
    ud.pointHands{ud.currentProbe, 2}(end+1) = 0;
    
end
set(f, 'UserData', ud);


function fh_wbmfcn(f, allData)
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
        updateOverlay(f, allData, ann)
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


function updateOverlay(f, allData, ann)
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
    set(ovIm, 'HitTest', 'off');
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        thisSlice = squeeze(allData.av(ud.currentSlice,:,:));
    else
        thisSlice = ud.im_annotation;
    end

    set(ovIm, 'CData', flipud(thisSlice));    
    plotAVoverlay(fliplr(thisSlice'), ann, ud.atlasAx, ud.overlayAx);

            
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
