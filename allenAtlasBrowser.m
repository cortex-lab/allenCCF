

function f = allenAtlasBrowser(templateVolume, annotationVolume, structureTree, slice_figure)
% Browser for the allen atlas ccf data in matlab.
%
% Inputs templateVolume, annotationVolume, and structureTree are the data describing the atlas.
% The annotation volume should be the "by_index" version
%
% Features:
% - Use scroll wheel to go between slices
% - Apply 1mm-spacing grid lines ('g' to enable/disable)
% - View overlay of the region hovered over ('o' to enable/disable)
% - Switch to viewing the atlas annotations ('v' to toggle)
% - Click to register points in 3D ('p' to toggle)
%   - use: 
%     >> ud = get(f, 'UserData'); myPoints = ud.pointList; 
%     to retrieve these points

% todo/future features:
% - switch to sagittal/horizontal slice types


fprintf(1, 'Controls: \n');
fprintf(1, '--------- \n');
fprintf(1, 'scroll: move between slices \n');
fprintf(1, 'g: add/remove Gridlines \n');
fprintf(1, 'o: add/remove Overlay of current region extent \n');
fprintf(1, 'h: add/remove overlay of current histology slice \n');
fprintf(1, 'v: switch to Viewing annotations (or switch back) \n');
fprintf(1, 'p: enable/disable mode where clicks are logged for probe \n');
fprintf(1, 't: enable/disable mode where clicks are logged for transform \n');
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
ud.getPoint = false; ud.pointList = zeros(0,3); ud.pointHands = [];
ud.getPoint_for_transform = false; ud.pointList_for_transform = zeros(0,2); ud.pointHands_for_transform = [];
ud.current_pointList_for_transform = zeros(0,2);
ud.showAtlas = false;
ud.histology_overlay = false; ud.histology = false;
ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);


ud.im = plotTVslice(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.ref = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.im_annotation = zeros(800,1140,'uint16');

set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k));


ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;

hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'off');

set(f, 'UserData', ud);

set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, slice_figure, k, allData));
set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData))
set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData)); % Set the motion detector.


function hotkeyFcn(f, slice_figure, keydata, allData)

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
    case 'c' % toggle Contour
        % contours are disabled for now - they are slow and buggy... 
%         if ud.showContour
%             delete(ud.oldContour); ud.oldContour = [];
%         end
%         ud.showContour = ~ud.showContour;
    case 'o' % toggle Overlay
        ud.showOverlay = ~ud.showOverlay;
        if ~ud.showOverlay
            % weird bug means that point clicking stops working when the
            % overlayAx exists, even though I turned all the hittests off??
            % hack solution is to just delete it every time... :/
%             ovIm = get(ud.overlayAx, 'Children');
%             set(ovIm, 'AlphaData', zeros(size(get(ovIm, 'CData'))));
            delete(ud.overlayAx); ud.overlayAx = [];
        end
    case 'p' % toggle mode to register clicks as Points
        ud.getPoint = ~ud.getPoint;
        if ud.getPoint; disp('probe point mode!'); ud.getPoint_for_transform = false; end
    case 't' % toggle mode to register clicks as Points
        ud.getPoint_for_transform = ~ud.getPoint_for_transform;
        if ud.getPoint_for_transform; disp('transform point mode!'); ud.getPoint = false; end        
    case 'v' % toggle View between template/annotation
        ud.showAtlas = ~ud.showAtlas;
        if ud.showAtlas
            set(ud.im, 'CData', squeeze(allData.av(ud.currentSlice,:,:)));
            cmap = allen_ccf_colormap('2017');
            colormap(ud.atlasAx, cmap); caxis(ud.atlasAx, [1 size(cmap,1)]);            
        else
            set(ud.im, 'CData', squeeze(allData.tv(ud.currentSlice,:,:)));
            colormap(ud.atlasAx, 'gray'); caxis(ud.atlasAx, [0 400]);
        end
        set(f, 'UserData', ud);
    case 'uparrow' % scroll angles along M/L axis
        ud.scrollMode = 1;
        disp('switch scroll mode!')
    case 'rightarrow' % scroll angles along A/P axis
        ud.scrollMode = 2;
        disp('switch scroll mode!')
    case 'downarrow' % scroll along A/P axis
        ud.scrollMode = 0;
        disp('switch scroll mode!')
    case 'h'
        ud.histology_overlay = ~ud.histology_overlay;
        if ud.histology_overlay
            ud_slice = get(slice_figure, 'UserData');
            reference_points = ud.current_pointList_for_transform;
            slice_points = ud_slice.pointList;
            current_slice_image = flip(get(ud_slice.im, 'CData'),1);
            current_reference_image = ud.ref;
            transform = fitgeotrans(slice_points,reference_points,'affine'); %can use 'affine', 'projective', or 'pwl'
            R = imref2d(size(current_reference_image));
            slice_image_transformed = imwarp(current_slice_image, transform, 'OutputView',R);
            image_blend =  imfuse(uint8(current_reference_image*.6), slice_image_transformed(:,:,:),'blend','Scaling','none');
            set(ud.im, 'CData', image_blend);
        else
            set(ud.im, 'CData', ud.ref);
        end
    case 'j'     
        ud.histology = ~ud.histology;
        if ud.histology
          ud_slice = get(slice_figure, 'UserData');
          current_slice_image = flip(get(ud_slice.im, 'CData'),1);
          set(ud.im, 'CData',current_slice_image);
          ud.histology_overlay = ~ud.histology_overlay;
        else
            set(ud.im, 'CData', ud.ref);
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
 

% any points that are at this AP position, make them visible
% if ~isempty(ud.pointHands)
%        arrayfun(@(x)pointVis(ud.pointHands(x), ud.currentSlice, ud.pointList(x,:)), 1:numel(ud.pointHands));
% end


% ---------------------------------
% if no angle, just do normal thing
% ---------------------------------
if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
  if ud.showAtlas
        set(ud.im, 'CData', squeeze(allData.av(ud.currentSlice,:,:)));
  else
        set(ud.im, 'CData', squeeze(allData.tv(ud.currentSlice,:,:)));
  end
 
   % update title/overlay with brain region
    pixel = getPixel(ud.atlasAx);
    updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText);
    [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);
    updateTitle(ud.atlasAx, name, acr)
    if ud.showOverlay    
        updateOverlay(f, allData, ann);
    end
    ud.ref = uint8(get(ud.im, 'CData'));
    ud.current_pointList_for_transform = zeros(0,2);
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
      if ud.showAtlas
          angle_slice(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = ...
              squeeze(allData.av(ud.currentSlice + cur_offset_AP + cur_offset_ML,start_index_AP:end_index_AP,start_index_ML:end_index_ML));
      else
         angle_slice(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = ...
             squeeze(allData.tv(ud.currentSlice + cur_offset_AP + cur_offset_ML,start_index_AP:end_index_AP,start_index_ML:end_index_ML));
      end
    
      ud.im_annotation(start_index_AP:end_index_AP,start_index_ML:end_index_ML) = squeeze(allData.av(ud.currentSlice + cur_offset_AP + cur_offset_ML,...
                                                            start_index_AP:end_index_AP,start_index_ML:end_index_ML));
      ud.offset_map(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = cur_offset_AP + cur_offset_ML;
      
      start_index_ML = end_index_ML + 1;
    end
      start_index_AP = end_index_AP + 1;
  end     
  set(ud.im, 'CData', angle_slice);
  ud.ref = uint8(get(ud.im, 'CData'));
  ud.current_pointList_for_transform = zeros(0,2);
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
  set(f, 'UserData', ud);
end
  
for probe_point = 1:size(ud.pointList,1)
    distance = abs(ud.pointList(probe_point,3) - (ud.currentSlice + ud.offset_map(ud.pointList(probe_point,2),ud.pointList(probe_point,1))));
    if distance < 50
    color = abs(([1 0 0] * (50 - distance) + [.33 .33 .33] * distance) / 50);
    else
    color = [.33 .33 .33];
    end
    set(ud.pointHands(probe_point), 'color', color);
end



% function pointVis(h, sliceInd, coord)
%     if sliceInd==coord(3)
%         set(h, 'Visible', 'on');
%     else
%         set(h, 'Visible', 'off'); 
%     end


function atlasClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');

if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    clickZ = ud.currentSlice + ud.offset_map(clickY,clickX);
    
    ud.pointList(end+1, :) = [clickX, clickY, clickZ];
    ud.pointHands(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro', 'color', 'red','linewidth',2);
    
elseif ud.getPoint_for_transform
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    
    ud.pointList_for_transform(end+1, :) = [clickX, clickY];
    ud.current_pointList_for_transform(end+1, :) = [clickX, clickY];
    ud.pointHands_for_transform(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro', 'color', 'green','linewidth',2);    
    
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
    if ud.showAtlas
        % if showing the atlas annotations, make the cdata 1302 so you get
        % a black overlay (it's the number for the parafloccular sulcus, which 
        % is unused) in current annotations so I stole its entry in the
        % colormap to be black.
        adat = ones(size(get(ovIm, 'CData')));
        set(ovIm, 'CData', 1302*adat);
        adat(thisSlice~=ann) = 0;
        set(ovIm, 'AlphaData', flipud(100*adat));        
    else
        % otherwise show in the color of the atlas
        set(ovIm, 'CData', flipud(thisSlice));    
        plotAVoverlay(fliplr(thisSlice'), ann, ud.atlasAx, ud.overlayAx);
    end
            
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
