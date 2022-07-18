

function f = allenAtlasBrowser_original(templateVolume, annotationVolume, structureTree)
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
fprintf(1, 'v: switch to Viewing annotations (or switch back) \n');
fprintf(1, 'p: enable/disable mode where clicks are logged as Points \n');

f = figure; 

ud.bregma = allenCCFbregma; 

ud.currentSlice = min(ud.bregma(1), size(templateVolume, 1)); 
ud.oldContour = [];
ud.showContour = false;
ud.showOverlay = false; ud.overlayAx = [];
ud.getPoint = false; ud.pointList = zeros(0,3); ud.pointHands = [];
ud.showAtlas = false;

ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);


ud.im = plotTVslice(squeeze(templateVolume(ud.currentSlice,:,:)));

set(ud.im, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k));

ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;

hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'off');

set(f, 'UserData', ud);

set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, k, allData));
set(f, 'WindowScrollWheelFcn', @(src,evt)updateSlice(f, evt, allData))
set(f, 'WindowButtonMotionFcn',@(f,k)fh_wbmfcn(f, allData)); % Set the motion detector.


function hotkeyFcn(f, keydata, allData)

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
            
end
set(f, 'UserData', ud);


function updateSlice(f, evt, allData)

ud = get(f, 'UserData');
ud.currentSlice = ud.currentSlice+evt.VerticalScrollCount*3;
if ud.currentSlice>size(allData.tv,1); ud.currentSlice = 1; end %wrap around
if ud.currentSlice<1; ud.currentSlice = size(allData.tv,1); end %wrap around
set(f, 'UserData', ud);

if ud.showAtlas
    set(ud.im, 'CData', squeeze(allData.av(ud.currentSlice,:,:)));
else
    set(ud.im, 'CData', squeeze(allData.tv(ud.currentSlice,:,:)));
end

if ~isempty(ud.pointHands)
    % any points that are at this AP position, make them visible
    arrayfun(@(x)pointVis(ud.pointHands(x), ud.currentSlice, ud.pointList(x,:)), 1:numel(ud.pointHands));
end

pixel = getPixel(ud.atlasAx);
updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText);

[name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);

updateTitle(ud.atlasAx, name, acr)

if ud.showOverlay    
    updateOverlay(f, allData, ann);
end

function pointVis(h, sliceInd, coord)
    if sliceInd==coord(3)
        set(h, 'Visible', 'on');
    else
        set(h, 'Visible', 'off'); 
    end


function atlasClickCallback(im, keydata)
f = get(get(im, 'Parent'), 'Parent');
ud = get(f, 'UserData');

if ud.getPoint
    clickX = round(keydata.IntersectionPoint(1));
    clickY = round(keydata.IntersectionPoint(2));
    clickZ = ud.currentSlice;
    
    ud.pointList(end+1, :) = [clickX, clickY, clickZ];
    ud.pointHands(end+1) = plot(ud.atlasAx, clickX, clickY, 'ro');
end
set(f, 'UserData', ud);


function fh_wbmfcn(f, allData)
% WindowButtonMotionFcn for the figure.

ud = get(f, 'UserData');

ax = ud.atlasAx;

pixel = getPixel(ax);

% show bregma coords
updateStereotaxCoords(ud.currentSlice, pixel, ud.bregma, ud.bregmaText);

% get annotation for this pixel
[name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);

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
    avo = plotAVoverlay(fliplr(squeeze(allData.av(ud.currentSlice,:,:))'), ann, ud.atlasAx);
    ud.overlayAx = avo;
    set(ud.overlayAx, 'HitTest', 'off');
    set(f, 'UserData', ud);
else
    ovIm = get(ud.overlayAx, 'Children');
    set(ovIm, 'HitTest', 'off');
    thisSlice = squeeze(allData.av(ud.currentSlice,:,:));
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