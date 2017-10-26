
function f = allenTilt(tv, av, st)

refPoint = round(size(tv)/2);
ext = 500; 

moveStep = 3; % voxel, =10µm
angStep = 1; % degrees

v = [-1 0 0];
vPerp = [0 0 1];

thisSlice = sliceByVector(tv, v, vPerp, refPoint, -ext:ext,-ext:ext);

f = figure; 
im = imagesc(thisSlice);
axis image
axis off
colormap gray
caxis([0 400]);

ud.im = im; 
ud.atlasAx = gca;
ud.v = v;
ud.vPerp = vPerp;
ud.refPoint = refPoint;
ud.tiltAngle = 0;
ud.rollAngle = 0;
ud.tv = tv;
ud.av = av;
ud.st = st;
ud.pars.moveStep = moveStep;
ud.pars.angStep = angStep;
ud.pars.ext = ext;
ud.pars.mode = 1; 
ud.pars.modeNames = {'normal', 'tilt', 'roll'};
ud.pointList = zeros(0,3);

title(ud.pars.modeNames{ud.pars.mode});

ud.angleText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

ud.areaText = annotation('textbox', [0 0.9 0.4 0.05], ...
    'String', '[area]', 'EdgeColor', 'none', 'Color', 'k');

set(f, 'UserData', ud);

updateAngleText(f);

set(f, 'KeyPressFcn', @(f,k)hotkeyFcn(f, k));
set(f, 'WindowScrollWheelFcn', @(src,evt)scrollFcn(f, evt))
set(f, 'WindowButtonMotionFcn',@(f,k)hoverFcn(f,k)); % Set the motion detector.
set(f, 'ButtonDownFcn', @(f,k)clickFcn(f, k));
set(ud.im, 'HitTest', 'off'); set(ud.atlasAx, 'HitTest', 'off');

function updateSlice(f)
ud = get(f, 'UserData');
ext = ud.pars.ext;
thisSlice = sliceByVector(ud.tv, ud.v, ud.vPerp, ud.refPoint, -ext:ext,-ext:ext);
set(ud.im, 'CData', thisSlice);
drawnow;

function updateAngleText(f)
ud = get(f, 'UserData');
ud.angleText.String = sprintf('tilt = %d, roll = %d, center = [%d, %d, %d]', ...
    ud.tiltAngle, ud.rollAngle, round(ud.refPoint(1)), round(ud.refPoint(2)), round(ud.refPoint(3)));


function moveByStep(f, stepDirection)
ud = get(f, 'UserData');
stepSize = stepDirection*ud.pars.moveStep;
% stepVec = cross(ud.v, ud.vPerp);
% ud.refPoint = ud.refPoint+stepVec*stepSize;
ud.refPoint(2) = ud.refPoint(2)+stepSize;
set(f, 'UserData', ud);
updateSlice(f);
updateAngleText(f)


function moveByTilt(f, stepDirection)
ud = get(f, 'UserData');
stepSize = stepDirection*ud.pars.angStep;
ud.tiltAngle = ud.tiltAngle+stepSize;
[x,y] = pol2cart(ud.tiltAngle/180*pi, 1); % this calculation probably isn't right either
ud.v = [-x y 0];
set(f, 'UserData', ud);
updateSlice(f);
updateAngleText(f)

function moveByRoll(f, stepDirection)
ud = get(f, 'UserData');
stepSize = stepDirection*ud.pars.angStep;
ud.rollAngle = ud.rollAngle+stepSize;
[x,y] = pol2cart(ud.rollAngle/180*pi, 1); % I don't think this calculation is correct
ud.vPerp = [0 y x];
set(f, 'UserData', ud);
updateSlice(f);
updateAngleText(f)

function hotkeyFcn(f, keydata)

ud = get(f, 'UserData');
switch lower(keydata.Key)    
    case 'm'
        ud.pars.mode = ud.pars.mode+1; 
        if ud.pars.mode>numel(ud.pars.modeNames)
            ud.pars.mode = 1;
        end
        title(ud.pars.modeNames{ud.pars.mode});
end
set(f, 'UserData', ud);

function scrollFcn(f, evt)
ud = get(f, 'UserData');
if ud.pars.mode==1
    moveByStep(f, evt.VerticalScrollCount);
elseif ud.pars.mode==2
    moveByTilt(f, evt.VerticalScrollCount);
elseif ud.pars.mode==3
    moveByRoll(f, evt.VerticalScrollCount);    
end

function hoverFcn(f,k)
ud = get(f, 'UserData');
ax = ud.atlasAx;

pixel = getPixel(ax); %(0,0) is upper left, (2*ext, 2*ext) is lower right

[ap, dv, lr] = coordsFromLocation(ud, pixel);

% get annotation for this pixel
[name, acr, ann] = getPixelAnnotation(ud, ap, dv, lr);
set(ud.areaText, 'String', sprintf('%s (%s)', name, acr));

function clickFcn(f,k)
ud = get(f, 'UserData');
ax = ud.atlasAx;

pixel = getPixel(ax);
[ap, dv, lr] = coordsFromLocation(ud, pixel);
[name, acr, ann] = getPixelAnnotation(ud, ap, dv, lr);

ud.pointList(end+1,:) = [ap dv lr];
set(f, 'UserData', ud);
fprintf(1, 'logged point at %d, %d, %d (%s)\n', round(ap), round(dv), round(lr), acr);

function pixel = getPixel(ax)
currPoint = get(ax,'currentpoint');  % The current point w.r.t the axis.
Cx = currPoint(1,1); Cy = currPoint(1,2);
pixel = [Cy Cx];

function [ap, dv, lr] = coordsFromLocation(ud, pixel)
% compute where this is in the original volume
offsetFromRef = pixel-ud.pars.ext*[1 1];
% the interpretation of this offset is that the x offset is the distance
% from the refPoint along v. So if you take the refPoint and add v*xOffs
% then you get the coordinate in the original volume
orig = offsetFromRef(2)*ud.v+offsetFromRef(1)*ud.vPerp+ud.refPoint;
ap = orig(1); 
dv = orig(2);
lr = orig(3);


function [name, acr, ann] = getPixelAnnotation(ud, ap, dv, lr)
if ap>0&&ap<size(ud.av,1) && dv>0&&dv<size(ud.av,2) && lr>0&&lr<size(ud.av,3)
    ann = ud.av(ceil(ap), ceil(dv), ceil(lr));
    name = ud.st.safe_name{ann};
    acr = ud.st.acronym{ann};
else
    ann = []; name = []; acr = [];
end