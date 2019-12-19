

function addAllenCtxOutlines(bregma, lambda, lineColor, pixSize)
% to get bregma/lambda, do:
% >> figure; imagesc(myImg); axis image; 
% then use the cursor to click on bregma, then swap the coords (if cursor
% says X=xx and Y=yy, then bregma is [yy, xx]).
% or leave blank to get prompted with ginput. 
% lambda just needs to be on the midline posterior to bregma.
if isempty(bregma)
    fprintf(1, 'click bregma\n');
    [x,y] = ginput(1);
    bregma = [y,x]
    fprintf(1, 'click lambda\n');
    [x,y] = ginput(1);
    lambda = [y,x]
end

apDir = lambda-bregma; apDir = apDir./norm(apDir);

if nargin<4
    pixSize = 0.0217; % mm/pix. This is for PCO edge 5.5 with 0.6x mag (as kilotrode)
end

ccfbregma = allenCCFbregma()/100/pixSize;


load(fullfile(fileparts(mfilename('fullpath')), 'ctxOutlines.mat'));

hold on; 
for q = 1:numel(coords) % coords is from ctxOutlines.mat
    
    % these are in 10um voxel coordinates, so first convert to mm, then to
    % pixels 
    cx = coords(q).x/100/pixSize;
    cy = coords(q).y/100/pixSize;
    
    % to do this transformation, first subtract bregma to zero, then
    % rotate, then add back the other bregma    
    cx = cx-ccfbregma(3); cy = cy-ccfbregma(1);
    
    T = affineMat.rot(-atan(apDir(2)/apDir(1)));
    
    newc = T*[cx cy ones(size(cx))]';
    cx = newc(1,:)'; cy = newc(2,:)';
    
    cx = cx+bregma(2); cy = cy+bregma(1);
    
    plot(cx,cy, 'LineWidth', 1.0, 'Color', lineColor);
    hold on;
end