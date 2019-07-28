

function [coords, coordsReg, h] = sliceOutlineWithRegionVec(avSlice, regInd, regColor, ax)

if nargin<2
    regInd = [];
    regColor = [];
end
if nargin<4
    ax = [];
end

if ~isempty(regInd) % specified a region to highlight
    % do this first so it's in the background of the outlines
    c = contourc(double(ismember(avSlice,regInd)), [0.5 0.5]);
    coordsReg = makeSmoothCoords(c);    
    
    if ~isempty(ax)
        for cidx = 1:numel(coordsReg)
            fill(ax, coordsReg(cidx).x,coordsReg(cidx).y, regColor); hold on;                   
        end
    end
else
    coordsReg = [];
end


im = sliceOutlineWithRegion(avSlice);

c = contourc(double(im(:,:,1)), [128 128]);

% now want to smooth the contours and drop the really small ones. 
coords = makeSmoothCoords(c);

if ~isempty(ax)
    clear h
    for cidx = 1:numel(coords)
        h(cidx) = plot(ax, coords(cidx).x,coords(cidx).y, 'k','LineWidth', 2.0); hold on;    
        
        % using a fill looks nicer but doesn't save to pdf very well
        %fill(x,y, cmap(uAnn(u),:), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
    end
    axis(ax,'image');
    axis(ax,'off');
    set(ax, 'YDir', 'reverse'); 
end



return;

function coords = makeSmoothCoords(c)

coords = [];
ii = 1;coordsInd = 1; 
buff = 10;
while ii<size(c,2)
    n = c(2,ii);
    if n>100
        x = c(1,ii+1:ii+n);
        y = c(2,ii+1:ii+n);
        x = [x(end-buff:end) x x(1:buff)]; % buffer makes ends meet
        y = [y(end-buff:end) y y(1:buff)];
        x = smooth(x,25,'loess');
        y = smooth(y,25,'loess');
        x = [x(buff+1:end-buff-1);x(buff+1)];
        y = [y(buff+1:end-buff-1);y(buff+1)];
        coords(coordsInd).x = x;
        coords(coordsInd).y = y;
        coordsInd = coordsInd+1;                                        
        
    end
    ii = ii+n+1;
end