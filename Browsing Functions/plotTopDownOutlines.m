

function plotHands = plotTopDownOutlines(st, av, apCoords, mlCoords)

addLabels = false;
addColors = false;
cortexOnly = false;
downsamp = 2; % for speed, doesn't affect graphics much

av = av(1:downsamp:end,1:downsamp:end,1:downsamp:end);

% find first point on the way down
[~,ii] = max(av>1, [], 2);
ii = squeeze(ii);

% now take the annotation of that point
[xx,yy] = meshgrid(1:size(ii,2), 1:size(ii,1));

tdAnn = reshape(av(sub2ind(size(av),yy(:),ii(:),xx(:))), size(av,1), size(av,3));

% tdAnn = tdAnn(1:downsamp:end,1:downsamp:end); % downsample for speed

% figure; imagesc(tdAnn);
% axis image; axis off;
% colormap(cmap); 
% caxis([1 size(cmap,1)]);
% title('first annotation going down');

figure;hold on;
axis image; 
% axis off; 
set(gca, 'YDir','reverse');
uAnn = unique(tdAnn(:));

if cortexOnly 
    isoctxID = st.id(strcmp(st.acronym, 'Isocortex'));
    uAnn = uAnn(cellfun(@(x)contains(x, ['/' num2str(isoctxID) '/']), st.structure_id_path(uAnn)));
else
    uAnn = uAnn(uAnn~=1 & uAnn~=1107 & uAnn~=334 & uAnn~=381 & uAnn~=1143); % exclude root, onl, RSPd, MOB, bic
end

se = strel('disk', 3);
coords = struct();
coordsInd = 1;
plotHands = [];
for u = numel(uAnn):-1:0
    if u==0
        indInSlice = tdAnn>1; 
    else    
        if uAnn(u)==380
            indInSlice = tdAnn==uAnn(u) | tdAnn==381 | tdAnn==1107; %special case for olfactory, which otherwise looks funny
        else
            indInSlice = tdAnn==uAnn(u);
        end
    end
    indInSlice = imerode(imdilate(imdilate(imerode(indInSlice, se),se),se),se); 

    % pad with zeros otherwise contours on the edge do a funny thing
    blankCanvas = zeros(size(indInSlice,1)+2, size(indInSlice,2)+2);
    blankCanvas(2:end-1,2:end-1) = indInSlice; indInSlice = blankCanvas;

    ml = mlCoords(1:downsamp:end); dml = mean(diff(ml)); ml = [ml(1)-dml ml ml(end)+dml];
    ap = apCoords(1:downsamp:end); dap = mean(diff(ap)); ap = [ap(1)-dap ap ap(end)+dap];
    [c,h] = contour(ml,ap, ...
        double(indInSlice), [0.5 0.5]); delete(h);
    ii = 1;
    while ii<size(c,2)
        n = c(2,ii);
        if n>100
            x = c(1,ii+1:ii+n);
            y = c(2,ii+1:ii+n);
            x = smooth(x,50,'loess');
            y = smooth(y,50,'loess');
            x = [x;x(1)];
            y = [y;y(1)];            
            coords(coordsInd).x = x;
            coords(coordsInd).y = y;
            coordsInd = coordsInd+1;
            if addColors
                plotHands(end+1) = plot(x,y, 'LineWidth', 2.0, 'Color', cmap(uAnn(u),:));
            else
                plotHands(end+1) = plot(x,y, 'LineWidth', 2.0, 'Color', 'k');
            end
            
            % using a fill looks nicer but doesn't save to pdf very well
            %fill(x,y, cmap(uAnn(u),:), 'FaceAlpha', 0.2, 'EdgeAlpha', 0);
            
            if addLabels
                mnx = mean(x)-20; mny = mean(y)-8;
                str = st.acronym{uAnn(u)}; 
                if str(end)=='1'; str = str(1:end-1); end
                ah = annotation('textbox', [0.1, 0.1, 0.1, 0.1], 'String', str, 'FontSize', 7, 'LineStyle', 'none');
                set(ah, 'Parent', gca); set(ah, 'Position', [mnx, mny, 0.1, 0.1]);
            end
            
        end
        ii = ii+n+1;
    end
    drawnow;
end