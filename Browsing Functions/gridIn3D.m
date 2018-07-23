

function contourHands = gridIn3D(volData, contourHeight, sliceSpacing, origin)

hold on;

% contours of a coronal slice, i.e. at one AP (x) location
% xSlices = round(sliceSpacing/2):sliceSpacing:size(volData,1);
xSlices = origin(1)+(-1000:1000)*sliceSpacing;
xSlices = xSlices(xSlices>0 & xSlices<=size(volData,1));
contourHands{1} = {};
for x = xSlices
    thisSlice = squeeze(volData(x,:,:));
    
    % compute the contours
    c = contourc(thisSlice, contourHeight*[1 1]);
    cP = parseContours(c);
    
    % plot the contours    
    contourHands{1}{end+1} = cellfun(@(yz)plot3(x*ones(1,size(yz,2)), yz(1,:), yz(2,:), 'Color', [0 0 0 0.3]), cP, 'uni', false);    
end

% contours at one ML location, i.e. sagittal slices
% xSlices = round(sliceSpacing/2):sliceSpacing:size(volData,3);
xSlices = origin(3)+(-1000:1000)*sliceSpacing;
xSlices = xSlices(xSlices>0 & xSlices<=size(volData,3));
contourHands{2} = {};
for x = xSlices
    thisSlice = squeeze(volData(:,:,x));
    
    % compute the contours
    c = contourc(thisSlice, contourHeight*[1 1]);
    cP = parseContours(c);
    
    % plot the contours    
    contourHands{2}{end+1} = cellfun(@(yz)plot3(yz(2,:), x*ones(1,size(yz,2)), yz(1,:),  'Color', [0 0 0 0.3]), cP, 'uni', false);    
end

% Have elected to leave out the third possible axis - looks worse including
% this
% for x = round(sliceSpacing/2):sliceSpacing:size(volData,2)
%     thisSlice = squeeze(volData(:,x,:));
%     
%     % compute the contours
%     c = contourc(thisSlice, contourHeight*[1 1]);
%     cP = parseContours(c);
%     
%     % plot the contours    
%     cellfun(@(yz)plot3(yz(2,:), yz(1,:), x*ones(1,size(yz,2)),  'Color', [0 0 0 0.3]), cP);    
% end

xlabel('x'); ylabel('y'); zlabel('z');
return

function cout = parseContours(c)
startInd = 1;
cInd = 1;
cout = {};
while startInd<size(c,2)
    nC = c(2,startInd);
    cout{cInd} = c(:,startInd+1:startInd+nC);
    cInd = cInd+1;
    startInd = startInd+nC+1;
end
return

% Want to extract out the data for quick plotting later? use this:
brainGridPlot = nan(0,3); downSamp = 8;
for i = 1:numel(contourHands)
    for j = 1:numel(contourHands{i})
        for k = 1:numel(contourHands{i}{j})
            nD = numel(contourHands{i}{j}{k}.XData(1:downSamp:end)); 
            brainGridPlot(end+1:end+nD,:) = [...
                contourHands{i}{j}{k}.XData(1:downSamp:end)' ...
                contourHands{i}{j}{k}.YData(1:downSamp:end)' ...
                contourHands{i}{j}{k}.ZData(1:downSamp:end)'];
            brainGridPlot(end+1,:) = [NaN NaN NaN]; 
        end
    end
end
