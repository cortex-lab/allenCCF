

function gridIn3D(volData, contourHeight, sliceSpacing, origin)

hold on;

% contours of a coronal slice, i.e. at one AP (x) location
% xSlices = round(sliceSpacing/2):sliceSpacing:size(volData,1);
xSlices = origin(1)+(-1000:1000)*sliceSpacing;
xSlices = xSlices(xSlices>0 & xSlices<=size(volData,1));
for x = xSlices
    thisSlice = squeeze(volData(x,:,:));
    
    % compute the contours
    c = contourc(thisSlice, contourHeight*[1 1]);
    cP = parseContours(c);
    
    % plot the contours    
    cellfun(@(yz)plot3(x*ones(1,size(yz,2)), yz(1,:), yz(2,:), 'Color', [0 0 0 0.3]), cP);    
end

% repeat for the other two directions
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

% contours at one ML location, i.e. sagittal slices
% xSlices = round(sliceSpacing/2):sliceSpacing:size(volData,3);
xSlices = origin(3)+(-1000:1000)*sliceSpacing;
xSlices = xSlices(xSlices>0 & xSlices<=size(volData,3));
for x = xSlices
    thisSlice = squeeze(volData(:,:,x));
    
    % compute the contours
    c = contourc(thisSlice, contourHeight*[1 1]);
    cP = parseContours(c);
    
    % plot the contours    
    cellfun(@(yz)plot3(yz(2,:), x*ones(1,size(yz,2)), yz(1,:),  'Color', [0 0 0 0.3]), cP);    
end

xlabel('x'); ylabel('y'); zlabel('z');

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