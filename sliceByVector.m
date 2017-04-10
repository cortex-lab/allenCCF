

function thisSlice = sliceByVector(vol, v, vPerp, refPoint, sampPointsParallel, sampPointsPerp)
% function thisSlice = sliceByVector(vol, v, vPerp, refPoint, sampPointsParallel, sampPointsPerp)
%
% returns a 2d slice of the volume "vol". The slice is along vectors v and
% vPerp, starting at refPoint and extending sampPointsParallel in the
% direction of v, and sampPointsPerp in the direction of vPerp. 
%
% You could imagine a version of this function that uses interp3 rather
% than indexing the matrix directly at closest points. However, this method
% works for volumes that you don't want to interpolate, like the allen
% atlas ccf annotation volume. Could be a parameter to choose the method.

% points along the vector
x = refPoint(1)+v(1)*sampPointsParallel;
y = refPoint(2)+v(2)*sampPointsParallel;
z = refPoint(3)+v(3)*sampPointsParallel;

oX = vPerp(1)*sampPointsPerp; %offset in X direction
oY = vPerp(2)*sampPointsPerp; %offset in Y direction
oZ = vPerp(3)*sampPointsPerp; %offset in Z direction

% slice coordinates
xx = bsxfun(@plus, x, oX');
yy = bsxfun(@plus, y, oY');
zz = bsxfun(@plus, z, oZ');

% worse, older method for getting slice coords
% xx = repmat(x, numel(sampPointsPerp), 1);
% yy = repmat(y, numel(sampPointsPerp), 1);
% zz = repmat(z, numel(sampPointsPerp), 1);
% 
% xx = xx+repmat(oX', 1, numel(sampPointsParallel));
% yy = yy+repmat(oY', 1, numel(sampPointsParallel));
% zz = zz+repmat(oZ', 1, numel(sampPointsParallel));


xx = round(xx); yy = round(yy); zz=round(zz); % we're just going to literally index at these points

xx(xx<1|xx>size(vol,1)) = 1; % the value of tv(1) is 0 so for points outside the volume we'll just index them there
yy(yy<1|yy>size(vol,2)) = 1; % in principle you might want a better solution here though... 
zz(zz<1|zz>size(vol,3)) = 1; % involving NaN's somehow?

% now index the volume at these points
inds = sub2ind(size(vol), xx, yy, zz);

thisSlice = vol(inds);

if nargout==0
    imagesc(thisSlice');
    colormap gray
    set(gca, 'YDir', 'normal');
    axis image; axis off
end
