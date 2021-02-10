%% oblique Allen Atlas cross-section script
% Michael Krumin, October 2020

% typically can be found in '..\allenCCF\Browsing Functions\' folder
brainGridData = readNPY('brainGridData.npy'); 

voxelSize = 0.01; % Allen data has voxel size of 10 microns
bp = double(brainGridData); 
bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were
bregma = allenCCFbregma; % [AP, DV, ML]
bregma = bregma([1, 3, 2]);
bp = bsxfun(@minus, bp, bregma) * -1 * voxelSize;

% typically can be found in '..\allenCCF\' folder
tv = readNPY('template_volume_10um.npy'); % grey-scale "background signal intensity"
av = readNPY('annotation_volume_10um_by_index.npy'); % the number at each pixel labels the area, see note below
st = loadStructureTree('structure_tree_safe_2017.csv'); % a table of what all the labels mean

%% Normalizing the coordinates to mm and shifting them to be relative to bregma

% rotation the data to have dimensions in the [AP, ML, DV] order for Matlab plotting
tvPerm = permute(tv, [1, 3, 2]);
avPerm = permute(av, [1, 3, 2]);


[nAP, nML, nDV] = size(tvPerm);
apAxis = (single(1:nAP) - bregma(1)) * -voxelSize;
mlAxis = (single(1:nML) - bregma(2)) * -voxelSize;
dvAxis = (single(1:nDV) - bregma(3)) * -voxelSize;
[ML, AP, DV] = meshgrid(mlAxis, apAxis, dvAxis);


[mlSurf, apSurf] = meshgrid(mlAxis, apAxis);

%% Defining the plane
% where do you want the plane?
use3points = false;
if ~use3points
    % one option is to define a normal to the plane and a single point:
    normal = [0, 1, 1]; % vector normal to the plane (ML, AP, DV)
    intersect = [0, -3, -2]; % a point the plain should go through (ML, AP, DV)
    
else
    % another option would be to define three points the plane should pass through:
    point1 = [-2, -3, -2];
    point2 = [2, -3, -2];
    point3 = [0, -2, -3];
    
    % and then the normal and the single intersect point can be calculated
    intersect = (point1 + point2 + point3)/3;
    normal = cross(point2 - point1, point3 - point1);
    normal = normal/norm(normal);
    [~, ind] = max(abs(normal));
    if normal(ind)<0
        normal = -normal;
    end
end

% plane equation:
% ax + by + cz + d = 0
% normal = [a, b, c]; intersect = [x0, y0, z0]
% solution:
% d = -(ax0 + by0 + cz0);
% z = (-d - ax - by)/c
dvSurf = (normal * intersect' - normal(1) * mlSurf - normal(2) * apSurf) / normal(3);

%% Actual plotting is done here

plotExtras = true; % bregma, intersection points, normal

figure;
hSlice = slice(ML, AP, DV, single(avPerm), mlSurf, apSurf, dvSurf, 'nearest');
hSlice.LineStyle = 'none';

ax = gca;

load allen_ccf_colormap_2017.mat
colormap(ax, cmap); 
caxis(ax, [1 size(cmap,1)]);

hSlice.AlphaData = ones(size(hSlice.XData));
hSlice.AlphaData(hSlice.CData == 1) = 0;
hSlice.AlphaData(isnan(hSlice.CData)) = 0;
hSlice.FaceAlpha = 'flat';

% ax.ZDir = 'reverse';
view(240, 20)

xlabel('ML [mm]');
ylabel('AP [mm]');
zlabel('DV [mm]');
axis equal

% add the mesh of the brain surface
hold on;
plot3(bp(:,2), bp(:, 1), bp(:,3), 'k');

xlim([min(mlAxis), max(mlAxis)]);
ylim([min(apAxis), max(apAxis)]);
zlim([min(dvAxis), max(dvAxis) + 0.5]);

if plotExtras
    % plot bregma
    plot3(0, 0, 0, 'r.', 'MarkerSize', 30);
    % plot the intersect point (or three points defining the plane)
    if use3points
        plot3(point1(1), point1(2), point1(3), 'b.', 'MarkerSize', 30);
        plot3(point2(1), point2(2), point2(3), 'b.', 'MarkerSize', 30);
        plot3(point3(1), point3(2), point3(3), 'b.', 'MarkerSize', 30);
    else
        plot3(intersect(1), intersect(2), intersect(3), 'k.', 'MarkerSize', 30);
    end
    % plot the normal to the plane
    q = quiver3(intersect(1), intersect(2), intersect(3), normal(1), normal(2), normal(3));
    q.MaxHeadSize = 10;
    q.LineWidth = 3;
    q.Color = [0 0 1];
    q.AutoScaleFactor = 3 / norm(normal);
end
