

function plotBrainGrid(brainGridData)
% function plotBrainGrid([brainGridData])
% 
% To plot the wire mesh data loaded from brainGridData.npy. 

if isempty(brainGridData)
    mf = mfilename('fullpath');
    brainGridData = readNPY(fullfile(fileparts(mf), 'brainGridData.npy'));
end

bp = double(brainGridData); 
bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were

figure;
plot3(bp(:,1), bp(:,2), bp(:,3), 'Color', [0 0 0 0.3]);

set(gca, 'ZDir', 'reverse')
axis equal
axis vis3d
axis off