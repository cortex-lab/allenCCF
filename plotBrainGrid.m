

function f = plotBrainGrid(brainGridData, ax)
% function plotBrainGrid([brainGridData], [ax])
% 
% To plot the wire mesh data loaded from brainGridData.npy. 

if isempty(brainGridData)
    mf = mfilename('fullpath');
    brainGridData = readNPY(fullfile(fileparts(mf), 'brainGridData.npy'));
end

bp = double(brainGridData); 
bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were

if isempty(ax)
    f = figure;
    ax = axes('Parent', f);
end
plot3(ax, bp(:,1), bp(:,2), bp(:,3), 'Color', [0 0 0 0.3]);    
set(ax, 'ZDir', 'reverse')
axis(ax, 'equal');
axis(ax, 'vis3d');
axis(ax, 'off');