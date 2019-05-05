function [f, h] = plotBrainGrid(brainGridData, ax, brain_figure, black_brain)
% function plotBrainGrid([brainGridData], [ax])
% 
% To plot the wire mesh data loaded from brainGridData.npy. 

if nargin<1 || isempty(brainGridData)
    mf = mfilename('fullpath');
    brainGridData = readNPY(fullfile(fileparts(mf), 'brainGridData.npy'));
end

bp = double(brainGridData); 
bp(sum(bp,2)==0,:) = NaN; % when saved to uint16, NaN's become zeros. There aren't any real vertices at (0,0,0) and it shouldn't look much different if there were

if nargin<2||isempty(ax)    
    if nargin<3||isempty(brain_figure)
        brain_figure = figure('Name','Brain View');
    end
    
    ax = axes('Parent', brain_figure);  
end

if nargin<4||~black_brain
    h = plot3(ax, bp(:,1), bp(:,2), bp(:,3), 'Color', [0 0 0 0.3]);
else
    h = plot3(ax, bp(:,1), bp(:,2), bp(:,3), 'Color', [.7 .7 .7 0.3]);
    set(get(ax, 'Parent'),'color','k')
end

set(ax, 'ZDir', 'reverse')
axis(ax, 'equal');
axis(ax, 'vis3d');
axis(ax, 'off');
f = get(ax, 'Parent');
    
   