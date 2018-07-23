
function im = plotAVslice(thisSlice)
% function plotAVslice(thisSlice)
% 
% useful way to plot a slice you get from sliceByVector, if the slice was
% made from the allen atlas ccf "template volume"

im = imagesc(thisSlice');
% set(gca, 'YDir','normal')
axis image
axis off
colormap(allen_ccf_colormap);
caxis([1 1305]);