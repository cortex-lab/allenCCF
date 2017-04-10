function cmap = allen_ccf_colormap()
% colormap for allen annotations "by index". The mat file with the colormap
% must be int he same directory as this function. 
p = mfilename('fullpath');
load(fullfile(fileparts(p), 'allen_ccf_colormap.mat'));