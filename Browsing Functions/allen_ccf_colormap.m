function cmap = allen_ccf_colormap(varargin)
% colormap for allen annotations "by index". The mat file with the colormap
% must be int he same directory as this function. 
p = mfilename('fullpath');

if nargin>0 && strcmp(varargin{1}, '2017')
    load(fullfile(fileparts(p), 'allen_ccf_colormap_2017.mat'));
else    
    load(fullfile(fileparts(p), 'allen_ccf_colormap.mat'));
end
