

function overlayAxis = plotAVoverlay(avSlice, indices, tvAxis, varargin)

if ~isempty(varargin)
    overlayAxis = varargin{1};
    im = get(overlayAxis, 'Children');
else
    overlayAxis = axes();
    set(overlayAxis, 'Position', get(tvAxis, 'Position'));

    im = imagesc(avSlice', 'Parent', overlayAxis);
    set(gca, 'YDir','normal')
    axis image
    axis off
    colormap(overlayAxis, allen_ccf_colormap);
    caxis([1 1305]);
end

% make the pixels matching the chosen annotation be semi-transparent, the
% rest completely transparent.
adat =  0.2*ones(size(get(im, 'CData')));
adat(~ismember(avSlice', double(indices))) = 0;
set(im, 'AlphaData', adat)

% if ~isempty(varargin)
%     colormap(tvAxis, 'gray');
% end