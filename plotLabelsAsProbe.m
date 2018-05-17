

function [borders, fD, fAllen] = plotLabelsAsProbe(m, p, xc, yc, av, st, scaleFactor)

showAllenColors = false;

% these are the query points
yc = yc*scaleFactor;
t = yc/10; % dividing by 10 accounts for the 10um resolution of the atlas
x = m(1)+p(1)*t;
y = m(2)+p(2)*t;
z = m(3)+p(3)*t;

ann = ones(size(t));
cm = zeros(numel(t),3);
allenCmap = allen_ccf_colormap();

% collect annotations along the track
for ind = 1:length(t)
    if x(ind)>0&&x(ind)<=size(av,1)&&...
            y(ind)>0&&y(ind)<=size(av,2)&&...
            z(ind)>0&&z(ind)<=size(av,3)
        ann(ind) = av(ceil(x(ind)), ceil(y(ind)), ceil(z(ind)));
        cm(ind,:) = allenCmap(ann(ind),:);
    end
end

if showAllenColors
    fAllen = figure;
    plotAsProbe([], xc, yc, cm, 16, 40*scaleFactor);    
    %set(fAllen, 'Position', [-1896         -68         116         980]);
    box off;
    ylim([0 3840*scaleFactor]);
    set(gca, 'YDir','reverse');
    xlabel(sprintf('Scale=%.2f',scaleFactor));
else
    fAllen = [];
end

% second version with distinguishable colors
uAnn = unique(ann);
nC = numel(unique(ann(ann>1)));
distinctCmap = distinguishable_colors(nC);
if any(uAnn==1)
    distinctCmap = vertcat([1 1 1], distinctCmap); % always make white be ann==1, outside the brain
end
dc = zeros(max(uAnn),3); 
dc(uAnn,:) = distinctCmap;
cmD = dc(ann,:);

fD = []; % fD = figure;
plotAsProbe([], xc, yc, cmD, 16, 40*scaleFactor)
box off;
ylim([0 3840*scaleFactor]);
xlabel(sprintf('Scale=%.2f',scaleFactor));

% algorithm for finding areas and labels
borders = [0; find(diff(ann)~=0); length(yc)];
midY = zeros(numel(borders)-1,1);
acr = {}; 
for b = 1:length(borders)-1        
    midY(b) = mean(yc(borders(b)+1:min(borders(b+1), length(yc))));
    acr{b} = st.acronym{ann(borders(b)+1)};
    name{b} = st.safe_name{ann(borders(b)+1)};
    annBySegment(b) = ann(borders(b)+1);
%     fprintf(1, '%d\t%d\t%s\t%s\n', round(yc(borders(b)+1)), ...
%         round(yc(min(borders(b+1), length(yc)))), ...
%         acr{b}, name{b});
end
borders = table(yc(borders(2:end-1)), yc(borders(3:end)), acr(2:end)', name(2:end)', annBySegment(2:end)', ...
    'VariableNames', {'upperBorder', 'lowerBorder', 'acronym', 'name', 'avIndex'})
set(gca, 'YTick', midY, 'YTickLabel', acr);
set(gca, 'YDir','reverse');
%set(fD, 'Position', [-1765         -68         146         980]);
