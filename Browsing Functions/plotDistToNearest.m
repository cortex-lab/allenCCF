function borders = plotDistToNearest(m, p, yc, av, st, scaleFactor, error_length)
% function borders = plotDistToNearest(m, p, yc, av, st, scaleFactor, error_length)
%
% TODO:
% - add also indications of how far to nearest place outside of higher
% level structures. E.g. you might be in DG-mo and close to another area,
% but if that other area is DG-sg then you don't really care, and you'd
% like to also see just the distance to nearest DG.
% - remove scaleFactor and error_length arguments - just calculate from
% root to root, and for a large error_length
% - comments!

% these are the query points along the probe tract
yc = yc*scaleFactor;
t = yc/10; % dividing by 10 accounts for the 10um resolution of the atlas
x = m(1)+p(1)*t;
y = m(2)+p(2)*t;
z = m(3)+p(3)*t;
ortho_plane = null(p);

projection_matrix = ortho_plane * ortho_plane'; %projection matrix onto plane orthogonal to probe tract
scaling_factor = sin(atan2(norm(cross(p,[1 0 1])), dot(p,[1 0 1]))); % sin of angle between x-z plane and probe tract


% collect annotations along the track
annotation_square = ones(size(t,1), error_length*2+1, error_length*2+1);
ann = ones(1,size(t,1));
cm = zeros(numel(t),3);
allenCmap = allen_ccf_colormap();

inds = -error_length:error_length;
dists = (inds.^2+inds'.^2).^(0.5);
dists = dists(:);
[dists, distOrder] = sort(dists);
otherDist = zeros(numel(t),1);
for ind = 1:length(t)
    if x(ind)>0 && x(ind)<=size(av,1) &&...
            y(ind)>0 && y(ind)<=size(av,2) &&...
            z(ind)>0 && z(ind)<=size(av,3)
        
        % add annotation to list of annotations
        ann(ind) = av(ceil(x(ind)), ceil(y(ind)), ceil(z(ind)));
        
        % go in square orthogonal to probe tract and get other regions
        for index1 = -error_length:error_length
            for index2 = -error_length:error_length
                % get index of square index projected onto plane orthogonal to probe
                project_index_onto_ortho_plane = round(projection_matrix * [index1 0 index2]' / scaling_factor);
                % use this index and register its annotation
                xx = ceil(x(ind)+project_index_onto_ortho_plane(1));
                yy = ceil(y(ind)+project_index_onto_ortho_plane(2));
                zz = ceil(z(ind)+project_index_onto_ortho_plane(3));
                if xx>0 && xx<=size(av,1) &&...
                        yy>0 && yy<=size(av,2) &&...
                        zz>0 && zz<=size(av,3)
                    annotation_square(ind,error_length+1+index1,error_length+1+index2) = av(xx,yy,zz);                        
                else
                    annotation_square(ind,error_length+1+index1,error_length+1+index2) = 1;
                end
            end
        end
        
        cur_annotation_square = squeeze(annotation_square(ind,:,:));
        
        currSqVec = cur_annotation_square(:);
        otherInd = find(currSqVec(distOrder)~=currSqVec(distOrder(1)),1);
        if isempty(otherInd)
            otherDist(ind) = dists(end);
        else
            otherDist(ind) = dists(otherInd);
        end
        
        % also add color map
        cm(ind,:) = allenCmap(ann(ind),:);
    end
end


ann(ann==0) = 1;
uAnn = unique(ann);
nC = numel(unique(ann(ann>1)));
distinctCmap = distinguishable_colors(nC);

if any(uAnn==1)
    distinctCmap = vertcat([1 1 1], distinctCmap); % always make white be ann==1, outside the brain
end
dc = zeros(max(uAnn),3); 
dc(uAnn,:) = distinctCmap;
cmD = dc(ann,:)*.8;


% algorithm for finding areas and labels
borders = [0; find(diff(ann)~=0)'; length(yc)];
midY = zeros(numel(borders)-1,1);
acr = {};
for b = 1:length(borders)-1
    
    ycInds = (borders(b):min(borders(b+1)-1, length(yc)))+1;
    theseYC = yc(ycInds);
    
    fill([0 0 otherDist(ycInds)'*10],...
        [max(theseYC) min(theseYC) theseYC], cmD(borders(b)+1,:),...
        'EdgeAlpha', 0);
    hold on;
    
    midY(b) = mean(theseYC);
    acr{b} = st.acronym{ann(borders(b)+1)};
    name{b} = st.safe_name{ann(borders(b)+1)};
    annBySegment(b) = ann(borders(b)+1);
end
% borders = table(yc(borders(2:end-1))', yc(borders(3:end))', acr(2:end)', name(2:end)', annBySegment(2:end)', ...
%     'VariableNames', {'upperBorder', 'lowerBorder', 'acronym', 'name', 'avIndex'})
set(gca, 'YTick', midY, 'YTickLabel', acr);
set(gca, 'YDir','reverse');
xlabel('dist to nearest');
ylim([0 max(yc)])
xlim([0 50]);
box off;





