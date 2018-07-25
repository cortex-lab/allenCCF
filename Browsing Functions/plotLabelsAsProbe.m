function [borders, fD] = plotLabelsAsProbe(m, p, av, st, rpl, error_length, active_site_start, probage_past_tip_to_plot)

showAllenColors = false;

% these are the query points along the probe tract
yc = 10*[0:(rpl + probage_past_tip_to_plot*100 -1)] - 20;

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
ann2 = ones(size(t,1),2); %top 1 and 2 most common regions without error_length square
ann_percents = zeros(size(t,1), 2); %percent area of each
cm = zeros(numel(t),3);
allenCmap = allen_ccf_colormap();

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
              annotation_square(ind,error_length+1+index1,error_length+1+index2) = av(ceil(x(ind)+project_index_onto_ortho_plane(1)), ...
                                                          ceil(y(ind)+project_index_onto_ortho_plane(2)), ...
                                                          ceil(z(ind)+project_index_onto_ortho_plane(3)));
            end
        end
                                                  
        cur_annotation_square = annotation_square(ind,:,:);
        [count, ann_id] = hist(cur_annotation_square(:), unique(cur_annotation_square(:)));
        [~, count_index] = sort(-count);
        ann2(ind,1:min(length(count_index),2)) = ann_id(count_index(1:min(length(count_index),2)));
        ann_percents(ind, 1:min(length(count_index),2)) = count(count_index(1:1:min(length(count_index),2))) / (error_length*2+1)^2;
        
        % also add color map
        cm(ind,:) = allenCmap(ann(ind),:);
    end
end


% version with distinguishable colors
ann = ann2(:,1);
uAnn = unique(ann);
nC = numel(unique(ann(ann>1)));
distinctCmap = distinguishable_colors(nC);

distinctCmap = vertcat([1 1 1], distinctCmap); % always make white be ann==1, outside the brain
dc = zeros(max(uAnn),3); dc(uAnn,:) = distinctCmap;
cmD = dc(ann,:)*.8;

% put in loop and search colors for second regions
cmD2 = zeros(size(cmD));
for region = unique(ann2(:,2))'
    ann_index = find(ann2(:,1)==region);
    ann2_index = find(ann2(:,2)==region);
    if isempty(ann_index)
        color_to_use = allenCmap(region+1,:); %use allen Cmap if not in ann
        color_to_use = [.8 .8 .8]; % or light grey
    else
    	color_to_use = cmD(ann_index(1),:); %use the same color as in ann
    end
    cmD2(ann2_index,:) = repmat(color_to_use,size(cmD2(ann2_index,:),1),1);
end
    
% plot labelled probe tract
fD = figure('Name','Probe Tract','Position',[1000 100 300 900]);
plotAsProbe([], ann_percents, yc, cmD, cmD2, 10, active_site_start, rpl)
box off;
% ylim([0 3840*scaleFactor]);
xlabel(['% of ' num2str(error_length*10) ' um radius occupation']);
ylim([0 (rpl+probage_past_tip_to_plot*100)*10])


% algorithm for finding areas and labels
borders = [0; find(diff(ann)~=0); length(yc)];
midY = zeros(numel(borders)-1,1);
acr = {}; 
for b = 1:length(borders)-1        
    midY(b) = mean(yc(borders(b)+1:min(borders(b+1), length(yc))));
    acr{b} = st.acronym{ann(borders(b)+1)};
    name{b} = st.safe_name{ann(borders(b)+1)};
    annBySegment(b) = ann(borders(b)+1);
end
borders = table(yc(borders(2:end-1))', yc(borders(3:end))', acr(2:end)', name(2:end)', annBySegment(2:end)', ...
    'VariableNames', {'upperBorder', 'lowerBorder', 'acronym', 'name', 'avIndex'})
set(gca, 'YTick', midY, 'YTickLabel', acr);
set(gca, 'YDir','reverse');

% do the same for secondary areas:




