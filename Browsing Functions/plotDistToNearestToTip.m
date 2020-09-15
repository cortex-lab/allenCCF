function borders = plotDistToNearestToTip(m, p, av, st, rpl, error_length, active_site_start, probage_past_tip_to_plot, show_parent_category, show_region_table, plane)


% these are the query points along the probe tract
yc = 10*[0:(rpl + probage_past_tip_to_plot*100 )] - 0;
t = yc/10; % dividing by 10 accounts for the 10um resolution of the atlas
% select the plane for the viewer
if strcmp(plane,'coronal')
    m = m;
    p = p;
elseif strcmp(plane,'sagittal')
    m = m([3 2 1]);
    p = p([3 2 1]);
elseif strcmp(plane,'transverse')
    m = m([2 3 1]);
    p = p([2 3 1]);
end
% get the x y and z coordinates of the track
x = m(1)+p(1)*t;
y = m(2)+p(2)*t;
z = m(3)+p(3)*t;


ortho_plane = null(p);

projection_matrix = ortho_plane * ortho_plane'; % projection matrix onto plane orthogonal to probe tract
scaling_factor = sin(atan2(norm(cross(p,[1 0 1])), dot(p,[1 0 1]))); % sin of angle between x-z plane and probe tract

% plot labelled probe tract
fD = figure('Name','Probe Tract');
try; screen_size = get(0,'ScreenSize'); screen_size = [max(screen_size(3:4)) min(screen_size(3:4))]./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end
    
set(fD,'Position', [600*screen_size(1) 50*screen_size(2) 200*screen_size(1) 1275*screen_size(2)])
movegui(fD,'onscreen')
box off;

for ann_type = 1:1+show_parent_category % loop between specific and parent region annotations
    
    if ann_type == 2
        disp('calculating confidences for the parent regions...');
    end
% collect annotations along the track
annotation_square = ones(size(t,2), error_length*2+1, error_length*2+1);
ann = ones(1,size(t,2));
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
        if ann_type == 1
            ann(ind) = av(ceil(x(ind)), ceil(y(ind)), ceil(z(ind)));
            cm(ind,:) = allenCmap(ann(ind),:);    % also add color map?         
        elseif ann_type == 2
            id_path = st(av(ceil(x(ind)), ceil(y(ind)), ceil(z(ind))),:).structure_id_path; % get annotation for parent region
            id_path_array = regexp(id_path, '/*/', 'split');
            if size(id_path_array{1},2) > 4
                ann(ind) = str2num(id_path_array{1}{size(id_path_array{1},2)-3});
            end            
            
        end
        
        % go in square orthogonal to probe tract and get other regions
        try
        for index1 = -error_length:error_length
            for index2 = -error_length:error_length
              % get index of square index projected onto plane orthogonal to probe  
              project_index_onto_ortho_plane = round(projection_matrix * [index1 0 index2]' / scaling_factor);
              % use this index and register its annotation
              if ann_type == 1
                 annotation_square(ind,error_length+1+index1,error_length+1+index2) = av(ceil(x(ind)+project_index_onto_ortho_plane(1)), ...
                                                          ceil(y(ind)+project_index_onto_ortho_plane(2)), ...
                                                          ceil(z(ind)+project_index_onto_ortho_plane(3)));
                % get parent annotations
              elseif ann_type == 2
                 id_path = st(av(ceil(x(ind)+project_index_onto_ortho_plane(1)), ...
                                                          ceil(y(ind)+project_index_onto_ortho_plane(2)), ...
                                                          ceil(z(ind)+project_index_onto_ortho_plane(3))),14).structure_id_path;   
                  
                 id_path_array = regexp(id_path, '/*/', 'split');
                 if size(id_path_array{1},2) > 4
                     annotation_square(ind,error_length+1+index1,error_length+1+index2) = str2num(id_path_array{1}{size(id_path_array{1},2)-3});
                 end
                 
              end
              
           end
        end
        catch; disp('you''re way out of the brain');end
        
        % square of annotations orthogonal to point
        cur_annotation_square = squeeze(annotation_square(ind,:,:));
        currSqVec = cur_annotation_square(:);
        
        % find where the annotation isn't equal to that at the curent point
        otherInd = find(currSqVec(distOrder)~=currSqVec(distOrder(1)),1); 
        if isempty(otherInd)
            otherDist(ind) = dists(end);
        else
            otherDist(ind) = dists(otherInd);
        end
    end
end

% get the color map
if (ann_type == 1 && ~show_parent_category) || (ann_type == 2 && show_parent_category)
    uAnn = unique(ann);
    nC = numel(unique(ann(ann>1)));
    distinctCmap = flip(distinguishable_colors(nC+1));

    if any(uAnn==1) || any(uAnn==0)
        distinctCmap = vertcat([1 1 1], distinctCmap); % always make white be ann==1, outside the brain
        uAnn(uAnn==0)=1; ann(ann==0)=1;
    end
    dc = zeros(max(uAnn),3); dc(uAnn,:) = distinctCmap(1:end-1,:);
    cmD = dc(ann,:);
    cmD = cmD*.3*(1+show_parent_category*.5);
else
    cmD = ones(length(ann),3) * .05;
end

% algorithm for finding areas and labels
region_borders = [0; find(diff(ann)~=0)'; length(yc)]';

if ann_type == 1
    midY = zeros(numel(region_borders)-1 - logical( sum(region_borders==round(rpl)) ) - logical(sum(region_borders==round(active_site_start))) ,1);
    acr = {}; 
    acr_for_table = {};
    shift_ind = 0;
end


% add active site start and rpl as pseudo-borders
% if active_site_start
    if ~sum(region_borders==round(active_site_start))
        region_borders = [region_borders(1:max(find(region_borders<active_site_start)) ) round(active_site_start)  ...
            region_borders(max(find(region_borders<active_site_start))+1:end)];
        active_site_start_is_boundary = 0;
    else
        active_site_start_is_boundary = 1;
    end
    if ~sum(region_borders==round(rpl))
        region_borders = [region_borders(1:find(region_borders>rpl,1)-1) round(rpl) region_borders(find(region_borders>rpl,1):end)];
        probe_tip_is_boundary = 0;
    else
        probe_tip_is_boundary = 1;
    end
% else; 
    
    
% end
borders = region_borders;



% plot a 2D shape at each border
for b = 1:length(borders)-1    
     
    ycInds = (borders(b):min(borders(b+1)-1, length(yc)))+1;
    theseYC = yc(ycInds);
    

    if max(theseYC) < active_site_start*10  || max(theseYC) > rpl*10
        cur_alpha = .5 / ann_type;
    else
        cur_alpha = .7 / ann_type;
    end   
    
    if size(theseYC,2) < 3
    fill([0 0 otherDist(ycInds(1))*10 otherDist(ycInds(end))'*10 otherDist(ycInds(end))'*10 0],...
        [max(theseYC)+5 min(theseYC)-5  min(theseYC)-5 max(theseYC)+5  max(theseYC)+5 max(theseYC)+5], cmD(borders(b)+1,:),...
        'EdgeAlpha', 0, 'FaceAlpha', cur_alpha); 
    else

    fill([0 0 otherDist(ycInds(1))*10 otherDist(ycInds)'*10 otherDist(ycInds(end))*10 0],...
        [max(theseYC)+5 min(theseYC)-5  min(theseYC)-5 theseYC  max(theseYC)+5 max(theseYC)+5], cmD(borders(b)+1,:),...
        'EdgeAlpha', 0, 'FaceAlpha', cur_alpha); 
    end
    hold on;

    if ann_type==1
        if (borders(b+1) == round(active_site_start) && ~active_site_start_is_boundary) || ( (borders(b) == round(rpl)) && ~probe_tip_is_boundary)
     	  shift_ind = shift_ind + 1;
        else
            midY(b - shift_ind) = mean(theseYC);
            if max(otherDist(ycInds)) == 1 % don't show label if region is one pixel away from another region for its full extent (plot gets crowded)
                acr{b - shift_ind} = '';
            else
                acr{b - shift_ind} = st.acronym{ann(borders(b)+1)};    
            end            
%             name{b - shift_ind} = st.safe_name{ann(borders(b)+1)};
%             annBySegment(b - shift_ind) = ann(borders(b)+1);     
        end
        acr_for_table{b} = st.acronym{ann(borders(b)+1)};
        name{b} = st.safe_name{ann(borders(b)+1)};
        annBySegment(b) = ann(borders(b)+1); 
    end
end
yc(borders(2)) = 0;
if show_region_table && ann_type==1
    borders_table = table(yc(borders(2:end-1))', yc(borders(3:end))', acr_for_table(2:end)', name(2:end)', annBySegment(2:end)', ...
     'VariableNames', {'upperBorder', 'lowerBorder', 'acronym', 'name', 'avIndex'})
end


end


 
if strcmp(acr{end}, 'root')
    acr{end} = 'end';
end

yyaxis left
set(gca, 'YTick', midY, 'YTickLabel', acr);
set(gca, 'YDir','reverse');
set(gca,'Color',[1 1 1]*.85);
xlabel('dist to nearest (\mum)','color','k');
set(gca,'fontsize',10)

% ylim([0 yc(borders(end-1))])
ylim([1 yc(end)+1])
xlim([0 error_length*10]);
set(fD,'Color',[1 1 1]*.85)
fD.InvertHardcopy = 'off';

yyaxis right
if active_site_start> 0
    set(gca, 'YTick', [1 active_site_start*10 rpl*10 yc(end)+(1-probage_past_tip_to_plot)*100], 'YTickLabel', {'0 ' [num2str(active_site_start*10) ' '] [num2str(rpl*10) ' '] [num2str(yc(end)+(1-probage_past_tip_to_plot)*100) ' ']});
else
    set(gca, 'YTick', [1 rpl*10 yc(end)+(1-probage_past_tip_to_plot)*100], 'YTickLabel', {'0' [num2str(rpl*10) ' '] [num2str(yc(end)+(1-probage_past_tip_to_plot)*100) ' ']});
end
set(gca, 'YDir','reverse');
set(gca,'YColor',[1 1 1]*.2)
ylim([1 yc(end)+1])

% plot line(s) indicating active site length
plot([0 100], [(active_site_start*10) (active_site_start*10)], 'color',[.1 .1 .1], 'LineStyle',':', 'linewidth',3);
plot([0 100], [(rpl)*10 (rpl)*10], 'color', [.1 .1 .1], 'LineStyle',':', 'linewidth',3);

box off;

