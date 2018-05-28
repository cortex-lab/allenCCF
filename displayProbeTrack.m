% ------------------------------------------------------------------------
%          Display Probe Track
% ------------------------------------------------------------------------



%% ENTER PARAMETERS AND FILE LOCATION

% probe parameters
probe_length = 5.6; % in mm -- how far into the brain did you go
active_probe_length = 3.84; % in mm
probe_radius = 100; % in um

reference_probe_length = false; % in mm -- how far from the surface of the brain the regions are plotted
                                        % set to false to use the deepest probe point selected as the bottom tip of the probe
probage_past_tip_to_plot = 0.0; %0.0; % in mm -- if the above is not set to false, this should be set to zero
                                        
% file location
processed_images_folder = 'C:\\Users\\Experiment\\Desktop\\brain volumes\\slices\\SS096\\processed\\';
probe_save_name_suffix = '';

probes_to_analyze = 1:3; %'all'; %'all'; % either set to 'all' or e.g. [2,3]








%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load points
probePoints = load([processed_images_folder  'probe_points' probe_save_name_suffix]);
ProbeColors = [1 1 1; 1 .75 0; .3 1 1; .7 0 .8; 1 0 0; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .8 .95 .5]; 
fwireframe = [];

if strcmp(probes_to_analyze,'all')
    probes = 1:size(probePoints.pointList.pointList,1);
else
    probes = probes_to_analyze;
end 

% do analysis for each probe
for selected_probe = probes
    
curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);

% if curr_probePoints(1,3)>570 % analyze all on same side (necessary?)
%         curr_probePoints(:,3) = 1140-curr_probePoints(:,3);
% end


% get line of best fit through points
% m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
[m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));


% ensure proper orientation: want 0 at the top of the brain 
% and positive distance goes down into the brain
if p(2)<0
    p = -p;
end

% determine "origin" at top of brain -- step upwards along tract direction
% until tip of brain / past cortex
ann = 10;
isoCtxId = num2str(st.id(strcmp(st.acronym, 'Isocortex')));
gotToCtx = false;
while ~(ann==1 && gotToCtx)
    m = m-p; % step 10um, backwards up the track
    ann = av(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
    if ~isempty(strfind(st.structure_id_path{ann}, isoCtxId))
        % if the track didn't get to cortex yet, keep looking...
        gotToCtx = true;
    end
end

% plot brain grid
fwireframe = plotBrainGrid([], [], fwireframe); hold on; 


% plot probe points
hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);

% plot brain entry point
plot3(m(1), m(3), m(2), 'k*','linewidth',1)

if ~reference_probe_length
    % find length of probe in reference atlas space
    [depth, tip_index] = max(curr_probePoints(:,2));
    reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - m).^2)); 
    shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;

    disp(' ');
    disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
    disp(['shrinkage factor of ' num2str(shrinkage_factor)]);

    % plot line the length of the probe in reference space
    rpl = round(reference_probe_length_tip);
    
    % plot line the length of the active probe sites in reference space
    percent_of_tract_with_active_sites = active_probe_length / probe_length;
    active_site_start = reference_probe_length_tip*(1-percent_of_tract_with_active_sites);
    apl = round([active_site_start  reference_probe_length_tip]);
    plot3(m(1)+p(1)*[apl(1) apl(2)], m(3)+p(3)*[apl(1) apl(2)], m(2)+p(2)*[apl(1) apl(2)], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);

else % use user-defined probe plotting length
    rpl = round(reference_probe_length * 100); active_site_start = false;
end

plot3(m(1)+p(1)*[1 rpl], m(3)+p(3)*[1 rpl], m(2)+p(2)*[1 rpl], ...
    'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');



%% GET AND PLOT BRAIN REGION LABELS ALONG EXTENT OF PROBE

error_length = round(probe_radius / 10); %microns error as first number
% [borders, fD] = plotLabelsAsProbe(m, p, av, st, rpl, error_length, active_site_start*10, probage_past_tip_to_plot); % plots the percent of surrounding area occupied by region along probe
borders = plotDistToNearestToTip(m, p, av, st, rpl, error_length, active_site_start, probage_past_tip_to_plot); % plots confidence score based on distance to nearest region along probe


title(['Brain Regions Along Probe ' num2str(selected_probe)])

% plot line(s) indicating active site length
plot([0 100], [(round(active_site_start)-2)*10 (round(active_site_start)-2)*10], 'color',[.1 .1 .1], 'LineStyle',':', 'linewidth',3);
if ~reference_probe_length
    plot([0 100], [(rpl-2)*10 (rpl-2)*10], 'color', [.1 .1 .1], 'LineStyle',':', 'linewidth',3);
end
end