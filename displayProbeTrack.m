% ------------------------------------------------------------------------
%          Display Probe Track
% ------------------------------------------------------------------------



%% ENTER PARAMETERS AND FILE LOCATION


% file location
processed_images_folder = 'C:\Drive\Histology\for tutorial - sample data\SS096\processed';
probe_save_name_suffix = '_tutorial';

probes_to_analyze = 'all'; %[2];  % either set to 'all' or e.g. [2,3]


% probe parameters
probe_lengths = 5.0; % in mm -- how far into the brain did you go, for all probes or just one num.
active_probe_length = 3.84; % in mm
probe_radius = 100; % in um -- error range queried for confidence metric
show_parent_category = false; %true; % overlay in gray distance between parent regions (takes a while)

probage_past_tip_to_plot = .3; % in mm -- plot this far or to the bottom of the brain, whichever is shorter

scaling_factor = 1.1; % set scaling e.g. based on lining up with ephys, or set to *false* to get scaling automatically from histology

show_region_table = true;


                                        
% directory of reference atlas files
annotation_volume_location = 'C:\Drive\Histology\for tutorial\annotation_volume_10um_by_index.npy';
structure_tree_location = 'C:\Drive\Histology\for tutorial\structure_tree_safe_2017.csv';





% close all

%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE


% load points
probePoints = load(fullfile(processed_images_folder, ['probe_points' probe_save_name_suffix]));
ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
fwireframe = [];

if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end


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
if length(probe_lengths) > 1
    probe_length = probe_lengths(selected_probe);
else
    probe_length = probe_lengths;
end

% get scaling-factor method
if scaling_factor
    use_tip_to_get_reference_probe_length = false;
    reference_probe_length = probe_length * scaling_factor;
    disp(['probe scaling of ' num2str(scaling_factor) ' determined by user input']);    
else
    use_tip_to_get_reference_probe_length = true;
    disp(['getting probe scaling from histology data...']);
end

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
fwireframe.InvertHardcopy = 'off';

% plot probe points
hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3), curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);

% plot brain entry point
plot3(m(1), m(3), m(2), 'k*','linewidth',1)

if use_tip_to_get_reference_probe_length
    % find length of probe in reference atlas space
    [depth, tip_index] = max(curr_probePoints(:,2));
    reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - m).^2)); 
    shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;

    
    disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
    disp(['probe scaling of ' num2str(shrinkage_factor)]);
    disp(' ');
    
    % plot line the length of the probe in reference space
    rpl = round(reference_probe_length_tip);
    
else % use user-defined probe plotting length or scaling factor
    rpl = round(reference_probe_length * 100); 
end

% plot line the length of the active probe sites in reference space
percent_of_tract_with_active_sites = min([active_probe_length / probe_length, 1.0]);
active_site_start = rpl*(1-percent_of_tract_with_active_sites);
apl = round([active_site_start  rpl]);
plot3(m(1)+p(1)*[apl(1) apl(2)], m(3)+p(3)*[apl(1) apl(2)], m(2)+p(2)*[apl(1) apl(2)], ...
    'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);

plot3(m(1)+p(1)*[1 rpl], m(3)+p(3)*[1 rpl], m(2)+p(2)*[1 rpl], ...
    'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');



%% GET AND PLOT BRAIN REGION LABELS ALONG EXTENT OF PROBE

error_length = round(probe_radius / 10); %microns error as first number
% [borders, fD] = plotLabelsAsProbe(m, p, av, st, rpl, error_length, active_site_start*10, probage_past_tip_to_plot); % plots the percent of surrounding area occupied by region along probe
borders_table = plotDistToNearestToTip(m, p, av, st, rpl, error_length, active_site_start, probage_past_tip_to_plot, show_parent_category, show_region_table); % plots confidence score based on distance to nearest region along probe

title(['Probe ' num2str(selected_probe)],'color',ProbeColors(selected_probe,:))


% plot line(s) indicating active site length
plot([0 100], [(active_site_start*10) (active_site_start*10)], 'color',[.1 .1 .1], 'LineStyle',':', 'linewidth',3);
plot([0 100], [(rpl)*10 (rpl)*10], 'color', [.1 .1 .1], 'LineStyle',':', 'linewidth',3);
    
pause(.05)
end