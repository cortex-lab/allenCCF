function [probes, fwireframe, fwireframe2, fig_probes, T_probes, Tapdvml_contacts, T_borders] ...
    = plot_and_compute_probe_positions_from2(av, st, subject_id, plane,...
    processed_images_folder, probe_save_name_suffix,probes_to_analyze, probe_lengths,...
    active_probe_length, probe_radius)
% Create 3D plot for probe positions, 2D plots for brain structures along
% each probe, and return tables for structure names and coordinates
%
%
% SYNTAX
% [] = plot_and_compute_probe_positions(fwireframe, black_brain, probes, subject_id, plane,...
%    use_tip_to_get_reference_probe_length, probe_insertion_direction)
% [D] = func1(A,B)
% [D] = func1(____,'Param',value)
%
% longer description may come here
%
% INPUT ARGUMENTS
% av          uint16
%             3 dimensional array of integers representing the different
%             brain regions in the model space
%
% st          table
%
% subject_id  char
%             Animal name
%
% plane       'coronal' | 'sagittal' | 'transverse'
%
% processed_images_folder
%             char
%
% probe_save_name_suffix
%             char
%             Needed to access stored data
%
% probes_to_analyze
%             'all' | a row vector of positive integers specifying probes
%             

% B           0 (default) | non-negative integers
%             (Optional) Description about B comes here.
% 
%
% OPTIONAL PARAMETER/VALUE PAIRS
% 'C'         'on' (default) | 'off'
%             (Optional) Description about 'C' comes here.
%
%
% OUTPUT ARGUMENTS
% D           cell array 
%             Description about D comes here.
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 09-Aug-2023 11:12:42
%
% See also
% doc



% load the reference brain annotations
% if ~exist('av','var') || ~exist('st','var')
%     disp('loading reference atlas...')
%     av = readNPY(annotation_volume_location);
%     st = loadStructureTree(structure_tree_location);
% end

arguments

    av (:,:,:) uint16
    st (:, 22) table
    subject_id (1, :) char
    plane (1, :) char {mustBeMember(plane,{'coronal','sagittal','transverse'})}
    processed_images_folder (1, :) char {mustBeFolder}
    probe_save_name_suffix (1, :) char
    probes_to_analyze % 'all' | A vector of integer indices for probes
    probe_lengths {mustBeVector}
    active_probe_length (1,1)
    probe_radius (1,1)

end


scaling_factor = false;
show_region_table = true;
black_brain = true;
probe_insertion_direction = 'down';
show_parent_category = false; 
distance_past_tip_to_plot = .5;



% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
end

% load probe points
probePoints = load(fullfile(processed_images_folder, ['probe_points' probe_save_name_suffix]));
ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','red','purple','orange'};
fwireframe = [];

% % scale active_probe_length appropriately
% active_probe_length = active_probe_length*100;

% determine which probes to analyze
if strcmp(probes_to_analyze,'all')
    probes = 1:size(probePoints.pointList.pointList,1);
else
    probes = probes_to_analyze;
end 



%%

% create a new figure with wireframe
fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on; 
fwireframe.InvertHardcopy = 'off';

fwireframe2 = plotBrainGrid([], [], fwireframe, black_brain); %TODO what is this one for?????????
hold on;
fwireframe2.InvertHardcopy = 'off';

fig_probes = gobjects(1,length(probes));
borders_table = cell(1,length(probes)) %TODO
t_size = [size(probePoints.pointList.pointList,1),1];
T_probes = table(repmat(string(subject_id), t_size),...
    strings(t_size), ...
    NaN(t_size),...
    strings(t_size), ...
    NaN(t_size),...
    NaN(t_size),...
    NaN(t_size),...
    NaN(t_size),...
    NaN(t_size),...
    NaN(t_size),...
    strings(t_size), ...    
    'VariableNames',{'subject_id','session_id','probe_id','probe_AB',...
        'depth_reported_um','depth_from_ACCF_um','scaling', ...
        'scaling2',...
        'tip_active_um', 'top_active_um', 'probe_note'});
T_probes.probe_id(probes') = probes';
T_probes.probe_lengths(probes) = probe_lengths(probes)' * 1000;


c_Tapdvml_contacts = cell(length(probes),1);

for selected_probe = probes

    % get the probe points for the currently analyzed probe
    if strcmp(plane,'coronal')
        curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [3 2 1]);
    elseif strcmp(plane,'sagittal')
        curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 2 3]);
    elseif strcmp(plane,'transverse')
        curr_probePoints = probePoints.pointList.pointList{selected_probe,1}(:, [1 3 2]);
    end

    assert(size(curr_probePoints,1) ==2, 'Only accepts 2 points (the tip and surface) per probe')

    % get user-defined probe length from experiment
    if length(probe_lengths) > 1
        probe_length = probe_lengths(selected_probe);
    else
        probe_length = probe_lengths;
    end

    % get the scaling-factor method to use
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
    [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3)); %TODO
    if isnan(m(1))
        disp(['no points found for probe ' num2str(selected_probe)])
        continue
    end

    % ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
    if p(2)<0
        p = -p;
    end


    % determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
    ann = 10;
    out_of_brain = false;
    while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
        m = m-p; % step 10um, backwards up the track
        ann = av(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
        if strcmp(st.safe_name(ann), 'root')
            % make sure this isn't just a 'root' area within the brain
            m_further_up = m - p*20; % is there more brain 200 microns up along the track?
            ann_further_up = av(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
            if strcmp(st.safe_name(ann_further_up), 'root')
                out_of_brain = true;
            end
        end
    end
    %NOTE m is now the highest point

    %% focus on wireframe plot
    figure(fwireframe);

    % plot probe points
    hp = plot3(curr_probePoints(:,1), curr_probePoints(:,3),...
        curr_probePoints(:,2), '.','linewidth',2, 'color',[ProbeColors(selected_probe,:) .2],'markers',10);

    % plot brain entry point
    plot3(m(1), m(3), m(2), 'r*','linewidth',1)

    % use the deepest clicked point as the tip of the probe, if no scaling provided (scaling_factor = false)
    if use_tip_to_get_reference_probe_length
        % find length of probe in reference atlas space
        if strcmp(probe_insertion_direction, 'down')
            [depth, tip_index] = max(curr_probePoints(:,2));
        elseif strcmp(probe_insertion_direction, 'up')
            [depth, tip_index] = min(curr_probePoints(:,2));
        end
        reference_probe_length_tip = sqrt(sum((curr_probePoints(tip_index,:) - m).^2));

        % and the corresponding scaling factor
        shrinkage_factor = (reference_probe_length_tip / 100) / probe_length;

        % display the scaling
        disp(['probe length of ' num2str(reference_probe_length_tip/100) ...
            ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
        disp(['probe scaling of ' num2str(shrinkage_factor)]); disp(' ');

        T_probes.depth_from_ACCF_um(selected_probe) = reference_probe_length_tip/100 * 1000; 
        T_probes.scaling(selected_probe) = shrinkage_factor; 

        % plot line the length of the probe in reference space
        probe_length_histo = round(reference_probe_length_tip);
        % probe_length_histo is the distance from the tip to the surface
        % measured in Allen CCF space.

        % if scaling_factor is user-defined as some number, use it to plot the length of the probe
    else
        probe_length_histo = round(reference_probe_length * 100);
    end

    % find the percent of the probe occupied by electrodes
    % percent_of_tract_with_active_sites = min([active_probe_length*100 / (probe_length*100), 1.0]); %TODO here probe_length is used, but this should be ignored
    % active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);

    %NOTE ignore the probe_length. Use probe_length_histo instead.
    percent_of_tract_with_active_sites = min([active_probe_length*100 / probe_length_histo, 1.0]); 
    active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);
    active_probe_position = round([active_site_start  probe_length_histo]); %TODO use this?

    % plot line the length of the active probe sites in reference space
    plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], ...
        m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], ...
        m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);
    % plot line the length of the entire probe in reference space
    plot3(m(1)+p(1)*[1 probe_length_histo], m(3)+p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');


    %% ----------------------------------------------------------------
    % Get and plot brain region labels along the extent of each probe
    % ----------------------------------------------------------------

    % convert error radius into mm
    error_length = round(probe_radius / 10);

    % find and regions the probe goes through, confidence in those regions, and plot them
    borders_table{selected_probe} = plotDistToNearestToTip(m, p, av_plot, st, probe_length_histo, ...
        error_length, active_site_start, distance_past_tip_to_plot, ...
        show_parent_category, show_region_table, plane); % plots confidence score based on distance to nearest region along probe
    title(['Probe ' num2str(selected_probe)],'color',ProbeColors(selected_probe,:))

    T_probes.tip_active_um(selected_probe) = probe_length_histo *10; %TODO
    T_probes.top_active_um(selected_probe) = active_site_start * 10; %TODO

    fig_probes(selected_probe) = gcf;


    %%%% coordinates, cf. Analyze_Clicked_Points.m
    %
    % m is the brain entry point

    Tapdvml_m = apdvml2info(m, av, st, plane);

    Tapdvml_tip = apdvml2info(curr_probePoints(tip_index,:), av, st, plane);

    % measure the distance 

    tip2surface_mm = sqrt((Tapdvml_tip.ap_mm - Tapdvml_m.ap_mm)^2 + ...
        (Tapdvml_tip.dv_mm - Tapdvml_m.dv_mm)^2 + ...
        (Tapdvml_tip.ml_mm - Tapdvml_m.ml_mm)^2);

    tip2surface_mm_paxinos = sqrt((Tapdvml_tip.ap_mm - Tapdvml_m.ap_mm)^2 + ...
        (Tapdvml_tip.dv_mm_paxinos - Tapdvml_m.dv_mm_paxinos)^2 + ...
        (Tapdvml_tip.ml_mm - Tapdvml_m.ml_mm)^2);
    
    top_active = (m * active_probe_length + ...
        curr_probePoints(tip_index,:) * (tip2surface_mm - active_probe_length))...
        /tip2surface_mm;

    % obtained the information for all the 384 channels


    a = [linspace(curr_probePoints(tip_index,1), top_active(1), 192)', ...
        linspace(curr_probePoints(tip_index,2), top_active(2), 192)', ...
        linspace(curr_probePoints(tip_index,3), top_active(3), 192)'];
    probe_contact_points = zeros(384, 3);
    for i = 1:192
        probe_contact_points(2*i-1,:) = a(i,:);
        probe_contact_points(2*i,:) = a(i,:);
    end
    clear a

    c_t_contacts = cell(384,1);

    theta = acos(p(2)); % Assuming p(2) corresponds to the DV direction

    for i = 1:384
        c_t_contacts{i} = apdvml2info(probe_contact_points(i,:), av, st, plane);
        c_t_contacts{i}.contact_id = repmat(i,height(c_t_contacts{i}));
        c_t_contacts{i}.probe_id = repmat(selected_probe,height(c_t_contacts{i}));
        c_t_contacts{i}.depth_mm = tip2surface_mm - 0.020 * floor((i-1)/2);
        
        % Project the depth along the line
        projected_depth_mm = c_t_contacts{i}.depth_mm * cos(theta);

        % Convert using the transformation for the chosen plane
        projected_depth_mm_paxinos = accf2pxs_mm(projected_depth_mm, plane, 'distance');

        c_t_contacts{i}.depth_mm_paxinos = projected_depth_mm_paxinos;
    end

    c_Tapdvml_contacts{selected_probe} = vertcat(c_t_contacts{:});
    clear c_t_contacts

    pause(.05)
end


Tapdvml_contacts = vertcat(c_Tapdvml_contacts{:});
Tapdvml_contacts.probe_AB = strings(height(Tapdvml_contacts), 1);
Tapdvml_contacts.session_id = strings(height(Tapdvml_contacts), 1);
Tapdvml_contacts.subject_id = strings(height(Tapdvml_contacts), 1);


T_borders = [];
for i = probes
    borders_table{i}.probe_id = repmat(i,height(borders_table{i}),1);
    T_borders = [T_borders; join(borders_table{i}, T_probes(i,:))];
end

i = T_borders.Properties.VariableNames == "upperBorder";
T_borders.Properties.VariableNames{i} = 'upperBorder_um';

i = T_borders.Properties.VariableNames == "lowerBorder";
T_borders.Properties.VariableNames{i} = 'lowerBorder_um';

