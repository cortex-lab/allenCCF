% ------------------------------------------------------------------------
%          Run Allen Atlas Browser
% ------------------------------------------------------------------------


%% ENTER FILE LOCATIONS AND PARAMETERS

subject = 'SS099';
date = '2018-03-12';
tag = 'ZO';
probe_length = 5.6; % in mm
active_probe_length = 3.84; % in mm
probe_radius = 100; % in um

% directory of histology
processed_images_folder = 'C:\\Users\\Experiment\\Desktop\\brain volumes\\slices\\SS096\\processed\\';


mark_points = true;
save_points = false;

%% GET PROBE TRAJECTORY POINTS

if ~exist('av')
    disp('loading reference...')
    av = readNPY('\\ZSERVER.cortexlab.net\Lab\Atlas\allenCCF\annotation_volume_10um_by_index.npy');
    st = loadStructureTree('structure_tree_safe_2017.csv');
end

% -----------------------
% either mark new points
% -----------------------
if mark_points
    if ~exist('tv')
        tv = readNPY('\\ZSERVER.cortexlab.net\Lab\Atlas\allenCCF\template_volume_10um.npy');
    end
        
    % show histology
    try; figure(slice_figure);
    catch; slice_figure = figure('Name','Slice Viewer'); end
    sliceBrowser(slice_figure, processed_images_folder);

    % use application
    f = allenAtlasBrowser(tv,av,st, slice_figure);

    %***for current version: put break point here  --> run these lines once probe points are selected***
    % get probe points (before closing the atlas)!
    tracedPoints = f.UserData.pointList(:, [3 2 1]);
    
    if tracedPoints(1,3)>570
        tracedPoints(:,3) = 1140-tracedPoints(:,3);
    end
    
    % save traced points to folder
    if save_points
        save(fullfile(processed_images_folder, sprintf('%s_%s_%s_tracedPoints.mat', subject, date, tag)), 'tracedPoints')
    end
% ------------------
% or load old points
% ------------------    
else
    load(fullfile(processed_images_folder, sprintf('%s_%s_%s_tracedPoints.mat', subject, date, tag)))
end

%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE
    
% get line of best fit through points
% m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
[m,p,s] = best_fit_line(tracedPoints(:,1), tracedPoints(:,2), tracedPoints(:,3));


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
fwireframe = plotBrainGrid([], []); hold on; 
% plot probe points
hp = plot3(tracedPoints(:,1), tracedPoints(:,3), tracedPoints(:,2), 'o','linewidth',2, 'color','b');
% plot brain entry point
plot3(m(1), m(3), m(2), 'k*','linewidth',2)

% find length of probe in reference atlas space
[depth, tip_index] = max(tracedPoints(:,2));
reference_probe_length = sqrt(sum((tracedPoints(tip_index,:) - m).^2)); 
shrinkage_factor = (reference_probe_length / 100) / probe_length;

disp(' ');
disp(['probe length of ' num2str(reference_probe_length/100) ' mm in reference atlas space compared to a reported ' num2str(probe_length) ' mm']);
disp(['shrinkage factor of ' num2str(shrinkage_factor)]);

% plot line the length of the probe in reference space
rpl = round(reference_probe_length);
plot3(m(1)+p(1)*[1 rpl], m(3)+p(3)*[1 rpl], m(2)+p(2)*[1 rpl], ...
    'Color', 'r', 'LineWidth', 2, 'LineStyle',':');

% plot line the length of the active probe sites in reference space
percent_of_tract_with_active_sites = active_probe_length / probe_length;
active_site_start = reference_probe_length*(1-percent_of_tract_with_active_sites);
apl = round([active_site_start  reference_probe_length]);
plot3(m(1)+p(1)*[apl(1) apl(2)], m(3)+p(3)*[apl(1) apl(2)], m(2)+p(2)*[apl(1) apl(2)], ...
    'Color', 'r', 'LineWidth', 3);


%% GET AND PLOT BRAIN REGION LABELS ALONG EXTENT OF PROBE
yc = 10*[0:(rpl-1)];   

error_length = round(probe_radius / 10); %microns error as first number
[borders, fD] = plotLabelsAsProbe(m, p, rpl, yc-20, av, st, 1, error_length, active_site_start*10);
title('Brain Regions Along Probe')

% plot line(s) indicating active site length
active_marker = plot([0 100], [active_site_start*10 active_site_start*10], 'color', 'white', 'LineStyle','--');

