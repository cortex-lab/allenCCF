% ------------------------------------------------------------------------
%      Analyze non-linear ROIs that were clicked as 'object points'
% ------------------------------------------------------------------------


%% ENTER PARAMETERS AND FILE LOCATION

% file location of object points
save_folder = 'C:\Drive\Histology\for tutorial\Richards\processed';

% directory of reference atlas files
annotation_volume_location = 'C:\Drive\Histology\annotation_volume_10um_by_index.npy'; % from the allen inst (see readme)
structure_tree_location = 'C:\Drive\Histology\structure_tree_safe_2017.csv'; % located in github repo
CCF_to_FP_location =  'C:\Drive\Histology\CCF_to_FP.csv'; % located in github repo
FP_table_location = 'C:\Drive\Histology\FP_table_Chon_2020.csv'; % located in github repo
chon_images_loc = 'C:\Drive\Histology\Suppl_File1_Labels'; % from chon et al (supplementary data 4, https://www.nature.com/articles/s41467-019-13057-w)

% name of the saved object points
object_save_name_suffix = '';

% either set to 'all' or a list of indices from the clicked objects in this file, e.g. [2,3]
objects_to_analyze = 'all';

% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% brain figure black or white
black_brain = true;


%% LOAD THE REFERENCE ANNOTATIONS AND PROBE POINTS

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
if ~exist('CCFtoFPtable','var') || ~exist('FPtable','var')
    CCFtoFPtable = loadCCFtoFP(CCF_to_FP_location);
    FPtable = loadFPtable(FP_table_location);
end

% load object points
objectPoints = load(fullfile(save_folder, ['probe_points' object_save_name_suffix]));

% determine which objects to analyze
if strcmp(objects_to_analyze,'all')
    objects = 1:size(objectPoints.pointList.pointList,1);
else
    objects = objects_to_analyze;
end 


%% BRING UP THE RELEVANT DATA FOR EACH PROBE POINTS, FOR FURTHER ANALYSIS

% initialize cell array containing info on each clicked point
if length(objects) > 1
    roi_annotation = cell(length(objects),1);
    roi_location = cell(length(objects),1);
end

% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm

% plot brain grid
ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 
% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';



for object_num = objects
    
    selected_object = objects(object_num);
        
    % get the object points for the currently analyzed object    
    if strcmp(plane,'coronal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [3 2 1]);
    elseif strcmp(plane,'sagittal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 2 3]);
    elseif strcmp(plane,'transverse')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 3 2]);
    end

    % plot points on the wire frame brain
    figure(fwireframe); hold on
    hp = plot3(curr_objectPoints(:,1), curr_objectPoints(:,3), curr_objectPoints(:,2), '.','linewidth',2, 'color',[ProbeColors(object_num,:) .2],'markers',10);   

    % use the point's position in the atlas to get the AP, DV, and ML coordinates
    ap = -(curr_objectPoints(:,1)-bregma(1))*atlas_resolution;
    dv = (curr_objectPoints(:,2)-bregma(2))*atlas_resolution;
    ml = (curr_objectPoints(:,3)-bregma(3))*atlas_resolution;

    roi_location_curr = [ap dv ml];
    
    % initialize array of region annotations
    roi_annotation_curr = cell(size(curr_objectPoints,1),3);    
    roi_annotation_curr_FP = cell(size(curr_objectPoints,1),3);    
    
    % loop through every point to get ROI locations and region annotations
    for point = 1:size(curr_objectPoints,1)

        % find the annotation, name, and acronym of the current ROI pixel
        ann = av(curr_objectPoints(point,1),curr_objectPoints(point,2),curr_objectPoints(point,3));
        name = st.safe_name{ann};
        acr = st.acronym{ann};

        roi_annotation_curr{point,1} = ann;
        roi_annotation_curr{point,2} = name;
        roi_annotation_curr{point,3} = acr;
        
        % find the annotation, name, and acronym of the current ROI pixel
        [ann_FP, name_FP, acr_FP] = CCF_to_FP(curr_objectPoints(point,1), curr_objectPoints(point,2), curr_objectPoints(point,3), CCFtoFPtable, FPtable, chon_images_loc);

        roi_annotation_curr_FP{point,1} = ann_FP;
        roi_annotation_curr_FP{point,2} = name_FP;
        roi_annotation_curr_FP{point,3} = acr_FP;

    end
    
    % save results in cell array
    if length(objects) > 1
        roi_annotation{object_num} = roi_annotation_curr;
        roi_location{object_num} = roi_location_curr;
    else
        roi_annotation = roi_annotation_curr;
        roi_location = roi_location_curr;
    end
 
    % display results in a table
    disp(['Clicked points for object ' num2str(selected_object)])
    roi_table = table(roi_annotation_curr(:,2),roi_annotation_curr(:,3), roi_annotation_curr_FP(:,2),roi_annotation_curr_FP(:,3),...
                        roi_location_curr(:,1),roi_location_curr(:,2),roi_location_curr(:,3), roi_annotation_curr(:,1), roi_annotation_curr_FP(:,1),...
         'VariableNames', {'CCF_name', 'CCF_acronym', 'FP_name', 'FP_acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex', 'fpIndex'});
     disp(roi_table)
    
end


% now, use roi_location and roi_annotation for your further analyses


