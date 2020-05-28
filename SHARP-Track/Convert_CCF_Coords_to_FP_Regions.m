% ----------------------------------------------------------------------------
%
%         Convert Allen CCF indices to Franklin-Paxinos labels
%
%   Based on data from Chon et al. Enhanced and unified anatomical labeling 
%   for a common mouse brain atlas (2020).
%
% ----------------------------------------------------------------------------


%% Start with coordinates within the Allen CCF mouse brain atlas
% in the form [AP1, DV1, ML1
%              AP2, DV2, ML2] 

% This example is of points along a neuropixels probe track through cortex, SC, PAG
brain_points =    [893   189   475
                   890    57   454
                   891   114   468
                   922   156   482
                   923   114   472
                   920   215   492
                   920   268   505
                   954   199   479
                   948   256   491
                   943   312   508
                   938   369   525
                   932   421   537
                   928   471   551
                   951   281   507
                   944   342   522
                   937   400   538
                   931   450   551];
               
% directory of reference files
annotation_volume_location = 'C:\Drive\Histology\annotation_volume_10um_by_index.npy'; % from the allen inst (see readme)
structure_tree_location = 'C:\Drive\Histology\structure_tree_safe_2017.csv'; % located in github repo
CCF_to_FP_location =  'C:\Drive\Histology\CCF_to_FP.csv'; % located in github repo
FP_table_location = 'C:\Drive\Histology\FP_table_Chon_2020.csv'; % located in github repo
chon_images_loc = 'C:\Drive\Histology\Suppl_File1_Labels'; % from chon et al (supplementary data 4, https://www.nature.com/articles/s41467-019-13057-w)

% generate values for pixel-to-coordinate transformation
bregma = allenCCFbregma(); % estimated bregma position in reference data space
atlas_resolution = 0.010; % pixels to mm

% should the brain image be dark or light
black_brain = true;
brain_points_color = [.5 .5 1];


%% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
if ~exist('CCFtoFPtable','var') || ~exist('FPtable','var')
    disp('loading CCF-FP lookup tables...')
    CCFtoFPtable = loadCCFtoFP(CCF_to_FP_location);
    FPtable = loadFPtable(FP_table_location);
end

% initialize array of region annotations
annotation_CCF = cell(size(brain_points,1),3);    
annotation_FP = cell(size(brain_points,1),3);  

%% process data

% loop through every point to get ROI locations and region annotations
for point = 1:size(brain_points,1)

    % find the annotation, name, and acronym of the current point from
    % Allen CCF data
    ann = av(brain_points(point,1),brain_points(point,2),brain_points(point,3));
    name = st.safe_name{ann};
    acr = st.acronym{ann};

    annotation_CCF{point,1} = ann;
    annotation_CCF{point,2} = name;
    annotation_CCF{point,3} = acr;

    % find the annotation, name, and acronym of the current ROI pixel
    % using Chon et al data synthesizing CCF and Franklin-Paxinos
    [ann_FP, name_FP, acr_FP] = CCF_to_FP(brain_points(point,1), brain_points(point,2), brain_points(point,3), ...
                                          CCFtoFPtable, FPtable, chon_images_loc);

    annotation_FP{point,1} = ann_FP;
    annotation_FP{point,2} = name_FP;
    annotation_FP{point,3} = acr_FP;
    
end

% get coordinates relative to bregm\
ap = -(brain_points(:,1)-bregma(1))*atlas_resolution;
dv = (brain_points(:,2)-bregma(2))*atlas_resolution;
ml = (brain_points(:,3)-bregma(3))*atlas_resolution;

% generate table
data_table = table(annotation_CCF(:,2),annotation_CCF(:,3), annotation_FP(:,2),annotation_FP(:,3),...
                        ap,dv,ml, annotation_CCF(:,1), annotation_FP(:,1),...
         'VariableNames', {'CCF_name', 'CCF_abbrv', 'FP_name', 'FP_abbrv', 'AP_location', 'DV_location', 'ML_location', 'CCF_index', 'FP_index'});



%% display results

% plot points on the wire frame brain
fwireframe = plotBrainGrid([], [], [], black_brain); hold on; 
fwireframe.InvertHardcopy = 'off';
figure(fwireframe); hold on
hp = plot3(brain_points(:,1), brain_points(:,3), brain_points(:,2), '.','linewidth',2, 'color',brain_points_color,'markers',10);   

% display table
disp(data_table)




