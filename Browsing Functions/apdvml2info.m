function Tapdvml = apdvml2info(apdvml_points, av, st, plane)
% apdvml_points has three columns for ap, dv, and ml
% bregma is 1 x 3 in size
% atlas_resolution in mm
% av  1320 ,800, 1140
% st
%
% SYNTAX
% Tapdvml = apdvml2info(apdvml_points, av, st, plane)
%
% INPUT ARGUMENTS
% apdvml_points
%             n by 3 array
%             [ap, dv, ml]
%             Allen CCF coordinates in 10 micrometers.
%
% av          uint16 (1320 x 800 x 1140)
%             3D array for structure delineation
%
% st          table
%             For structures
%
% plane       'coronal' | 'sagittal' | 'transverse'
%
% OUTPUT ARGUMENTS
% Tapdvml     table
%             Including output ap_mm, dv_mm, and ml_mm
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 06-Nov-2023 14:41:14
%
% See also
% apdvml_mm2ccf

% see also
% Analyze_Clicked_Points.m, accf2pxs_mm

arguments
    apdvml_points (:,3)
    av
    st
    plane (1,1) string {mustBeMember(plane, {'coronal','sagittal','transverse'})}
end


% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm


ap = -(apdvml_points(:,1)-bregma(1))*atlas_resolution;
dv = (apdvml_points(:,2)-bregma(2))*atlas_resolution;
ml = (apdvml_points(:,3)-bregma(3))*atlas_resolution;

% roi_location_curr = [ap dv ml];

% initialize array of region annotations
roi_annotation_curr = cell(size(apdvml_points,1),3);

% loop through every point to get ROI locations and region annotations
tf_rows = true(size(apdvml_points,1),1);
for point = 1:size(apdvml_points,1)

    % find the annotation, name, and acronym of the current ROI pixel
    depth = ceil(apdvml_points(point,2));
    if depth >= 1

        ann = av(ceil(apdvml_points(point,1)), ...
            depth, ... %NOTE this can take the value of 0 and cause an error (must be >= 1)
            ceil(apdvml_points(point,3)));
        name = st.safe_name{ann};
        acr = st.acronym{ann};

        roi_annotation_curr{point,1} = ann;
        roi_annotation_curr{point,2} = name;
        roi_annotation_curr{point,3} = acr;

    else
        % ignore this point, too dorsal
        roi_annotation_curr{point,1} = NaN;
        roi_annotation_curr{point,2} = '';
        roi_annotation_curr{point,3} = '';

    end

end

Tapdvml = [table(ap, dv, accf2pxs_mm(dv, plane), ml, 'VariableNames',{'ap_mm','dv_mm', 'dv_mm_paxinos', 'ml_mm'}), ...
    cell2table(roi_annotation_curr, 'VariableNames', {'annotation', 'name', 'acronym'})];

end