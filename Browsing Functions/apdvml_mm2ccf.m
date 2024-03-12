function Tccf = apdvml_mm2ccf(apdvml_mm)%, av, st, plane)
% Opposite of apdvml2info.
% Convert coordinates in mm back to Allen CCF coodinates in 10 micrometers
%
% SYNTAX
% v = apdvml_mm2ccf(apdvml_mm)
%
%
% INPUT ARGUMENTS
% apdvml_mm   n by 3 array
%             [ap_mm, dv_mm, ml_mm]
%             a coordinate in mm
%
% % av          uint16 (1320 x 800 x 1140)
% %             3D array for structure delineation
% %
% % st          table
% %             For structures
% %
% % plane       'coronal' | 'sagittal' | 'transverse'
%
% OUTPUT ARGUMENTS
% Tccf        table
%             including Allen CCF coordinates in 10 micrometers.
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 06-Nov-2023 14:35:03
%
% See also
% apdvml2info


arguments
    apdvml_mm (:,3)
    % av
    % st
    % plane (1,1) string {mustBeMember(plane, {'coronal','sagittal','transverse'})}
end

% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm


%TODO

ap = (- apdvml_mm(:,1))/atlas_resolution + bregma(1);
dv = (apdvml_mm(:,2))/atlas_resolution + bregma(2);
ml = (apdvml_mm(:,3))/atlas_resolution + bregma(3);

% roi_location_curr = [ap dv ml];

% % initialize array of region annotations
% roi_annotation_curr = cell(size(apdvml_points,1),3);
% 
% % loop through every point to get ROI locations and region annotations
% tf_rows = true(size(apdvml_mm,1),1);
% for point = 1:size(apdvml_mm,1)
% 
%     % find the annotation, name, and acronym of the current ROI pixel
%     depth = ceil(apdvml_points(point,2));
%     if depth >= 1
% 
%         ann = av(ceil(apdvml_points(point,1)), ...
%             depth, ... %NOTE this can take the value of 0 and cause an error (must be >= 1)
%             ceil(apdvml_points(point,3)));
%         name = st.safe_name{ann};
%         acr = st.acronym{ann};
% 
%         roi_annotation_curr{point,1} = ann;
%         roi_annotation_curr{point,2} = name;
%         roi_annotation_curr{point,3} = acr;
% 
%     else
%         % ignore this point, too dorsal
%         roi_annotation_curr{point,1} = NaN;
%         roi_annotation_curr{point,2} = '';
%         roi_annotation_curr{point,3} = '';
% 
%     end
% 
% end

Tccf = table(ap, dv, ml, 'VariableNames',{'ap_mm','dv_mm', 'ml_mm'});%, ...
    %cell2table(roi_annotation_curr, 'VariableNames', {'annotation', 'name', 'acronym'})];

end