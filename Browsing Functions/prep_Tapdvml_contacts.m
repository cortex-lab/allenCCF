function Tapdvml_contacts = prep_Tapdvml_contacts(T_probes, Tapdvml_contacts)
% short description comes here
%
% SYNTAX
% Tapdvml_contacts = prep_Tapdvml_contacts(T_probes)
%
% longer description may come here
%
% INPUT ARGUMENTS
% T_probes     table
%             Rows for probes (+ optic fibers)
%             With the following variables:
%                probe_id          
%                probe_AB          
%                session_id        
%                subject_id        
%                probe_note        
%                upward_from_tip_um
%                upward_from_tim_um
%
% Tapdvml_contacts 
%             table
%
%
% OUTPUT ARGUMENTS
% Tapdvml_contacts
%             table
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 17-Aug-2023 09:51:42
%
% See also
% doc

arguments
    T_probes table
    Tapdvml_contacts table
end

for i = 1:height(T_probes)

    tf = Tapdvml_contacts.probe_id == i;
    Tapdvml_contacts.session_id(tf) = string(T_probes.session_id(i));
    Tapdvml_contacts.subject_id(tf) = string(T_probes.subject_id(i));
    Tapdvml_contacts.probe_AB(tf) = string(T_probes.probe_AB(i));
    Tapdvml_contacts.probe_note(tf) = string(T_probes.probe_note(i));
    
end


% Bregma's coordinates in grid units

bregma = allenCCFbregma();
bregma_ml = bregma(1);
bregma_dv = bregma(2);
bregma_ap = bregma(3);

%load boundary_mask
load('\\ettina\Magill_lab\Kouichi Nakamura\Analysis\allenCCFdata\boundary_mask_20um_radius.mat', 'boundary_mask') % boundary_mask

% Add a new column 'name_with_margin' to Tapdvml_contacts
Tapdvml_contacts.name_with_margin = Tapdvml_contacts.name;

% Iterate through each row in Tapdvml_contacts
for i = 1:height(Tapdvml_contacts)
    % Convert the coordinates to grid units and adjust based on bregma
    ml = round(Tapdvml_contacts.ml_mm(i) * 100 + bregma_ml);
    dv = round(Tapdvml_contacts.dv_mm(i) * 100 + bregma_dv);
    ap = round(Tapdvml_contacts.ap_mm(i) * 100 + bregma_ap);
    
    if dv <= 0
        dv = 1;
    end

    if ap <= 0
        ap = 1;
    end
    % Check the boundary_mask array
    if boundary_mask(ml, dv, ap)
        Tapdvml_contacts.name_with_margin(i) = {'boundary'};
    end
end
