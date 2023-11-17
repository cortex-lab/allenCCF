function updown_probe_with_slider(av, st, session_id, probeAB, plane, sessions_dir,...
    task_dir_name, imaging_session_dir, image_folder_name, depth_level, ...
    active_probe_length, probe_insertion_direction)
% Move up and down a probe using a slider to match the anatomical borders
% to the changes in the firing rates of neurons
%
% SYNTAX
% updown_probe_with_slider(av, st, session_id, probeAB, plane, sessions_dir,...
%    task_dir_name, imaging_session_dir, image_folder_name, depth_level, ...
%    active_probe_length, probe_insertion_direction)
%
%
% INPUT ARGUMENTS
% av          uint16
%             part of Allen CCF data    
%
% st          table
%             part of Allen CCF data
%
% session_id  string
%             eg. "kms058-2023-03-25-184034"
%
% probeAB     "ProbeA" | "ProbeB"
%
% plane       "coronal" | "sagittal" | "transverse"
%
% sessions_dir
%             folder path
%             "\\ettina\Magill_Lab\Julien\Data\head-fixed\by_sessions"
%
% task_dir_name
%             folder name
%             "reaching_go_spout_bar_nov22"
%
% imaging_session_dir
%             folder_path
%             "\\ettina\Magill_Lab\Kouichi Nakamura\Analysis\Images from Otto\20230406 kms058"
%
% image_folder_name
%             "RGB" (default)
%             folder_name
%
% depth_level
%             6 (default)
%             The hierarchical level of the anatomical structures to be
%             shown.
%
% active_probe_length
%             3.840 (default)    
%
% probe_insertion_direction 
%             "down" (default) | "up"
%
%
%
%
% EXAMPLE USAGE
%
% session_id = "kms058-2023-03-25-184034"
% probeAB ="ProbeA";
% plane = "sagittal";
% sessions_dir = "\\ettina\Magill_Lab\Julien\Data\head-fixed\by_sessions";
% task_dir_name = "reaching_go_spout_bar_nov22"
% imaging_session_dir = "\\ettina\Magill_Lab\Kouichi Nakamura\Analysis\Images from Otto\20230406 kms058";
% image_folder_name = 'RGB_ignore';
% 
% % load the reference brain and region annotations
% if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
%     disp('loading reference atlas...')
%     av = readNPY(annotation_volume_location);
%     st = loadStructureTree(structure_tree_location);
%     tv = readNPY(template_volume_location);
% end
%
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 17-Nov-2023 09:07:25
%
% See also
% plot_and_compute_probe_positions_from_ccf

arguments
    av uint16
    st table
    session_id (1,1) string % eg. "kms058-2023-03-25-184034"
    probeAB (1,1) string {mustBeMember(probeAB,["ProbeA","ProbeB"])} % eg. "ProbeA"
    plane (1,1) string {mustBeMember(plane,["coronal","sagittal", "transverse"])} % eg. "sagittal"
    sessions_dir (1,1) string  {mustBeFolder} % eg. "\\ettina\Magill_Lab\Julien\Data\head-fixed\by_sessions"
    task_dir_name (1,1) string % eg. "reaching_go_spout_bar_nov22"
    imaging_session_dir (1,1) string  {mustBeFolder} % eg. "\\ettina\Magill_Lab\Kouichi Nakamura\Analysis\Images from Otto\20230406 kms058"
    image_folder_name (1,1) string = "RGB"
    depth_level (1,1) double {mustBePositive, mustBeInteger} = 6
    active_probe_length (1,1) double = 3.84 % in mm
    probe_insertion_direction (1,1) string {mustBeMember(probe_insertion_direction, ["down","up"])} = "down"
end

%TODO add button to go next and previous probe

%TODO what if optic fiber data is chosen?




sorter_output_dir = fullfile(sessions_dir, task_dir_name,  ...
    session_id, "processed\kilosort", probeAB, "sorter_output");

recording_cell_metrics_path = fullfile(sorter_output_dir, "recording.cell_metrics.cellinfo.mat");
s = load(recording_cell_metrics_path);
recording_cell_metrics = s.cell_metrics;

image_folder =fullfile(imaging_session_dir, image_folder_name);

probe_save_name_suffix = '_probe';

Tapdvml_contacts_path = fullfile(imaging_session_dir, "Tapdvml_contacts.xlsx");

original_backup_path = regexprep(Tapdvml_contacts_path, '\.xlsx$' , "_original.xlsx");

% borders_table_path= fullfile(imaging_session_dir, "borders_table.xlsx");

Tprobes_path= fullfile(imaging_session_dir, "T_probes.xlsx");

%YLim_shift_mm_path = fullfile(imaging_session_dir,"YLim_shift_mm.mat"); %TODO
Y_shift_path = fullfile(imaging_session_dir,"T_Y_shift.mat"); %TODO





% Tborders_table = readtable(borders_table_path);

T_probes=  readtable(Tprobes_path);

if ~ ismember('p_1',T_probes.Properties.VariableNames) || ...
        ~ ismember('p_2',T_probes.Properties.VariableNames) ||...
        ~ ismember('p_3',T_probes.Properties.VariableNames)
    warning('The eigen vectors p is not found in T_probes. You cannot update the Tapdvml_contacts.')
end

switch probeAB
    case "ProbeA"
        probeAB_short = "A";
    case "ProbeB"
        probeAB_short = "B";
    otherwise

end

tf = T_probes.session_id == session_id & T_probes.probe_AB == probeAB_short;

probe_id = T_probes{tf, "probe_id"};

%% deal with the cache file
opts = detectImportOptions(Tapdvml_contacts_path);
opts = setvartype(opts, 'session_id', 'string');

YLim_shift_mm = [];
if isfile(original_backup_path)

    Tapdvml_contacts = readtable(original_backup_path, opts );
    Tapdvml_contacts_this = Tapdvml_contacts(Tapdvml_contacts.probe_id == probe_id, :);

    if isfile(Y_shift_path) % use the cache
        S = load(Y_shift_path); % YLim_shift_mm

        assert(istable(S.T_Y_shift))
        T_Y_shift = S.T_Y_shift;
        YLim_shift_mm = S.T_Y_shift.YLim_shift_mm(probe_id);
        Tanc = S.T_Y_shift.Tanc{probe_id};
        if isempty(Tanc)
            Tanc = compute_Tanc(Tapdvml_contacts_this, st, depth_level);
            T_Y_shift.Tanc{probe_id} = Tanc;
            save(Y_shift_path, 'T_Y_shift');
        end

        clear S
    else
        % cache file is not found
        Tanc_ = cell(height(T_probes), 1);
        for i = 1:size(Tanc_,1)
            Tanc_{i} = table();
        end

        T_Y_shift = table(zeros(height(T_probes), 1), Tanc_,...
            'VariableNames', ["YLim_shift_mm","Tanc"]); % YLim_shift_mm and Tanc
        clear Tanc_
        Tanc = compute_Tanc(Tapdvml_contacts_this, st, depth_level);
        T_Y_shift.Tanc{probe_id} = Tanc;
        save(Y_shift_path, 'T_Y_shift');
    end
    Tapdvml_contacts_this = [Tapdvml_contacts_this, Tanc];

else
    % the very first run

    Tapdvml_contacts = readtable(Tapdvml_contacts_path, opts );
    Tapdvml_contacts_this = Tapdvml_contacts(Tapdvml_contacts.probe_id == probe_id, :);

    Tanc_ = cell(height(T_probes), 1);
    for i = 1:size(Tanc_,1)
        Tanc_{i} = table();
    end

    T_Y_shift = table(zeros(height(T_probes), 1), Tanc_,...
        'VariableNames', ["YLim_shift_mm","Tanc"]); % YLim_shift_mm and Tanc
    clear Tanc_
    Tanc = compute_Tanc(Tapdvml_contacts_this, st, depth_level);
    T_Y_shift.Tanc{probe_id} = Tanc;
    save(Y_shift_path, 'T_Y_shift');

    Tapdvml_contacts_this = [Tapdvml_contacts_this, Tanc];

end


%% Figure

ufig = uifigure(Position=[800, 150, 720, 1080]);

uax1 = uiaxes(ufig);
uax1.PositionConstraint = 'innerposition';
uax1.Position=[100, 150, 350, 800];
ylabel(uax1, 'Dsitance from the tip (mm)')
uax1.XTick =[];

if ~isempty(YLim_shift_mm)
    assert(isequal(size(YLim_shift_mm), [1 1]))
    uax1.UserData.YLim_shift_mm = YLim_shift_mm;
end

uax2 = uiaxes(ufig);
uax2.PositionConstraint = 'innerposition';
uax2.Position=[448, 150, 80, 800];
uax2.YTickLabel =[];
uax2.XTick =[];


uax3 = uiaxes(ufig);
uax3.PositionConstraint = 'innerposition';
uax3.Position=[525, 150, 120, 800];
uax3.YTickLabel =[];
uax3.XTick =[];


tickdir([uax1,uax2, uax3], 'out')
ticklengthcm([uax1,uax2, uax3], 0.2)

% usl1 = uislider(ufig, Orientation='vertical' , Limits=[-50, 50]);
% usl1.Position = [50, 200, 3, 300];

bg1 = uibuttongroup(ufig, Title='Move');
bg1.Position = [20, 150, 70, 200];

bt1 = uibutton(bg1, "Text","+5", "ButtonPushedFcn", @(src, event) move_vertically(uax1, 5));
bt1.Position = [10, 140, 50, 30];

bt2 = uibutton(bg1, "Text","+1", "ButtonPushedFcn", @(src, event) move_vertically(uax1, 1));
bt2.Position = [10, 110, 50, 30];

bt3 = uibutton(bg1, "Text","-1", "ButtonPushedFcn", @(src, event) move_vertically(uax1, -1));
bt3.Position = [10, 40, 50, 30];

bt4 = uibutton(bg1, "Text","-5", "ButtonPushedFcn", @(src, event) move_vertically(uax1, -5));
bt4.Position = [10, 10, 50, 30];

bt6 = uibutton(bg1, "Text","0", "ButtonPushedFcn", @(src, event) move_back_to_original(uax1));
bt6.Position = [10, 75, 50, 30];


bt7 = uibutton(ufig, "Text","Save", "ButtonPushedFcn", @(src, event) save_updown_to_table(uax1));
bt7.Position = [150, 30, 300, 50];


%% Plot cells
n = length(recording_cell_metrics.firingRate);

x_val = rand(n,1);

y = recording_cell_metrics.maxWaveformCh * 0.010; % mm

c = recording_cell_metrics.firingRate;

sc1 = scatter(uax3, x_val, y, 36, c, LineWidth=1.5);
uax3.Color = 'k';

ylim(uax3, [-0.100, 3.940])

cb1 =colorbar(uax3,'eastoutside');
ylabel(cb1,'Mean firing rates (spikes/s)' )
cb1.Label.FontSize = 14;
cb1.TickDirection = 'out';
ticklengthcm(cb1, 0.2)


%% Moving average of five channels
% 
% SDF like plot with Gaussian kernel convolution is better
cla(uax2,'reset')

movingaverageval = zeros(384,1);
celldensity = zeros(384,1);

for i = 1:384

    % moving average of five channels
    chans = i - 2 : i + 2;
    
    ind = ismember(recording_cell_metrics.maxWaveformCh, chans);
    celldensity(i) = nnz(ind)/(length(chans)*10); % cells per micrometers
    movingaverageval(i) = mean(recording_cell_metrics.firingRate(ind));

end

for i = 1:384

    % ten channels
    chans = i - 4 : i + 5;
    
    ind = ismember(recording_cell_metrics.maxWaveformCh, chans);
    celldensity(i) = nnz(ind)/(length(chans)*10); % cells per micrometers

end

isc1 = imagesc(uax2, [0 1], (1:384)'*0.01, movingaverageval);

uax2.YDir = 'normal';
ylim(uax2, [-0.100, 3.940])
uax2.XTick =[];

uax2.Color = 'k';

uax3.CLim = uax2.CLim; % same color scale

ln1 = line(uax2, celldensity/max(celldensity), (1:384)'*0.01, 'Color', 'w');
xlim(uax2,[0 1])

tx1 = annotation(ufig, 'textbox', [0.2 0.85 0.8 0.1],  'String', ...
    sprintf("%s/%s (%d)", session_id, probeAB, probe_id),...
    HorizontalAlignment='center',FontSize=18,LineStyle='none');

%% 
anc_id = Tapdvml_contacts_this.anc_id;

anc_name = Tapdvml_contacts_this.anc_name;

anc_color = Tapdvml_contacts_this.anc_color;


anc_id_diff = diff(anc_id);

start_ind = [1; find(anc_id_diff) + 1];

end_ind = [start_ind(2:end) - 1; length(anc_id)];

block_names = anc_name(start_ind);

block_colors = anc_color(start_ind, :);


Tblocks = table(block_names, start_ind, end_ind, block_colors, ...
    VariableNames={'block_names', 'start_ind','end_ind', 'block_colors'});


%% Plot structures

cla(uax1)

uax1.Color = 'k';

for i = 1:height(Tblocks)
    X = [0 1 1 0];
    Y = [Tblocks{i,'start_ind'}-0.5, Tblocks{i,'start_ind'}-0.5,...
         Tblocks{i,'end_ind'}+0.5, Tblocks{i,'end_ind'}+0.5]/100; % in mm
    patch(uax1, X, Y, Tblocks.block_colors(i,:));
    text(uax1, 0.5, mean([Tblocks{i,'start_ind'}, Tblocks{i,'end_ind'}])/100, ...
        Tblocks.block_names(i),...
        'FontSize',14,...
        'HorizontalAlignment','center', Clipping='on');
end

horzline = line(uax1, [0 1], [0 0], 'Color','r'); % distance from the tip 0 mm

ylim(uax1, [-0.100, 3.940])

uax1.YLabel.FontSize = 14;

if isempty(uax1.UserData) || ~isfield(uax1.UserData, 'YLim_shift_mm')
    tx2 = annotation(ufig, 'textbox', [0.2 0.03 0.3 0.1],  'String', ...
        regexprep(sprintf("%+.2f mm", 0), '-', char(8722)),...
        HorizontalAlignment='left',FontSize=12, LineStyle='none');
else
    tx2 = annotation(ufig, 'textbox', [0.2 0.03 0.3 0.1],  'String', ...
        regexprep(sprintf("%+.2f mm", uax1.UserData.YLim_shift_mm), '-', char(8722)),...
        HorizontalAlignment='left',FontSize=12, LineStyle='none');
end

% move items after initialiation
if isfield(uax1.UserData, 'YLim_shift_mm')
    % uax1.UserData.YLim_shift_mm

    ch = uax1.Children;

    for i = 1:length(ch)
        if isa(ch(i),'matlab.graphics.primitive.Text')
            ch(i).Position(2) = ch(i).Position(2) + YLim_shift_mm;

        elseif isa(ch(i),'matlab.graphics.primitive.Patch')
            ch(i).YData = ch(i).YData + YLim_shift_mm;

        elseif isa(ch(i),'matlab.graphics.primitive.Line')
            ch(i).YData = ch(i).YData + YLim_shift_mm;            
        
        end
    end
end


    function move_back_to_original(uax1)

        if isfield(uax1.UserData, 'YLim_shift_mm')
            ch = uax1.Children;

            for j = 1:length(ch)
                step = - uax1.UserData.YLim_shift_mm;

                if isa(ch(j),'matlab.graphics.primitive.Text')
                    ch(j).Position(2) = ch(j).Position(2) + step;

                elseif isa(ch(j),'matlab.graphics.primitive.Patch')
                    ch(j).YData = ch(j).YData + step;

                elseif isa(ch(j),'matlab.graphics.primitive.Line')
                    ch(j).YData = ch(j).YData + step;
                    
                end

            end

            uax1.UserData.YLim_shift_mm = 0;

            tx2.String = regexprep(sprintf("%+.2f mm", uax1.UserData.YLim_shift_mm), '-', char(8722));
        else
            %nothing to do
        end

    end


    function move_vertically(uax1, step)

        % arguments
        %     uax1(1,1)
        %     step (1,1) {mustBeInteger}
        % end

        ch = uax1.Children;

        for j = 1:length(ch)

            if isa(ch(j),'matlab.graphics.primitive.Text')
                ch(j).Position(2) = ch(j).Position(2) + step*0.01;

            elseif isa(ch(j),'matlab.graphics.primitive.Patch')
                ch(j).YData = ch(j).YData + step*0.01;

            elseif isa(ch(j),'matlab.graphics.primitive.Line')
                ch(j).YData = ch(j).YData + step*0.01;

            end

        end

        if isfield(uax1.UserData, 'YLim_shift_mm')
            uax1.UserData.YLim_shift_mm = uax1.UserData.YLim_shift_mm + step*0.01;
        else
            uax1.UserData.YLim_shift_mm = step*0.01;
        end

        tx2.String = regexprep(sprintf("%+.2f mm", uax1.UserData.YLim_shift_mm), '-', char(8722));

    end


    function save_updown_to_table(uax1, ~)

        if isempty(uax1.UserData)
            disp('The probe has not moved on GUI yet. No change made.')
            return
        else
            if isfield(uax1.UserData, 'YLim_shift_mm')
                % move
            else
                disp('YLim_shift_mm is not found. No change made.')
                return
            end
        end


        %% save the backup

        assert(endsWith(Tapdvml_contacts_path,'.xlsx'))

        assert(isfile(Tapdvml_contacts_path))

        f_info = dir(Tapdvml_contacts_path);

        dt = datetime(f_info.date, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');

        dt.Format = 'yyyy-MM-dd_HH-mm-ss';
        suffix = string(dt);

        new_file_path = regexprep(Tapdvml_contacts_path, '\.xlsx$' , "_" + suffix + ".xlsx");

        if isfile(original_backup_path)
            % save a backup file for the latest

            copyfile(Tapdvml_contacts_path, new_file_path)

        else
            % this is the first time you change the Excel file
            % save a backup file

            copyfile(Tapdvml_contacts_path, original_backup_path)

        end

        %% Job
        % see plot_and_compute_probe_positions_from_ccf.m
        % I need the eigen vector, to be saved in T_probes?

        switch probeAB
            case "ProbeA"
                probeAB_ = "A";
            case "ProbeB"
                probeAB_ = "B";
            case "optic fiber"
                probeAB_ = "optic fiber";
        end

        % eigen vector
        tf = T_probes.session_id == session_id & ...
            T_probes.probe_AB == probeAB_;

        assert(nnz(tf) == 1)

        p = T_probes{T_probes.session_id == session_id & ...
            T_probes.probe_AB == probeAB_, ["p_1","p_2","p_3"]};

        m = T_probes{T_probes.session_id == session_id & ...
            T_probes.probe_AB == probeAB_, ["m_1","m_2","m_3"]};

        processed_images_folder = fullfile(image_folder, 'processed');

        probePoints = load(fullfile(processed_images_folder, ['probe_points' probe_save_name_suffix]));

        % get the probe points for the currently analyzed probe

        if strcmp(plane,'coronal')
            curr_probePoints = probePoints.pointList.pointList{probe_id,1}(:, [3 2 1]);
        elseif strcmp(plane,'sagittal')
            curr_probePoints = probePoints.pointList.pointList{probe_id,1}(:, [1 2 3]);
        elseif strcmp(plane,'transverse')
            curr_probePoints = probePoints.pointList.pointList{probe_id,1}(:, [1 3 2]);
        end

        % m is the brain entry point

        % modifiy the m and curr_probePoints(tip_index,:) by YLim_shift_mm
        % p is in 10 µm

        if strcmp(probe_insertion_direction, 'down')
            [depth, tip_index] = max(curr_probePoints(:,2));
        elseif strcmp(probe_insertion_direction, 'up')
            [depth, tip_index] = min(curr_probePoints(:,2));
        end

        % curr_probePoints(tip_index,:) either to get it from the ProbePoints or
        % from Tapdvml_contacts.xlsx (convert ap_mm, dv_mm, ml_mm back to allen coordinates)
        % opposite of apdvml2info


        % Tapdvml_contacts.session_id

        % tip_mm =

        YLim_shift_mm = uax1.UserData.YLim_shift_mm;
        new_cp = curr_probePoints(tip_index,:) - YLim_shift_mm*100; % minus to move up
        %NOTE always read the original Excel and saved value of YLim_shift_mm
        % YLim_shift_mm should represent the gap from the original

        Tapdvml_m = apdvml2info(m, av, st, plane);

        Tapdvml_tip = apdvml2info(new_cp, av, st, plane);



        %% measure the distance

        tip2surface_mm = sqrt((Tapdvml_tip.ap_mm - Tapdvml_m.ap_mm)^2 + ...
            (Tapdvml_tip.dv_mm - Tapdvml_m.dv_mm)^2 + ...
            (Tapdvml_tip.ml_mm - Tapdvml_m.ml_mm)^2);

        % tip2surface_mm_paxinos = sqrt((Tapdvml_tip.ap_mm - Tapdvml_m.ap_mm)^2 + ...
        %     (Tapdvml_tip.dv_mm_paxinos - Tapdvml_m.dv_mm_paxinos)^2 + ...
        %     (Tapdvml_tip.ml_mm - Tapdvml_m.ml_mm)^2);

        top_active = (m * active_probe_length + ...
            new_cp * (tip2surface_mm - active_probe_length))...
            /tip2surface_mm;

        %% obtain the information for all the 384 channels

        a = [linspace(new_cp(1), top_active(1), 192)', ...
            linspace(new_cp(2), top_active(2), 192)', ...
            linspace(new_cp(3), top_active(3), 192)'];
        probe_contact_points = zeros(384, 3);
        for j = 1:192
            probe_contact_points(2*j-1,:) = a(j,:);
            probe_contact_points(2*j,:) = a(j,:);
        end
        clear a

        c_t_contacts = cell(384,1);

        theta = acos(p(2)); % Assuming p(2) corresponds to the DV direction

        for j = 1:384
            c_t_contacts{j} = apdvml2info(probe_contact_points(j,:), av, st, plane);
            c_t_contacts{j}.contact_id = repmat(j,height(c_t_contacts{j}));
            c_t_contacts{j}.probe_id = repmat(probe_id,height(c_t_contacts{j}));
            c_t_contacts{j}.depth_mm = tip2surface_mm - 0.020 * floor((j-1)/2);

            % Project the depth along the line
            projected_depth_mm = c_t_contacts{j}.depth_mm * cos(theta);

            % Convert using the transformation for the chosen plane
            projected_depth_mm_paxinos = accf2pxs_mm(projected_depth_mm, plane, 'distance');

            c_t_contacts{j}.depth_mm_paxinos = projected_depth_mm_paxinos;
        end

        Tapdvml_contacts_new = vertcat(c_t_contacts{:});
        clear c_t_contacts

        %% update the relevant rows of Tapdvml_contacts

        %% update values
        for j = 1:height(Tapdvml_contacts_new)
            probe_id = Tapdvml_contacts_new{j,"probe_id"};
            contact_id = Tapdvml_contacts_new{j,"contact_id"};

            cols = ["ap_mm","dv_mm","dv_mm_paxinos","ml_mm","annotation",...
                "name","acronym","contact_id","probe_id","depth_mm","depth_mm_paxinos"];
            Tapdvml_contacts(Tapdvml_contacts.probe_id == probe_id ...
                & Tapdvml_contacts.contact_id == contact_id,  cols) = ...
                Tapdvml_contacts_new(j, cols);
        end

        writetable(Tapdvml_contacts, fullfile(imaging_session_dir,"Tapdvml_contacts.xlsx"),'FileType','spreadsheet')
        fprintf("Tapdvml_contacts.xlsx has been updated with %.3f mm shift from the original for %s of %s (%d).\n", YLim_shift_mm, probeAB, session_id, probe_id)
        % save the value of YLim_shift_mm
        T_Y_shift.YLim_shift_mm(probe_id) = YLim_shift_mm;
        save(Y_shift_path, 'T_Y_shift');


    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [anc, this] = find_name_for_depth(st, depth, name)

arguments
    st table
    depth (1,1) double {mustBeInteger,mustBeNonnegative}
    name (1,1) string
end

st.name = regexprep(st.name, '["'']','');

id = st.id(st.name == name);


cstr = arrayfun(@(x) strsplit(x, "/"), string(st.structure_id_path),'UniformOutput',false);

% remove empty string ""

for j = 1:size(cstr,1) 

    if cstr{j}(1) == ""
        cstr{j}(1) = [];
    end
    if cstr{j}(end) == ""
        cstr{j}(end) = [];
    end    
    cstr{j} = arrayfun(@(x) int32(str2double(x)), cstr{j});%SUPER SLOW

end


st.structure_id_path_vec = cstr;


anc_id_path = st.structure_id_path_vec{st.id == id};
if length(anc_id_path) <= depth
    anc.anc_id = anc_id_path(end);
else
    anc.anc_id = anc_id_path(depth + 1);
end

anc.anc_name = string(st.name(st.id == anc.anc_id));
anc.anc_color_hex = string(st.color_hex_triplet(st.id == anc.anc_id,:));
anc.anc_color_255 = [hex2dec(anc.anc_color_hex{1}(1:2)), hex2dec(anc.anc_color_hex{1}(3:4)), hex2dec(anc.anc_color_hex{1}(5:6))];
anc.anc_color = anc.anc_color_255/255;
anc.index = st.index(st.id == anc.anc_id,:);

this.id = id;
this.name = name;
this.color_hex = string(st.color_hex_triplet(st.id == id,:));
this.color_255 = [hex2dec(this.color_hex{1}(1:2)), hex2dec(this.color_hex{1}(3:4)), hex2dec(this.color_hex{1}(5:6))];
this.color = this.color_255/255;


end


function Tanc = compute_Tanc(Tapdvml_contacts_this, st, depth_level)

anc = preallocatestruct(["index", "anc_id", "anc_name", "anc_color_hex", "anc_color_255", "anc_color"], [height(Tapdvml_contacts_this), 1]);

st.name = regexprep(st.name, '["'']',''); % remove clutters

f = waitbar(0,'1','Name','Analysing structure hierarchy...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

for i = 1:height(Tapdvml_contacts_this) % SLOW
    waitbar( i/height(Tapdvml_contacts_this), f, ...
        sprintf("%d of %d", i, height(Tapdvml_contacts_this)));
    name = Tapdvml_contacts_this.name{i};
    [anc(i)] = find_name_for_depth(st, depth_level, name);
end

delete(f)
beep();

Tanc = struct2table(anc, AsArray=true);

end