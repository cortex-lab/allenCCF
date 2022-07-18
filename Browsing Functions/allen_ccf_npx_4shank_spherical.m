function allen_ccf_npx_4shank_spherical(tv,av,st)
% allen_ccf_npx(tv,av,st)
% Andy Peters (peters.andrew.j@gmail.com)
% bunch of 4-shank modifications: Flora 

% GUI for planning Neuropixels trajectories with the Allen CCF
% Part of repository: https://github.com/cortex-lab/allenCCF
% - directions for installing atlas in repository readme
% - some dependent functions from that repository
%
% (optional inputs - if CCF path written in Line 22, loaded automatically)
% tv, av, st = CCF template volume, annotated volume, structure tree

% Initialize gui_data structure
gui_data = struct;

% Allen CCF-bregma transform (estimated from eyeballing Paxinos->CCF)
% [AP,DV,ML]
bregma = [540,0,570];

% If not already loaded in, load in atlas
if nargin < 3
    allen_atlas_path = 'C:\Users\Andrew\OneDrive for Business\Documents\Atlases\AllenCCF';
    if isempty(allen_atlas_path)
        error('Enter path where Allen CCF is stored at Line 23');
    end
    tv = readNPY([allen_atlas_path filesep 'template_volume_10um.npy']); % grey-scale "background signal intensity"
    av = readNPY([allen_atlas_path filesep 'annotation_volume_10um_by_index.npy']); % the number at each pixel labels the area, see note below
    st = loadStructureTree([allen_atlas_path filesep 'structure_tree_safe_2017.csv']); % a table of what all the labels mean    
end

% Load the colormap (located in the repository, find by associated fcn)
allenCCF_path = fileparts(which('allenCCFbregma'));
cmap_filename = [allenCCF_path filesep 'allen_ccf_colormap_2017.mat'];
load(cmap_filename);

% Set up the gui
probe_atlas_gui = figure('Toolbar','none','Menubar','none','color','w', ...
    'Name','Atlas-probe viewer','Units','normalized','Position',[0.2,0.2,0.7,0.7]);

% Set up the atlas axes
axes_atlas = subplot(1,5,1); 
[~, brain_outline] = plotBrainGrid([],axes_atlas);
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])


% 4 shank parameters 
shank_spacing=0:20:60; 
shank_no=4;

for i=1:shank_no
% Set up the probe area axes                    % add 4 of these      
    axes_probe_areas(i) = subplot(1,5,i+1);
    axes_probe_areas(i).ActivePositionProperty = 'position';
    set(axes_probe_areas(i),'FontSize',11);
    yyaxis(axes_probe_areas(i),'left');
    probe_areas_plot(i) = image(0);
    set(axes_probe_areas(i),'XTick','','YLim',[0,2880],'YColor','k','YDir','reverse');
    ylabel(axes_probe_areas(i),'');
    
    colormap(axes_probe_areas(i),cmap);
    caxis([1,size(cmap,1)])
    yyaxis(axes_probe_areas(i),'right');
    set(axes_probe_areas(i),'XTick','','YLim',[0,2880],'YColor','k','YDir','reverse');
    set(axes_probe_areas(i),'YTicklabels',[])
    title(axes_probe_areas(i),sprintf('shank %d',i));
end
yyaxis(axes_probe_areas(1),'left');
ylabel(axes_probe_areas(1),'Depth (\mum)');

% Position the axes              % change positioning so that there are 4
set(axes_atlas,'Position',[-0.15,-0.1,.8,1.2]);
set(axes_probe_areas(1),'Position',[.55,0.05,0.01,0.8]); 
set(axes_probe_areas(2),'Position',[.67,0.05,0.01,0.8]); 
set(axes_probe_areas(3),'Position',[.79,0.05,0.01,0.8]); 
set(axes_probe_areas(4),'Position',[.91,0.05,0.01,0.8]); 


% Set the current axes to the atlas (dirty, but some gca requirements)
axes(axes_atlas);


% Set up the probe reference/actual % again set up the reference so that 4
% shans are plotted
probe_ref_top = [bregma(1),bregma(3),0];
probe_ref_bottom = [bregma(1),bregma(3),size(tv,2)];
probe_ref_vector = [probe_ref_top',probe_ref_bottom'];
shank_spacing_vector=repmat(shank_spacing,2,1)'; 
allshanks_ref_vector=repmat(probe_ref_vector(2,:),shank_no,1)+shank_spacing_vector;

probe_ref_line = line(probe_ref_vector(1,:),allshanks_ref_vector,probe_ref_vector(3,:), ...
    'linewidth',1.5,'color','r','linestyle','--');

% probe_ref_line = line(probe_ref_vector(1,:),probe_ref_vector(2,:),probe_ref_vector(3,:), ...
%     'linewidth',1.5,'color','r','linestyle','--');

probe_length = 288.8; % IMEC phase 2.0 (in 10 ums)
probe_vector = [probe_ref_vector(:,1),diff(probe_ref_vector,[],2)./ ...
    norm(diff(probe_ref_vector,[],2))*probe_length + probe_ref_vector(:,1)];

allshanks_vector=repmat(probe_vector(2,:),shank_no,1)+shank_spacing_vector;
probe_line = line(probe_vector(1,:),allshanks_vector,probe_vector(3,:),'linewidth',3,'color','b','linestyle','-');
set(probe_line(1),'color','g'); % set shank 0 to be a different color 

%probe_line = line(probe_vector(1,:),probe_vector(2,:),probe_vector(3,:),'linewidth',3,'color','b','linestyle','-');

% Set up the text to display coordinates
probe_coordinates_text = uicontrol('Style','text','String','', ...
    'Units','normalized','Position',[0,0.95,1,0.05], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontSize',12);

% Store data
gui_data.tv = tv; % Intensity atlas
gui_data.av = av; % Annotated atlas
gui_data.st = st; % Labels table
gui_data.cmap = cmap; % Atlas colormap
gui_data.bregma = bregma; % Bregma for external referencing
gui_data.probe_length = probe_length; % Length of probe
gui_data.structure_plot_idx = []; % Plotted structures
gui_data.probe_angle = [0;0]; % Probe angles in azimuth/elevation
gui_data.shank_no=shank_no; 
gui_data.spin=0;
gui_data.shank_spacing=shank_spacing;

%Store handles
gui_data.handles.cortex_outline = brain_outline;
gui_data.handles.structure_patch = []; % Plotted structures
gui_data.handles.axes_atlas = axes_atlas; % Axes with 3D atlas
gui_data.handles.axes_probe_areas = axes_probe_areas; % Axes with probe areas
gui_data.handles.slice_plot = surface('EdgeColor','none'); % Slice on 3D atlas
gui_data.handles.slice_volume = 'tv'; % The volume shown in the slice
gui_data.handles.probe_ref_line = probe_ref_line; % Probe reference line on 3D atlas
gui_data.handles.probe_line = probe_line; % Probe reference line on 3D atlas
gui_data.handles.probe_areas_plot = probe_areas_plot; % Color-coded probe regions
gui_data.probe_coordinates_text = probe_coordinates_text; % Probe coordinates text

% Make 3D rotation the default state (toggle on/off with 'r')
h = rotate3d(axes_atlas);
h.Enable = 'on';
% Update the slice whenever a rotation is completed
h.ActionPostCallback = @update_slice;

% Set functions for key presses
hManager = uigetmodemanager(probe_atlas_gui);
[hManager.WindowListenerHandles.Enabled] = deal(false);
set(probe_atlas_gui,'KeyPressFcn',@key_press);
set(probe_atlas_gui,'KeyReleaseFcn',@key_release);

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

% Display the first slice and update the probe position
update_slice(probe_atlas_gui);
update_probe_coordinates(probe_atlas_gui);

% Display controls
display_controls;

end

function key_press(probe_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(probe_atlas_gui);

switch eventdata.Key
    
    case 'uparrow'
        if isempty(eventdata.Modifier)
            % Up: move probe anterior
            ap_offset = -10;
            
            for myshank=1:numel(gui_data.handles.probe_ref_line)
            set(gui_data.handles.probe_ref_line(myshank,1),'XData',get(gui_data.handles.probe_ref_line(myshank,1),'XData') + ap_offset);
            set(gui_data.handles.probe_line(myshank,1),'XData',get(gui_data.handles.probe_line(myshank,1),'XData') + ap_offset);
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-up: increase DV angle
            angle_change=-1;
            update_rotation(probe_atlas_gui,angle_change,'Elevation'); 
            % update gui value
            new_angle = gui_data.probe_angle(2) + angle_change;
            gui_data.probe_angle(2) = new_angle;
            
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % Alt-up: raise probe
            probe_offset = -10;
            
            for myshank=1:numel(gui_data.handles.probe_line)
            old_probe_vector = cell2mat(get(gui_data.handles.probe_line(myshank,1),{'XData','YData','ZData'})');
            
            move_probe_vector = diff(old_probe_vector,[],2)./ ...
                norm(diff(old_probe_vector,[],2))*probe_offset;
            
            new_probe_vector = bsxfun(@plus,old_probe_vector,move_probe_vector);
            
            set(gui_data.handles.probe_line(myshank,1),'XData',new_probe_vector(1,:), ...
                'YData',new_probe_vector(2,:),'ZData',new_probe_vector(3,:)); 
            end
        end
        
    case 'downarrow'
        if isempty(eventdata.Modifier)
            % Down: move probe posterior
            ap_offset = 10;
            for myshank=1:numel(gui_data.handles.probe_ref_line)
            set(gui_data.handles.probe_ref_line(myshank,1),'XData',get(gui_data.handles.probe_ref_line(myshank,1),'XData') + ap_offset);
            set(gui_data.handles.probe_line(myshank,1),'XData',get(gui_data.handles.probe_line(myshank,1),'XData') + ap_offset);
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-down: decrease DV angle
         
            angle_change=1;
            update_rotation(probe_atlas_gui,angle_change,'Elevation');
            new_angle = gui_data.probe_angle(2) + angle_change;
            gui_data.probe_angle(2) = new_angle;
            
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % Alt-down: lower probe
            probe_offset = 10;
            for myshank=1:numel(gui_data.handles.probe_line)  % this could be yoked together 
            old_probe_vector = cell2mat(get(gui_data.handles.probe_line(myshank,1),{'XData','YData','ZData'})');
            
            move_probe_vector = diff(old_probe_vector,[],2)./ ...
                norm(diff(old_probe_vector,[],2))*probe_offset;
            
            new_probe_vector = bsxfun(@plus,old_probe_vector,move_probe_vector);
            
            set(gui_data.handles.probe_line(myshank,1),'XData',new_probe_vector(1,:), ...
                'YData',new_probe_vector(2,:),'ZData',new_probe_vector(3,:)); 
            end 
        end
        
    case 'rightarrow'
        if isempty(eventdata.Modifier)
            % Right: move probe right
            ml_offset = 10;
            for myshank=1:numel(gui_data.handles.probe_ref_line)
            set(gui_data.handles.probe_ref_line(myshank,1),'YData',get(gui_data.handles.probe_ref_line(myshank,1),'YData') + ml_offset);
            set(gui_data.handles.probe_line(myshank,1),'YData',get(gui_data.handles.probe_line(myshank,1),'YData') + ml_offset);
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-right: increase vertical angle
%             angle_change = [1;0];
%             new_angle = gui_data.probe_angle + angle_change;
%             gui_data.probe_angle = new_angle;
%             update_probe_angle(probe_atlas_gui);
            angle_change=-1;
            update_rotation(probe_atlas_gui,angle_change,'Azimuth'); 
            new_angle = gui_data.probe_angle(1) + angle_change;
            gui_data.probe_angle(1) = new_angle;
            
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % alt-right : rotate clockwise 
            angle_change=-1;
            update_rotation(probe_atlas_gui,angle_change,'spin'); 
            new_angle=gui_data.spin+angle_change;
            gui_data.spin=new_angle; 
        end
        
    case 'leftarrow'
        if isempty(eventdata.Modifier)
            % Left: move probe left
            ml_offset = -10;
            for myshank=1:numel(gui_data.handles.probe_ref_line)
            set(gui_data.handles.probe_ref_line(myshank,1),'YData',get(gui_data.handles.probe_ref_line(myshank,1),'YData') + ml_offset);
            set(gui_data.handles.probe_line(myshank,1),'YData',get(gui_data.handles.probe_line(myshank,1),'YData') + ml_offset);
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-right: increase vertical angle
%             angle_change = [-1;0];
%             new_angle = gui_data.probe_angle + angle_change;
%             gui_data.probe_angle = new_angle;
%             update_probe_angle(probe_atlas_gui);
            angle_change=1;
            update_rotation(probe_atlas_gui,angle_change,'Azimuth');
            new_angle = gui_data.probe_angle(1) + angle_change;
            gui_data.probe_angle(1) = new_angle;
            
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % alt-left : rotate anticlockwise 
            angle_change=1;
            update_rotation(probe_atlas_gui,angle_change,'spin');
            new_angle=gui_data.spin+angle_change;
            gui_data.spin=new_angle; 
        end
        
    case 'c'
        % Bring up controls again
        display_controls;
        
    case 'b'
        % Toggle brain outline visibility
        current_visibility = gui_data.handles.cortex_outline.Visible;
        switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
        set(gui_data.handles.cortex_outline,'Visible',new_visibility);
        
    case 'a'
        % Toggle plotted structure visibility
        if ~isempty(gui_data.structure_plot_idx)
            current_visibility = get(gui_data.handles.structure_patch(1),'Visible');
            switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
            set(gui_data.handles.structure_patch,'Visible',new_visibility);
        end
        
    case 's'
        % Toggle slice volume/visibility
        slice_volumes = {'tv','av','none'};        
        new_slice_volume = slice_volumes{circshift( ...
            strcmp(gui_data.handles.slice_volume,slice_volumes),[0,1])};
        
        if strcmp(new_slice_volume,'none')
            set(gui_data.handles.slice_plot,'Visible','off');
        else
            set(gui_data.handles.slice_plot,'Visible','on');
        end
        
        gui_data.handles.slice_volume = new_slice_volume;
        guidata(probe_atlas_gui, gui_data);
        
        update_slice(probe_atlas_gui);
        
    case 'p'
        % Toggle probe visibility
        current_visibility = gui_data.handles.probe_ref_line.Visible;
        switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
        set(gui_data.handles.probe_ref_line,'Visible',new_visibility);
        set(gui_data.handles.probe_line,'Visible',new_visibility);
        
    case 'r'
        % Toggle 3D rotation
        h = rotate3d(gui_data.handles.axes_atlas);
        switch h.Enable
            case 'off'
                h.Enable = 'on';
                % Update the slice whenever a rotation is completed
                h.ActionPostCallback = @update_slice;
                %(need to restore key-press functionality with rotation)
                hManager = uigetmodemanager(probe_atlas_gui);
                [hManager.WindowListenerHandles.Enabled] = deal(false);
                set(probe_atlas_gui,'KeyPressFcn',@key_press);
            case 'on'
                h.Enable = 'off';
        end
        
    case 'm'
        % Set probe angle
        set_probe_position(probe_atlas_gui);
        % Get updated guidata
        gui_data = guidata(probe_atlas_gui);
        
    case {'equal','add'}
        % Add structure(s) to display
        slice_spacing = 10;
        
        % Prompt for which structures to show (only structures which are
        % labelled in the slice-spacing downsampled annotated volume)
        
        if ~any(strcmp(eventdata.Modifier,'shift'))
            % (no shift: list in native CCF order)
            parsed_structures = unique(reshape(gui_data.av(1:slice_spacing:end, ...
                1:slice_spacing:end,1:slice_spacing:end),[],1));
            plot_structures_parsed = listdlg('PromptString','Select a structure to plot:', ...
                'ListString',gui_data.st.safe_name(parsed_structures),'ListSize',[520,500]);
            plot_structures = parsed_structures(plot_structures_parsed);
            
            
            if ~isempty(plot_structures)
                for curr_plot_structure = reshape(plot_structures,1,[])
                    % If this label isn't used, don't plot
                    if ~any(reshape(gui_data.av( ...
                            1:slice_spacing:end,1:slice_spacing:end,1:slice_spacing:end),[],1) == curr_plot_structure)
                        disp(['"' gui_data.st.safe_name{curr_plot_structure} '" is not parsed in the atlas'])
                        continue
                    end
                    
                    gui_data.structure_plot_idx(end+1) = curr_plot_structure;
                    
                    plot_structure_color = hex2dec(reshape(gui_data.st.color_hex_triplet{curr_plot_structure},2,[])')./255;
                    structure_3d = isosurface(permute(gui_data.av(1:slice_spacing:end, ...
                        1:slice_spacing:end,1:slice_spacing:end) == curr_plot_structure,[3,1,2]),0);
                    
                    structure_alpha = 0.2;
                    gui_data.handles.structure_patch(end+1) = patch('Vertices',structure_3d.vertices*slice_spacing, ...
                        'Faces',structure_3d.faces, ...
                        'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
                end
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % (shift: use hierarchy search)
            plot_structures = hierarchicalSelect(gui_data.st);
            
            if ~isempty(plot_structures) % will be empty if dialog was cancelled
                % get all children of this one
                thisID = gui_data.st.id(plot_structures);
                idStr = sprintf('/%d/', thisID);
                theseCh = find(cellfun(@(x)contains(x,idStr), gui_data.st.structure_id_path));
                
                % plot the structure
                slice_spacing = 5;
                plot_structure_color = hex2dec(reshape(gui_data.st.color_hex_triplet{plot_structures},3,[]))./255;
                structure_3d = isosurface(permute(ismember(gui_data.av(1:slice_spacing:end, ...
                    1:slice_spacing:end,1:slice_spacing:end),theseCh),[3,1,2]),0);
                
                structure_alpha = 0.2;
                gui_data.structure_plot_idx(end+1) = plot_structures;
                gui_data.handles.structure_patch(end+1) = patch('Vertices',structure_3d.vertices*slice_spacing, ...
                    'Faces',structure_3d.faces, ...
                    'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
                
            end
            
            
        end
        
    case {'hyphen','subtract'}
        % Remove structure(s) already plotted
        if ~isempty(gui_data.structure_plot_idx)
            remove_structures = listdlg('PromptString','Select a structure to remove:', ...
                'ListString',gui_data.st.safe_name(gui_data.structure_plot_idx));
            delete(gui_data.handles.structure_patch(remove_structures))
            gui_data.structure_plot_idx(remove_structures) = [];
            gui_data.handles.structure_patch(remove_structures) = [];
        end

    case 'x'
        % Export the probe coordinates in Allen CCF to the workspace
        probe_vector = cell2mat(get(gui_data.handles.probe_line,{'XData','YData','ZData'})');
        probe_vector_ccf = round(probe_vector([1,3,2],:))';
        assignin('base','probe_vector_ccf',probe_vector_ccf)
        disp('Copied probe vector coordinates to workspace');
        
    case 'h'
        % Load probe histology points, plot line of best fit
        probe_histology(probe_atlas_gui);
end

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function key_release(probe_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(probe_atlas_gui);

switch eventdata.Key
    case {'rightarrow','leftarrow','uparrow','downarrow'}
        % Update the probe info/slice on arrow release 
        update_probe_coordinates(probe_atlas_gui);
        update_slice(probe_atlas_gui);
end

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function update_slice(probe_atlas_gui,varargin)

% Get guidata
gui_data = guidata(probe_atlas_gui);

% Only update the slice if it's visible
if strcmp(gui_data.handles.slice_plot(1).Visible,'on')
    
    % Get current position of camera
    curr_campos = campos;
    
    % Get probe vector
    probe_ref_top = [gui_data.handles.probe_ref_line(1,1).XData(1), ...
        gui_data.handles.probe_ref_line(1,1).YData(1),gui_data.handles.probe_ref_line(1,1).ZData(1)];
    probe_ref_bottom = [gui_data.handles.probe_ref_line(1).XData(2), ...
        gui_data.handles.probe_ref_line(1,1).YData(2),gui_data.handles.probe_ref_line(1,1).ZData(2)];
    probe_vector = probe_ref_top - probe_ref_bottom;
    
    % Get probe-camera vector
    probe_camera_vector = probe_ref_top - curr_campos;
    
    % Get the vector to plot the plane in (along with probe vector)
    plot_vector = cross(probe_camera_vector,probe_vector);
    
    % Get the normal vector of the plane
    normal_vector = cross(plot_vector,probe_vector);
    
    % Get the plane offset through the probe
    plane_offset = -(normal_vector*probe_ref_top');
    
    % Define a plane of points to index
    % (the plane grid is defined based on the which cardinal plan is most
    % orthogonal to the plotted plane. this is janky but it works)
    slice_px_space = 3;
    %[~,cam_plane] = max(abs((campos - camtarget)./norm(campos - camtarget)));
    
    [~,cam_plane] = max(abs(normal_vector./norm(normal_vector)));
    
    switch cam_plane
        
        case 1
            [plane_y,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,3),1:slice_px_space:size(gui_data.tv,2));
            plane_x = ...
                (normal_vector(2)*plane_y+normal_vector(3)*plane_z + plane_offset)/ ...
                -normal_vector(1);
            
        case 2
            [plane_x,plane_z] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,2));
            plane_y = ...
                (normal_vector(1)*plane_x+normal_vector(3)*plane_z + plane_offset)/ ...
                -normal_vector(2);
            
        case 3
            [plane_x,plane_y] = meshgrid(1:slice_px_space:size(gui_data.tv,1),1:slice_px_space:size(gui_data.tv,3));
            plane_z = ...
                (normal_vector(1)*plane_x+normal_vector(2)*plane_y + plane_offset)/ ...
                -normal_vector(3);
            
    end
    
    % Get the coordiates on the plane
    x_idx = round(plane_x);
    y_idx = round(plane_y);
    z_idx = round(plane_z);
    
    % Find plane coordinates in bounds with the volume
    use_xd = x_idx > 0 & x_idx < size(gui_data.tv,1);
    use_yd = y_idx > 0 & y_idx < size(gui_data.tv,3);
    use_zd = z_idx > 0 & z_idx < size(gui_data.tv,2);
    use_idx = use_xd & use_yd & use_zd;
    
    curr_slice_idx = sub2ind(size(gui_data.tv),x_idx(use_idx),z_idx(use_idx),y_idx(use_idx));
    
    % Find plane coordinates that contain brain
    curr_slice_isbrain = false(size(use_idx));
    curr_slice_isbrain(use_idx) = gui_data.av(curr_slice_idx) > 1;
    
    % Index coordinates in bounds + with brain
    grab_pix_idx = sub2ind(size(gui_data.tv),x_idx(curr_slice_isbrain),z_idx(curr_slice_isbrain),y_idx(curr_slice_isbrain));
    
    % Grab pixels from (selected) volume
    curr_slice = nan(size(use_idx));
    switch gui_data.handles.slice_volume
        case 'tv'
            curr_slice(curr_slice_isbrain) = gui_data.tv(grab_pix_idx);
            colormap(gui_data.handles.axes_atlas,'gray');
            caxis([0,255]);
        case 'av'
            curr_slice(curr_slice_isbrain) = gui_data.av(grab_pix_idx);
            colormap(gui_data.handles.axes_atlas,gui_data.cmap);
            caxis([1,size(gui_data.cmap,1)]);
    end
    
    % Update the slice display
    set(gui_data.handles.slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);
    
    % Upload gui_data
    guidata(probe_atlas_gui, gui_data);
    
end

end

function set_probe_position(probe_atlas_gui,varargin)

% Get guidata
gui_data = guidata(probe_atlas_gui);

% Prompt for angles
prompt_text = { ...
    'AP position (mm from bregma)', ...
    'ML position (mm from bregma)', ...
    'ML angle (relative to midline)', ....
    'DV angle (relative to coronal cut)',...
    'probe spin'};

new_probe_position = cellfun(@str2num,inputdlg(prompt_text,'Set probe position',1));

% Convert probe position: mm->CCF and degrees->radians
probe_ccf_coordinates = round(gui_data.bregma([1,3])' - new_probe_position(1:2)*100);
%probe_angle_rad = (new_probe_position(3:4)/360)*2*pi;

% Update the probe and trajectory reference
[ap_max,dv_max,ml_max] = size(gui_data.tv);

max_ref_length = sqrt(sum(([ap_max,dv_max,ml_max].^2)));

%[~,y,z] = sph2cart(probe_angle_rad(1),probe_angle_rad(2),max_ref_length);

%[~,x,~] = sph2cart(probe_angle_rad(2),probe_angle_rad(1),max_ref_length);

probe_ref_top = [probe_ccf_coordinates(1),probe_ccf_coordinates(2),0];
probe_ref_bottom = [probe_ccf_coordinates(1),probe_ccf_coordinates(2),max_ref_length];
probe_ref_vector = [probe_ref_top;probe_ref_bottom]';
probe_vector = [probe_ref_vector(:,1),diff(probe_ref_vector,[],2)./ ...
    norm(diff(probe_ref_vector,[],2))*gui_data.probe_length + probe_ref_vector(:,1)];

 
shank_spacing_vector=reshape([gui_data.shank_spacing;gui_data.shank_spacing],[],1)'; 

probe_ref_vector=repmat(probe_ref_vector,1,gui_data.shank_no);
probe_ref_vector(2,:)=probe_ref_vector(2,:)+shank_spacing_vector;

probe_vector=repmat(probe_vector,1,gui_data.shank_no);
probe_vector(2,:)=probe_vector(2,:)+shank_spacing_vector;


gui_data.probe_angle(1) =new_probe_position(3);
gui_data.probe_angle(2) =new_probe_position(4);
gui_data.spin=new_probe_position(5); 
% % spin
% rotation_matrix=get_rotation_matrix(gui_data.spin,'spin');
% [probe_ref_vector,probe_vector]=please_rotate(probe_ref_vector,probe_vector,rotation_matrix);
% 
% % rotate ML 
% rotation_matrix=get_rotation_matrix(gui_data.probe_angle(1),'ML');
% [probe_ref_vector,probe_vector]=please_rotate(probe_ref_vector,probe_vector,rotation_matrix);
% 
% % rotate AP 
% rotation_matrix=get_rotation_matrix(gui_data.probe_angle(2),'AP');
% [probe_ref_vector,probe_vector]=please_rotate(probe_ref_vector,probe_vector,rotation_matrix);



for myshank=1:gui_data.shank_no
    st_ix=2*myshank-1;
    set(gui_data.handles.probe_ref_line(myshank),'XData',probe_ref_vector(1,st_ix:(st_ix+1)), ...
            'YData',probe_ref_vector(2,st_ix:(st_ix+1)), ...
            'ZData',probe_ref_vector(3,st_ix:(st_ix+1)));

    set(gui_data.handles.probe_line(myshank),'XData',probe_vector(1,st_ix:(st_ix+1)), ...
        'YData',probe_vector(2,st_ix:(st_ix+1)), ...
        'ZData',probe_vector(3,st_ix:(st_ix+1)));    
end 

% Upload gui_data
guidata(probe_atlas_gui, gui_data);

% Update the slice and probe coordinates
update_slice(probe_atlas_gui);
update_probe_coordinates(probe_atlas_gui);

end

function update_probe_angle(probe_atlas_gui,varargin)
% not used right now 
% Get guidata
gui_data = guidata(probe_atlas_gui);

% Get the positions of the probe and trajectory reference
for myshank=1:numel(gui_data.handles.probe_ref_line)
    probe_ref_vector = cell2mat(get(gui_data.handles.probe_ref_line(myshank,1),{'XData','YData','ZData'})');
    probe_vector = cell2mat(get(gui_data.handles.probe_line(myshank,1),{'XData','YData','ZData'})');

    % Update the probe trajectory reference (rotate about trajectory top)
    [ap_max,dv_max,ml_max] = size(gui_data.tv);

    max_ref_length = sqrt(sum(([ap_max,dv_max,ml_max].^2)));

    probe_angle_rad = (gui_data.probe_angle./360)*2*pi;
    [x,y,z] = sph2cart(pi-probe_angle_rad(1),probe_angle_rad(2),max_ref_length);

    new_probe_ref_top = [probe_ref_vector(1,1),probe_ref_vector(2,1),0];
    new_probe_ref_bottom = new_probe_ref_top + [x,y,z];
    new_probe_ref_vector = [new_probe_ref_top;new_probe_ref_bottom]';

    set(gui_data.handles.probe_ref_line(myshank,1),'XData',new_probe_ref_vector(1,:), ...
        'YData',new_probe_ref_vector(2,:), ...
        'ZData',new_probe_ref_vector(3,:));

    % Update probe (retain depth)
    new_probe_vector = [new_probe_ref_vector(:,1),diff(new_probe_ref_vector,[],2)./ ...
        norm(diff(new_probe_ref_vector,[],2))*gui_data.probe_length + new_probe_ref_vector(:,1)];

    probe_depth = sqrt(sum((probe_ref_vector(:,1) - probe_vector(:,1)).^2));
    new_probe_vector_depth = (diff(new_probe_vector,[],2)./ ...
        norm(diff(new_probe_vector,[],2))*probe_depth) + new_probe_vector;

    set(gui_data.handles.probe_line(myshank,1),'XData',new_probe_vector_depth(1,:), ...
        'YData',new_probe_vector_depth(2,:),'ZData',new_probe_vector_depth(3,:));
end
% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function update_probe_coordinates(probe_atlas_gui,varargin)

% Get guidata
gui_data = guidata(probe_atlas_gui);

for myshank=1:numel(gui_data.handles.probe_ref_line)
    % Get the positions of the probe and trajectory reference
    probe_ref_vector = cell2mat(get(gui_data.handles.probe_ref_line(myshank,1),{'XData','YData','ZData'})');
    probe_vector = cell2mat(get(gui_data.handles.probe_line(myshank,1),{'XData','YData','ZData'})');

    trajectory_n_coords = max(abs(diff(probe_ref_vector,[],2)));
    [trajectory_xcoords,trajectory_ycoords,trajectory_zcoords] = deal( ...
        linspace(probe_ref_vector(1,1),probe_ref_vector(1,2),trajectory_n_coords), ...
        linspace(probe_ref_vector(2,1),probe_ref_vector(2,2),trajectory_n_coords), ...
        linspace(probe_ref_vector(3,1),probe_ref_vector(3,2),trajectory_n_coords));

    probe_n_coords = sqrt(sum(diff(probe_vector,[],2).^2));
    [probe_xcoords,probe_ycoords,probe_zcoords] = deal( ...
        linspace(probe_vector(1,1),probe_vector(1,2),probe_n_coords), ...
        linspace(probe_vector(2,1),probe_vector(2,2),probe_n_coords), ...
        linspace(probe_vector(3,1),probe_vector(3,2),probe_n_coords));

    % Get brain labels across the probe and trajectory, and intersection with brain
    pixel_space = 5;
    trajectory_areas = interp3(single(gui_data.av(1:pixel_space:end,1:pixel_space:end,1:pixel_space:end)), ...
        round(trajectory_zcoords/pixel_space),round(trajectory_xcoords/pixel_space),round(trajectory_ycoords/pixel_space),'nearest');
    trajectory_brain_idx = find(trajectory_areas > 1,1);
    trajectory_brain_intersect = ...
        [trajectory_xcoords(trajectory_brain_idx),trajectory_ycoords(trajectory_brain_idx),trajectory_zcoords(trajectory_brain_idx)]';

    probe_areas = interp3(single(gui_data.av(1:pixel_space:end,1:pixel_space:end,1:pixel_space:end)), ...
        round(probe_zcoords/pixel_space),round(probe_xcoords/pixel_space),round(probe_ycoords/pixel_space),'nearest')';
    probe_area_boundaries = intersect(unique([find(~isnan(probe_areas),1,'first'); ...
        find(diff(probe_areas) ~= 0);find(~isnan(probe_areas),1,'last')]),find(~isnan(probe_areas)));
    probe_area_centers = probe_area_boundaries(1:end-1) + diff(probe_area_boundaries)/2;
    probe_area_labels = gui_data.st.safe_name(probe_areas(round(probe_area_centers)));

    % Get position of brain intersect relative to bregma
    probe_bregma_coordinate = round((gui_data.bregma([1,3])' - trajectory_brain_intersect(1:2))*10);

    % Get the depth of the bottom of the probe (sign: hack by z offset)
    probe_depth = round(sqrt(sum((trajectory_brain_intersect - probe_vector(:,2)).^2))*10)* ...
        sign(probe_vector(3,2)-trajectory_brain_intersect(3));

    % Update the text
    if myshank==1 
        probe_text = ['Shank 1 insertion: ' ....
            num2str(probe_bregma_coordinate(1)) ' AP, ', ...
            num2str(-probe_bregma_coordinate(2)) ' ML, ', ...
            num2str(probe_depth) ' Probe-axis, ' ...
            num2str(gui_data.probe_angle(1)+90) char(176) ' from midline, ' ...
            num2str(-(gui_data.probe_angle(2)+90)) char(176) ' from horizontal,' ...
            num2str(gui_data.spin) char(176) ' spin '];
        set(gui_data.probe_coordinates_text,'String',probe_text);
    end

    % Update the probe areas
    yyaxis(gui_data.handles.axes_probe_areas(myshank),'right');
    set(gui_data.handles.probe_areas_plot(myshank),'YData',[1:length(probe_areas)]*10,'CData',probe_areas); 
    set(gui_data.handles.axes_probe_areas(myshank),'YTick',probe_area_centers*10,'YTickLabels',probe_area_labels','YTickLabelRotation',45);

end
% Upload gui_data
guidata(probe_atlas_gui, gui_data);

end

function probe_histology(probe_atlas_gui)
% Load histology points
% UNDER CONSTRUCTION
% (used to use SHARP-Track, now using mine)

% Get guidata
gui_data = guidata(probe_atlas_gui);

[probe_file,probe_path] = uigetfile('*.mat','Choose probe coordinate file');
load([probe_path,probe_file]);

if exist('pointList','var')
    histology_points = pointList.pointList{1};
elseif exist('probe_ccf','var')
    histology_points = probe_ccf(1).points; % only use first probe
end

r0 = mean(histology_points,1);
xyz = bsxfun(@minus,histology_points,r0);
[~,~,V] = svd(xyz,0);
histology_probe_direction = V(:,1);

probe_eval_points = [-1000,1000];
probe_line_endpoints = bsxfun(@plus,bsxfun(@times,probe_eval_points',histology_probe_direction'),r0);

% Philip's GUI: not saved in native CCF order?
% plot3(histology_points(:,3),histology_points(:,1),histology_points(:,2),'.b','MarkerSize',20);
% line(P(:,3),P(:,1),P(:,2),'color','k','linewidth',2)

% % Mine: saved in native CCF order [AP,DV,ML]
plot3(histology_points(:,1),histology_points(:,3),histology_points(:,2),'.b','MarkerSize',20);
line(probe_line_endpoints(:,1),probe_line_endpoints(:,3),probe_line_endpoints(:,2),'color','k','linewidth',2)

% Place the probe on the histology best-fit axis
[ap_max,dv_max,ml_max] = size(gui_data.tv);

probe_ref_top = probe_line_endpoints(1,[1,3,2]);
probe_ref_bottom = probe_line_endpoints(2,[1,3,2]);
probe_ref_vector = [probe_ref_top;probe_ref_bottom]';

set(gui_data.handles.probe_ref_line,'XData',probe_ref_vector(1,:), ...
    'YData',probe_ref_vector(2,:), ...
    'ZData',probe_ref_vector(3,:));

probe_vector = [probe_ref_vector(:,1),diff(probe_ref_vector,[],2)./ ...
    norm(diff(probe_ref_vector,[],2))*gui_data.probe_length + probe_ref_vector(:,1)];
set(gui_data.handles.probe_line,'XData',probe_vector(1,:), ...
    'YData',probe_vector(2,:),'ZData',probe_vector(3,:));

% Upload gui_data
[theta,phi] = cart2sph(diff(probe_ref_vector(1,:)),diff(probe_ref_vector(2,:)),diff(probe_ref_vector(3,:)));
gui_data.probe_angle = ([theta,phi]/(2*pi))*360;
guidata(probe_atlas_gui, gui_data);

% Update the slice and probe coordinates
update_slice(probe_atlas_gui);
update_probe_coordinates(probe_atlas_gui);
end

function display_controls

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Probe: \rm' ...
    'Arrow keys : translate probe' ...
    'Alt/Option up/down : raise/lower probe' ...
    'Alt/Option left/right : spin probe' ...
    'Shift arrow keys : rotate probe (around top)' ...
    'm : set probe location manually', ...
    '\bf 3D brain areas: \rm' ...
    ' =/+ : add (list selector)' ...
    ' Shift =/+ : add (hierarchy selector)' ...
    ' - : remove', ...
    '\bf Visibility: \rm' ...
    's : atlas slice (toggle tv/av/off)' ...
    'b : brain outline' ...
    'p : probe' ...
    'a : 3D brain areas' ...
    '\bf Other: \rm' ...
    'r : toggle clickable rotation' ...
    'x : export probe coordinates to workspace' ...
    'h : load and plot histology-defined trajectory', ...
    'c : bring up controls box'}, ...
    'Controls',CreateStruct);

end

function [rotation_matrix]=get_rotation_matrix(probe_vector,change_angle,current_angle,direction)
% Zhiwen modified 
% Convert Axis-Angle Rotation to Rotation Matrix
angle_rad=deg2rad(change_angle); 
current_angle_rad = deg2rad(current_angle); 
probeV = probe_vector(:,2)-probe_vector(:,1);
probehV = probe_vector(:,3)-probe_vector(:,1);
    switch direction
        case 'Azimuth'        
            % Azimuth
            % Convert Axis-Angle Rotation to Rotation Matrix
            % DV line as rotation axis
            axang = [0 0 288.8 angle_rad];
        case 'Elevation'
            % Elevation
            % Convert Axis-Angle Rotation to Rotation Matrix
            % plane of 4-shank probe is the plane to rotate on;
            % which means it rotates around vector orthogonal to this
            % plane;
            axang1 = [0 0 1 current_angle_rad];
            rotation_matrix1 = axang2rotm1(axang1);
            rotated_azimuth = rotation_matrix1*[1; 0; 0];
            %%
            othorg = cross(probeV,rotated_azimuth);           
            axang = [othorg(1) othorg(2) othorg(3) angle_rad];
        case 'spin'
            % spin 
            % Convert Axis-Angle Rotation to Rotation Matrix
            % spin around reference vector;
            axang = [probeV(1,1) probeV(2,1) probeV(3,1) angle_rad];
    end 
    rotation_matrix = axang2rotm1(axang);
end


function [rotated_ref_bregma,rotated_probe_bregma]=please_rotate(probe_ref_vector,probe_vector,rotation_matrix)
    zeropoint=probe_vector(:,1);
    zerocentered_ref=probe_ref_vector-zeropoint;


    zerocentered_probe=probe_vector-zeropoint;


    rotated_ref=rotation_matrix*zerocentered_ref; 
    rotated_ref_bregma=rotated_ref+zeropoint; 

    rotated_probe=rotation_matrix*zerocentered_probe; 
    rotated_probe_bregma=rotated_probe+zeropoint;
end 

function update_rotation(probe_atlas_gui,angle_change,mydir)
% this is the function when we update the rotation of the probe either in
% ML or AP 
    gui_data = guidata(probe_atlas_gui);
    probe_ref_vector = cell2mat(get(gui_data.handles.probe_ref_line,{'XData','YData','ZData'})');
    probe_vector = cell2mat(get(gui_data.handles.probe_line,{'XData','YData','ZData'})');
    current_angle = gui_data.probe_angle(1)+90;
   
    rotation_matrix=get_rotation_matrix(probe_vector,angle_change,current_angle,mydir);
    [rotated_ref_bregma,rotated_probe_bregma]=please_rotate(probe_ref_vector,probe_vector,rotation_matrix);

    for myshank=1:gui_data.shank_no
        st_ix=2*myshank-1;
        set(gui_data.handles.probe_ref_line(myshank),'XData',rotated_ref_bregma(1,st_ix:(st_ix+1)), ...
                'YData',rotated_ref_bregma(2,st_ix:(st_ix+1)), ...
                'ZData',rotated_ref_bregma(3,st_ix:(st_ix+1)));

        set(gui_data.handles.probe_line(myshank),'XData',rotated_probe_bregma(1,st_ix:(st_ix+1)), ...
            'YData',rotated_probe_bregma(2,st_ix:(st_ix+1)), ...
            'ZData',rotated_probe_bregma(3,st_ix:(st_ix+1)));    
    end 
end

function R = axang2rotm1( axang )
%AXANG2ROTM Convert axis-angle rotation representation to rotation matrix
%   R = AXANG2ROTM(AXANG) converts a 3D rotation given in axis-angle form,
%   AXANG, to an orthonormal rotation matrix, R. AXANG is an N-by-4
%   matrix of N axis-angle rotations. The first three elements of every
%   row specify the rotation axis and the last element defines the rotation
%   angle (in radians).
%   The output, R, is an 3-by-3-by-N matrix containing N rotation matrices. 
%   Each rotation matrix has a size of 3-by-3 and is orthonormal.
%
%   Example:
%      % Convert a rotation from axis-angle to rotation matrix
%      axang = [0 1 0 pi/2];
%      R = axang2rotm(axang)

% For a single axis-angle vector [ax ay az theta] the output rotation
% matrix R can be computed as follows:
% R =  [t*x*x + c	  t*x*y - z*s	   t*x*z + y*s
%       t*x*y + z*s	  t*y*y + c	       t*y*z - x*s
%       t*x*z - y*s	  t*y*z + x*s	   t*z*z + c]
% where,
% c = cos(theta)
% s = sin(theta)
% t = 1 - c
% x = normalized axis ax coordinate
% y = normalized axis ay coordinate
% z = normalized axis az coordinate

% Normalize the axis
v = bsxfun(@times, axang(:,1:3), 1./sqrt(sum(axang(:,1:3).^2,2)));
% Extract the rotation angles and shape them in depth dimension
numInputs = size(axang,1);
theta = zeros(1,1,numInputs);
theta(1,1,:) = axang(:,4);

% Compute rotation matrices
cth = cos(theta);
sth = sin(theta);
vth = (1 - cth);

% Preallocate input vectors
vx = zeros(1,1,numInputs,'like',axang);
vy = vx;
vz = vx;

% Shape input vectors in depth dimension
vx(1,1,:) = v(:,1);
vy(1,1,:) = v(:,2);
vz(1,1,:) = v(:,3);

% Explicitly specify concatenation dimension
tempR = cat(1, vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth, ...
               vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth, ...
               vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth);

R = reshape(tempR, [3, 3, length(vx)]);
R = permute(R, [2 1 3]);

end

