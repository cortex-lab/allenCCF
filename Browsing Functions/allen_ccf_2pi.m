function allen_ccf_2pi(tv,av,st)
% allen_ccf_2pi(tv,av,st)
% written by Samuel Picard (samuel.picard@ucl.ac.uk)
% based on original allen_ccf_npx tool written by Andy Peters (peters.andrew.j@gmail.com)
%
% GUI for planning 2pi chronic window implant with the Allen CCF
% Part of repository: https://github.com/cortex-lab/allenCCF
% - directions for installing atlas in repository readme
% - some dependent functions from that repository
%
% (optional inputs - if CCF path written in Line 22, loaded automatically)
% tv, av, st = CCF template volume, annotated volume, structure tree
% tv = readNPY('template_volume_10um.npy');
% av = readNPY('annotation_volume_10um_by_index.npy');
% st = loadStructureTree('structure_tree_safe_2017.csv');

% Initialize gui_data structure
gui_data = struct;

% Allen CCF-bregma transform (estimated from eyeballing Paxinos->CCF)
% [AP,DV,ML]
bregma = [540,0,570];

% If not already loaded in, load in atlas
if nargin < 3
    allen_atlas_path = 'C:\Users\Samuel\Documents\GitHub\allenCCF'; %put path to atlas files here
    if isempty(allen_atlas_path)
        error('Enter path where Allen CCF is stored at Line 26');
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
window_atlas_gui = figure('Toolbar','none','Menubar','none','color','w', ...
    'Name','Atlas-window viewer','Units','normalized','Position',[0.2,0.2,0.7,0.7]);

% Set up the atlas axes
axes_atlas = subplot(2,2,[1 3]);
[~, brain_outline] = plotBrainGrid([],axes_atlas);
hold(axes_atlas,'on');
axis vis3d equal off manual
view([-30,25]);
caxis([0 300]);
[ap_max,dv_max,ml_max] = size(tv);
xlim([-10,ap_max+10])
ylim([-10,ml_max+10])
zlim([-10,dv_max+10])

% Set up the window area axes
%top of window
axes_window_areas = subplot(2,2,2);
axes_window_areas.ActivePositionProperty = 'position';
window_areas_plot = image(nan(166,166)); %initialized for 5mm window at 30um resolution
axis off auto
window_areas_text = text(axes_window_areas,zeros(1,100),zeros(1,100),repmat({''}, 1, 100));
set(axes_window_areas,'FontSize',11,'XLimSpec', 'Tight','YLimSpec', 'Tight');
colormap(axes_window_areas,cmap);
caxis([1,size(cmap,1)])
%bottom of window
axes_window_areas_bottom = subplot(2,2,4);
axes_window_areas_bottom.ActivePositionProperty = 'position';
window_areas_bottom_plot = image(nan(166,166)); %initialized for 5mm window at 30um resolution
axis off auto
window_areas_bottom_text = text(axes_window_areas_bottom,zeros(1,100),zeros(1,100),repmat({''}, 1, 100));
set(axes_window_areas_bottom,'FontSize',11,'XLimSpec', 'Tight','YLimSpec', 'Tight');
colormap(axes_window_areas_bottom,cmap);
caxis([1,size(cmap,1)])

% Position the axes
set(axes_atlas,'Position',[-0.15,-0.1,1,1.2]);
set(axes_window_areas,'Position',[0.7,0.55,0.2,0.34]);
set(axes_window_areas_bottom,'Position',[0.7,0.1,0.2,0.34]);

% Set the current axes to the atlas (dirty, but some gca requirements)
axes(axes_atlas);

% Set up the window reference/actual
window_centre = [bregma(1),bregma(3),0];
window_depth = 60; % effective z-range of 2p microscope (in 10 ums)
window_diameter = 500; % diameter in 10 ums
window_vector = [window_centre',[window_centre(1),window_centre(2),window_depth]'];
window_vector_line = line(window_vector(1,:),window_vector(2,:),window_vector(3,:),'linewidth',3,'color','b','linestyle','-');
window_rads = linspace(-pi,pi,100);
window_circle = line(0.5*window_diameter*cos(window_rads)+window_centre(1),...
    0.5*window_diameter*sin(window_rads)+window_centre(2),...
    zeros(size(window_rads))+window_centre(3),...
    'linewidth',1.5,'color','r','linestyle','--');

% add titles for window plots
title(axes_window_areas,'Top of window');
title(axes_window_areas_bottom,[sprintf('Bottom (%i ',10*window_depth),'\mu','m)']);

% Set up the text to display coordinates
window_coordinates_text = uicontrol('Style','text','String','', ...
    'Units','normalized','Position',[0,0.95,1,0.05], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontSize',12);

% Store data
gui_data.tv = tv; % Intensity atlas
gui_data.av = av; % Annotated atlas
gui_data.st = st; % Labels table
gui_data.cmap = cmap; % Atlas colormap
gui_data.bregma = bregma; % Bregma for external referencing
gui_data.window_depth = window_depth; % Effective z-depth under window
gui_data.structure_plot_idx = []; % Plotted structures
gui_data.window_angle = [0;0]; % window angles in ML/DV
gui_data.window_centre = window_centre; % Window reference centre on 3D atlas
gui_data.window_diameter = window_diameter;
gui_data.ref_centre = window_centre;

%Store handles
gui_data.handles.cortex_outline = brain_outline;
gui_data.handles.structure_patch = []; % Plotted structures
gui_data.handles.axes_atlas = axes_atlas; % Axes with 3D atlas
gui_data.handles.axes_window_areas = axes_window_areas; % Axes with window areas
gui_data.handles.axes_window_areas_bottom = axes_window_areas_bottom; % Axes with window areas
gui_data.handles.slice_plot = surface('EdgeColor','none'); % Slice on 3D atlas
gui_data.handles.slice_volume = 'tv'; % The volume shown in the slice
gui_data.handles.window_vector = window_vector_line; % Window centre vector on 3D atlas
gui_data.handles.window_circle = window_circle; % Window circle on 3D atlas
gui_data.handles.window_areas_plot = window_areas_plot; % Color-coded window regions
gui_data.handles.window_areas_bottom_plot = window_areas_bottom_plot; % Color-coded window regions
gui_data.handles.window_areas_text = window_areas_text; % Labels for window regions
gui_data.handles.window_areas_bottom_text = window_areas_bottom_text; % Labels for window regions
gui_data.window_coordinates_text = window_coordinates_text; % Window coordinates text

% Make 3D rotation the default state (toggle on/off with 'r')
h = rotate3d(axes_atlas);
h.Enable = 'on';
% Update the slice whenever a rotation is completed
%h.ActionPostCallback = @update_slice;

% Set functions for key presses
hManager = uigetmodemanager(window_atlas_gui);
[hManager.WindowListenerHandles.Enabled] = deal(false);
set(window_atlas_gui,'KeyPressFcn',@key_press);
set(window_atlas_gui,'KeyReleaseFcn',@key_release);

% Upload gui_data
guidata(window_atlas_gui, gui_data);

% Display the first slice and update the window position
update_slice(window_atlas_gui);
update_window_coordinates(window_atlas_gui);

% Display controls
display_controls;

end

function key_press(window_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(window_atlas_gui);

switch eventdata.Key
    
    case 'uparrow'
        if isempty(eventdata.Modifier)
            % Up: move window anterior
            ap_offset = -10;
            set(gui_data.handles.window_circle,'XData',get(gui_data.handles.window_circle,'XData') + ap_offset);
            set(gui_data.handles.window_vector,'XData',get(gui_data.handles.window_vector,'XData') + ap_offset);
            gui_data.window_centre = gui_data.window_centre + [ap_offset,0,0];
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-up: increase DV angle
            angle_change = [1;0];
            gui_data = update_window_angle(window_atlas_gui,angle_change);
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % Alt-up: raise window
            dv_offset = -10;
            old_window_vector = cell2mat(get(gui_data.handles.window_vector,{'XData','YData','ZData'})');
            old_window_circle = cell2mat(get(gui_data.handles.window_circle,{'XData','YData','ZData'})');
            move_vector = diff(old_window_vector,[],2)./ ...
                norm(diff(old_window_vector,[],2))*dv_offset;
            new_window_vector = bsxfun(@plus,old_window_vector,move_vector);
            new_window_circle = bsxfun(@plus,old_window_circle,move_vector);
            set(gui_data.handles.window_vector,'XData',new_window_vector(1,:), ...
                'YData',new_window_vector(2,:),'ZData',new_window_vector(3,:));
            set(gui_data.handles.window_circle,'XData',new_window_circle(1,:), ...
                'YData',new_window_circle(2,:),'ZData',new_window_circle(3,:));
            gui_data.window_centre = gui_data.window_centre + move_vector';
            
        end
        
    case 'downarrow'
        if isempty(eventdata.Modifier)
            % Down: move window posterior
            ap_offset = 10;
            set(gui_data.handles.window_circle,'XData',get(gui_data.handles.window_circle,'XData') + ap_offset);
            set(gui_data.handles.window_vector,'XData',get(gui_data.handles.window_vector,'XData') + ap_offset);
            gui_data.window_centre = gui_data.window_centre + [ap_offset,0,0];
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-down: decrease DV angle
            angle_change = [-1;0];
            gui_data = update_window_angle(window_atlas_gui,angle_change);
        elseif any(strcmp(eventdata.Modifier,'alt'))
            % Alt-down: lower window
            dv_offset = 10;
            old_window_vector = cell2mat(get(gui_data.handles.window_vector,{'XData','YData','ZData'})');
            old_window_circle = cell2mat(get(gui_data.handles.window_circle,{'XData','YData','ZData'})');
            move_vector = diff(old_window_vector,[],2)./ ...
                norm(diff(old_window_vector,[],2))*dv_offset;
            new_window_vector = bsxfun(@plus,old_window_vector,move_vector);
            new_window_circle = bsxfun(@plus,old_window_circle,move_vector);
            set(gui_data.handles.window_vector,'XData',new_window_vector(1,:), ...
                'YData',new_window_vector(2,:),'ZData',new_window_vector(3,:));
            set(gui_data.handles.window_circle,'XData',new_window_circle(1,:), ...
                'YData',new_window_circle(2,:),'ZData',new_window_circle(3,:));
            gui_data.window_centre = gui_data.window_centre + move_vector';
            
        end
        
    case 'rightarrow'
        if isempty(eventdata.Modifier)
            % Right: move window right
            ml_offset = 10;
            set(gui_data.handles.window_circle,'YData',get(gui_data.handles.window_circle,'YData') + ml_offset);
            set(gui_data.handles.window_vector,'YData',get(gui_data.handles.window_vector,'YData') + ml_offset);
            gui_data.window_centre = gui_data.window_centre + [0,ml_offset,0];
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-right: increase vertical angle
            angle_change = [0;1];
            gui_data = update_window_angle(window_atlas_gui,angle_change);
        end
        
    case 'leftarrow'
        if isempty(eventdata.Modifier)
            % Left: move window left
            ml_offset = -10;
            set(gui_data.handles.window_circle,'YData',get(gui_data.handles.window_circle,'YData') + ml_offset);
            set(gui_data.handles.window_vector,'YData',get(gui_data.handles.window_vector,'YData') + ml_offset);
            gui_data.window_centre = gui_data.window_centre + [0,ml_offset,0];
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % Ctrl-left: decrease vertical angle
            angle_change = [0;-1];
            gui_data = update_window_angle(window_atlas_gui,angle_change);
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
            current_alpha = get(gui_data.handles.structure_patch(1),'FaceAlpha');
            switch current_alpha
                case 0
                    new_alpha = 0.2;
                case 0.2
                    new_alpha = 1;
                case 1
                    new_alpha = 0;
            end
            set(gui_data.handles.structure_patch,'FaceAlpha',new_alpha);
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
        guidata(window_atlas_gui, gui_data);
        
        update_slice(window_atlas_gui);
        
    case 'w'
        % Toggle window visibility
        current_visibility = gui_data.handles.window_circle.Visible;
        switch current_visibility; case 'on'; new_visibility = 'off'; case 'off'; new_visibility = 'on'; end;
        set(gui_data.handles.window_circle,'Visible',new_visibility);
        
    case 'm'
        % Set window position manually and find window angle automatically
        set_window_position(window_atlas_gui);
        % Get updated guidata
        gui_data = guidata(window_atlas_gui);
        
    case {'equal','add'}
        % Add structure(s) to display
        slice_spacing = 10;
        
        % Prompt for which structures to show (only structures which are
        % labelled in the slice-spacing downsampled annotated volume)
        
        if ~any(strcmp(eventdata.Modifier,'shift'))
            % (no shift: list in native CCF order)
            
            parsed_structures = unique(reshape(gui_data.av(1:slice_spacing:end, ...
                1:slice_spacing:end,1:slice_spacing:end),[],1));
            if ~any(strcmp(eventdata.Modifier,'alt'))
                % (no alt: list all)
                plot_structures_parsed = listdlg('PromptString','Select a structure to plot:', ...
                    'ListString',gui_data.st.safe_name(parsed_structures),'ListSize',[520,500]);
                plot_structures = parsed_structures(plot_structures_parsed);
            else
                % (alt: search list)
                structure_search = lower(inputdlg('Search structures'));
                structure_match = find(contains(lower(gui_data.st.safe_name),structure_search));
                list_structures = intersect(parsed_structures,structure_match);
                if isempty(list_structures)
                    error('No structure search results')
                end
                
                plot_structures_parsed = listdlg('PromptString','Select a structure to plot:', ...
                    'ListString',gui_data.st.safe_name(list_structures),'ListSize',[520,500]);
                plot_structures = list_structures(plot_structures_parsed);
            end
            
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
                    
                    if isempty(gui_data.handles.structure_patch)
                        structure_alpha = 0.2;
                    else
                        structure_alpha = get(gui_data.handles.structure_patch(1),'FaceAlpha');
                    end
                    gui_data.handles.structure_patch(end+1) = patch('Vertices',structure_3d.vertices*slice_spacing, ...
                        'Faces',structure_3d.faces, ...
                        'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
                end
            end
            
        elseif any(strcmp(eventdata.Modifier,'shift'))
            % (shift: use hierarchy search)
            if ~any(strcmp(eventdata.Modifier,'control'))
                % (no ctrl, choose between all structures)
                plot_structures = hierarchicalSelect(gui_data.st);
            else
                % (ctrl: take pre-defined regions: isocortex, cerebellum, SC & IC)
                plot_structures = [7,13,31,108,115,122,158,221,239,246,253,279,298,340,361,368,374,381,809,813,1016,1017];
            end
            
            if ~isempty(plot_structures) % will be empty if dialog was cancelled
                % get all children of this one
                for plot_structure = plot_structures
                    thisID = gui_data.st.id(plot_structure);
                    idStr = sprintf('/%d/', thisID);
                    theseCh = find(cellfun(@(x)contains(x,idStr), gui_data.st.structure_id_path));
                    
                    % plot the structure
                    slice_spacing = 5;
                    plot_structure_color = hex2dec(reshape(gui_data.st.color_hex_triplet{plot_structure},3,[]))./255;
                    structure_3d = isosurface(permute(ismember(gui_data.av(1:slice_spacing:end, ...
                        1:slice_spacing:end,1:slice_spacing:end),theseCh),[3,1,2]),0);
                    
                    if isempty(gui_data.handles.structure_patch)
                        structure_alpha = 0.2;
                    else
                        structure_alpha = get(gui_data.handles.structure_patch(1),'FaceAlpha');
                    end
                    gui_data.structure_plot_idx(end+1) = plot_structure;
                    gui_data.handles.structure_patch(end+1) = patch('Vertices',structure_3d.vertices*slice_spacing, ...
                        'Faces',structure_3d.faces, ...
                        'FaceColor',plot_structure_color,'EdgeColor','none','FaceAlpha',structure_alpha);
                    
                end
                
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
        % Export the window vector coordinates in Allen CCF to the workspace
        window_vector = cell2mat(get(gui_data.handles.window_vector_line,{'XData','YData','ZData'})');
        window_vector_ccf = round(window_vector([1,3,2],:))';
        assignin('base','window_vector_ccf',window_vector_ccf)
        disp('Copied window vector coordinates to workspace');
        
end

% Upload gui_data
guidata(window_atlas_gui, gui_data);

end

function key_release(window_atlas_gui,eventdata)

% Get guidata
gui_data = guidata(window_atlas_gui);

switch eventdata.Key
    case {'rightarrow','leftarrow','uparrow','downarrow'}
        % Update the window info/slice on arrow release
        update_window_coordinates(window_atlas_gui);
        update_slice(window_atlas_gui);
end

% Upload gui_data
guidata(window_atlas_gui, gui_data);

end


function update_slice(window_atlas_gui,varargin)

% Get guidata
gui_data = guidata(window_atlas_gui);

% Only update the slice if it's visible
if strcmp(gui_data.handles.slice_plot(1).Visible,'on')
    
    offsets = [0,gui_data.window_depth];
    plot_handles = {'window_areas','window_areas_bottom'};
    
    for islice = 1:length(offsets)
        
        % Get current window outline coordinate
        curr_window_circle = cell2mat(get(gui_data.handles.window_circle,{'XData','YData','ZData'})');
        window_vector = cell2mat(get(gui_data.handles.window_vector,{'XData','YData','ZData'})');
        move_vector = diff(window_vector,[],2)./ ...
            norm(diff(window_vector,[],2))*offsets(islice);
        window_circle = bsxfun(@plus,curr_window_circle,move_vector);
        window_centre = gui_data.window_centre + move_vector';
        
        
        % Get two vectors on the window plane
        vlength = size(window_circle,2);
        window_vector_1 = window_circle(:,1)' - window_circle(:,round(0.4*vlength))';
        window_vector_2 = window_circle(:,1)' - window_circle(:,round(0.7*vlength))';
        
        %get the normal vector of the plane
        normal_vector = cross(window_vector_1,window_vector_2);
        
        % Get the plane offset through the window centre
        window_top = window_centre;
        plane_offset = -(normal_vector*window_top');
        
        % Define a plane of points to index
        % (the plane grid is defined based on the which cardinal plan is most
        % orthogonal to the plotted plane. this is janky but it works)
        slice_px_space = 3;
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
        
        % Get the coordinates on the plane
        x_idx = round(plane_x);
        y_idx = round(plane_y);
        z_idx = round(plane_z);
        
        % make a mask for the window (in the horizontal plane)
        window_mask = uint8(((plane_x-window_centre(1)).^2 + (plane_y-window_centre(2)).^2 + (plane_z-window_centre(3)).^2)<(gui_data.window_diameter/2)^2);
        
        % Find plane coordinates in bounds with the volume and the window
        use_xd = x_idx > 0 & x_idx < size(gui_data.tv,1);
        use_yd = y_idx > 0 & y_idx < size(gui_data.tv,3);
        use_zd = z_idx > 0 & z_idx < size(gui_data.tv,2);
        use_idx = use_xd & use_yd & use_zd & window_mask;
        
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
        if islice==2 %1 = top, 2 = bottom
            set(gui_data.handles.slice_plot,'XData',plane_x,'YData',plane_y,'ZData',plane_z,'CData',curr_slice);
        end
        
        % Update the title with new depth
        if islice==2
            set(eval(['gui_data.handles.axes_',plot_handles{islice},'.Title']),'String',[sprintf('Bottom (%i ',10*offsets(islice)),'\mu','m)']);
        end
        
        % Update the window overview
        curr_slice_av = nan(size(use_idx));
        curr_slice_av(curr_slice_isbrain) = gui_data.av(grab_pix_idx);
        x_edges = [find(sum(window_mask,1),1,'first'),find(sum(window_mask,1),1,'last')];
        y_edges = [find(sum(window_mask,2),1,'first'),find(sum(window_mask,2),1,'last')];
        cropped_window = curr_slice_av(y_edges(1):y_edges(2),x_edges(1):x_edges(2))';
        set(eval(['gui_data.handles.',plot_handles{islice},'_plot']),'CData',cropped_window);
        
        %update area labels based on which areas are contained in window
        window_areas = cropped_window;
        unique_areas = unique(window_areas(~isnan(window_areas(:))));
        unique_parent_struct = unique(gui_data.st.parent_structure_id(unique_areas));
        for i_struct = 1:100
            if i_struct<=length(unique_parent_struct)
                areas_in_parent_struct = find(gui_data.st.parent_structure_id==unique_parent_struct(i_struct));
                area_bm = ismember(window_areas,areas_in_parent_struct);
                area_bm_largest = zeros(size(area_bm));
                cc = bwconncomp(area_bm);
                [~,largest_blob] = max(cellfun('length',cc.PixelIdxList));
                area_bm_largest(cc.PixelIdxList{largest_blob}) = 1;
                [area_x,area_y] = find(area_bm_largest==1);
                area_centre = [mean(area_y),mean(area_x)];
                area_string = gui_data.st.acronym{find(gui_data.st.id==unique_parent_struct(i_struct))};
            else
                area_centre = [0,0];
                area_string = '';
            end
            set(eval(['gui_data.handles.',plot_handles{islice},'_text(',num2str(i_struct),')']),...
                'Position',area_centre,'String',area_string,...
                'VerticalAlignment','middle','HorizontalAlignment','center');
        end
    end

    % Upload gui_data
    guidata(window_atlas_gui, gui_data);
    
end

end


function set_window_position(window_atlas_gui,varargin) %this function should place window at tangent plane

% Get guidata
gui_data = guidata(window_atlas_gui);

prompt_text = { ...
    'AP position (mm from bregma)', ...
    'ML position (mm from bregma)', ...
    'window diameter (mm)', ...
    'maximum depth (microns)'};

new_window_position = cellfun(@str2num,inputdlg(prompt_text,'Set window position',1));

% Convert centre position: mm->CCF
window_centre_ccf_coords = round(gui_data.bregma([1,3])' - 100*[new_window_position(1);-new_window_position(2)]);%not sure why we need this sign change...
window_diameter = new_window_position(3)*100;
window_depth = new_window_position(4)/10;

%load surface data or generate it from scratch
if ~exist('surf_coords','var')
    try
        load('dorsalsurface_nomidline_ccf_coords.mat')
    catch
        %CODE TO GENERATE SURFACE MAP
        surf_coords = nan(size(gui_data.av,1),size(gui_data.av,3));
        for i=1:size(gui_data.av,1)
            for j=1:size(gui_data.av,3)
                idxs = find(squeeze(gui_data.av(i,:,j))>1,1,'first');
                if ~isempty(idxs)
                    surf_coords(i,j) = idxs;
                end
            end
            idx_c = round(size(gui_data.av,3)/2);
            [m,idx_m] = min(surf_coords(i,:));
            idx_m_c = idx_m-idx_c;
            surf_coords(i,idx_c-abs(idx_m_c):idx_c+abs(idx_m_c)) = m;
        end
        save('dorsalsurface_nomidline_ccf_coords.mat',surf_coords);
    end
end

%figure out angle of tangent plane
surf_depth_coords = surf_ccf_coords(window_centre_ccf_coords(1),window_centre_ccf_coords(2));
centre_surf_ccf_coords = [window_centre_ccf_coords(1); surf_depth_coords; window_centre_ccf_coords(2)];
centre_patch = surf_ccf_coords(window_centre_ccf_coords(1)-50:window_centre_ccf_coords(1)+50,window_centre_ccf_coords(2)-50:window_centre_ccf_coords(2)+50);
centre_patch_smooth = imgaussfilt(centre_patch,'FilterSize',9);
[fx,fy] = gradient(centre_patch_smooth,1);
avg_gradient = [-nanmean(fy(:));nanmean(fx(:))];
window_angle_rad = atan(avg_gradient);
window_angle_deg = (window_angle_rad/(2*pi))*360;

%compute new window vector and circle (horizontal)
%window_diameter = gui_data.window_diameter;
window_rads = linspace(-pi,pi,100);
window_centre = [window_centre_ccf_coords', surf_depth_coords-gui_data.ref_centre(3)];
window_circle = [0.5*window_diameter*cos(window_rads)+window_centre(1);...
    0.5*window_diameter*sin(window_rads)+window_centre(2);...
    zeros(size(window_rads))+window_centre(3)];
window_vector = [window_centre',(window_centre+[0,0,window_depth])'];

%rotate these by the tangent angle
eul = [0;window_angle_rad]; %[yaw; pitch; roll]
R = eul2rotm(eul'); %find rotation matrix
window_circle_centred = window_circle - window_centre'; %translate to origin by subtracting its centre coords
window_vector_centred = window_vector - window_centre';
new_window_circle_centred = R*window_circle_centred; %rotate around origin
new_window_vector_centred = R*window_vector_centred;
new_window_circle = new_window_circle_centred + window_centre'; %translate back to its position
new_window_vector = new_window_vector_centred + window_centre';

% update gui_data
set(gui_data.handles.window_circle,'XData',new_window_circle(1,:), ...
    'YData',new_window_circle(2,:), ...
    'ZData',new_window_circle(3,:));
set(gui_data.handles.window_vector,'XData',new_window_vector(1,:), ...
    'YData',new_window_vector(2,:), ...
    'ZData',new_window_vector(3,:));
gui_data.window_angle = window_angle_deg;
gui_data.window_centre = [window_centre_ccf_coords', surf_depth_coords];
gui_data.window_diameter = window_diameter;
gui_data.window_depth = window_depth;

% Upload gui_data
guidata(window_atlas_gui, gui_data);

% Update the slice and window coordinates
update_slice(window_atlas_gui);
update_window_coordinates(window_atlas_gui);

end

function gui_data = update_window_angle(window_atlas_gui,angle_change)

% Get guidata
gui_data = guidata(window_atlas_gui);

% Set new angle
new_angle = gui_data.window_angle + angle_change;
gui_data.window_angle = new_angle;

% Get the positions of the window circle and vector
window_circle = cell2mat(get(gui_data.handles.window_circle,{'XData','YData','ZData'})');
window_vector = cell2mat(get(gui_data.handles.window_vector,{'XData','YData','ZData'})');

% Compute the rotation matrix
eul = [0;(angle_change/360)*2*pi];
R = eul2rotm(eul');

% temporarily translate window position to origin
window_circle_centred = window_circle - gui_data.window_centre';
window_vector_centred = window_vector - gui_data.window_centre';

% rotate the vector around the origin using rotation matrix
new_window_circle_centred = R*window_circle_centred;
new_window_vector_centred = R*window_vector_centred;

% translate rotated window back to the original centre point
new_window_circle = new_window_circle_centred + gui_data.window_centre';
new_window_vector = new_window_vector_centred + gui_data.window_centre';

% update gui_data
set(gui_data.handles.window_circle,'XData',new_window_circle(1,:), ...
    'YData',new_window_circle(2,:), ...
    'ZData',new_window_circle(3,:));
set(gui_data.handles.window_vector,'XData',new_window_vector(1,:), ...
    'YData',new_window_vector(2,:), ...
    'ZData',new_window_vector(3,:));

% Upload gui_data
guidata(window_atlas_gui, gui_data);

end

function update_window_coordinates(window_atlas_gui,varargin)

% Get guidata
gui_data = guidata(window_atlas_gui);

window_bregma_coordinate = round((gui_data.bregma([1,3,2]) - gui_data.window_centre)*10);

% Update the text
window_text = ['Window position: ' ....
    num2str(window_bregma_coordinate(1)) ' AP, ', ...
    num2str(-window_bregma_coordinate(2)) ' ML, ', ...
    num2str(-window_bregma_coordinate(3)) ' DV, ', ...
    num2str(round(gui_data.window_angle(1))) char(176) ' pitch, ' ...
    num2str(round(gui_data.window_angle(2))) char(176) ' roll'];
set(gui_data.window_coordinates_text,'String',window_text);

% Upload gui_data
guidata(window_atlas_gui, gui_data);

end


function display_controls

% Print controls
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
msgbox( ...
    {'\fontsize{12}' ...
    '\bf Window: \rm' ...
    'Arrow keys : translate window' ...
    'Alt/Option up/down : raise/lower window' ...
    'Shift arrow keys : change window angle' ...
    'm : set location manually (tangent to surface)', ...
    '\bf 3D brain areas: \rm' ...
    ' =/+ : add (list selector)' ...
    ' Alt/Option =/+ : add (search)' ...
    ' Shift =/+ : add (hierarchy selector)' ...
    ' Ctrl Shift =/+ : add all regions on surface' ...
    ' - : remove', ...
    '\bf Visibility: \rm' ...
    's : atlas slice (toggle tv/av/off)' ...
    'b : brain outline' ...
    'w : window outline' ...
    'a : 3D brain areas (toggle transparency)' ...
    '\bf Other: \rm' ...
    'x : export window coordinates to workspace' ...
    'c : bring up controls box'}, ...
    'Controls',CreateStruct);

end







