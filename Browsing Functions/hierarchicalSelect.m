

function selIdx = hierarchicalSelect(st)

selID = 567; % Cerebrum, default to start

[boxList, idList] = makeBoxList(st, selID); 

ud.idList = idList; ud.st = st;

% make figure
f = figure; set(f, 'KeyPressFcn', @doOK);

% selector box
ud.hBox = uicontrol(f, 'Style', 'listbox', 'String', boxList, ...
    'Callback', @updateSel, 'Value', find(idList==selID),...
    'Units', 'normalized', 'Position', [0.1 0.2 0.8 0.7],...
    'KeyPressFcn', @doOK); 

titleStr = boxList{idList==selID}; titleStr = titleStr(find(titleStr~=' ', 1):end);
ud.hSelTitle = uicontrol(f, 'Style', 'text', ...
    'String', sprintf('Selected: %s', titleStr), ...
    'Units', 'normalized', 'Position', [0.1 0.9 0.8 0.1]); 

ud.hCancel = uicontrol(f, 'Style', 'pushbutton', ...
    'String', 'Cancel', 'Callback', @doCancel, ...
    'Units', 'normalized', 'Position', [0.1 0.1 0.2 0.1]); 

ud.hOK = uicontrol(f, 'Style', 'pushbutton', ...
    'String', 'OK', 'Callback', @doOK, ...
    'Units', 'normalized', 'Position', [0.3 0.1 0.2 0.1]); 



set(f, 'UserData', ud);
drawnow;

try
    c = matlab.ui.internal.dialog.DialogUtils.disableAllWindowsSafely();
    uiwait(f);
    delete(c);
catch
    if ishghandle(f)
        delete(f)
    end
end

if ishghandle(f)
    ud = get(f, 'UserData');
    idList = ud.idList;

    if ud.hBox.Value>1
        selID = idList(get(ud.hBox, 'Value'));

        selIdx = find(st.id==selID);
    else
        selIdx = [];
    end
    delete(f)
    drawnow; 
else
    selIdx = [];
end



end

function [boxList, idList] = makeBoxList(st, selID)

idList = selID;
while idList(end)~=997 % root
    idList(end+1) = st.parent_structure_id(st.id==idList(end));
end
idList = idList(end:-1:1); % this is the tree of parents down to the selected one

% make the parent string representation
for q = 1:numel(idList)
    boxList{q} = sprintf('%s%s (%s)', repmat('  ', 1, q-1), ...
        st.acronym{st.id==idList(q)}, ...
        st.safe_name{st.id==idList(q)});
end
np = numel(idList);

% now add children
idList = [idList st.id(st.parent_structure_id==selID)'];

% make the parent string representation
for q = np+1:numel(idList)
    boxList{q} = sprintf('%s%s (%s)', repmat('  ', 1, np), ...
        st.acronym{st.id==idList(q)}, ...
        st.safe_name{st.id==idList(q)});
end

end

function updateSel(src, ~)

f = get(src, 'Parent'); 
ud = get(f, 'UserData');
st = ud.st; idList = ud.idList;

selID = idList(get(src, 'Value'));

[boxList, idList] = makeBoxList(st, selID); 

ud.idList = idList;
set(f, 'UserData', ud);
set(src, 'String', boxList, 'Value', find(idList==selID));

titleStr = boxList{idList==selID}; titleStr = titleStr(find(titleStr~=' ', 1):end);
set(ud.hSelTitle, 'String', sprintf('Selected: %s', titleStr));

end



% OK callback
function doOK(~, ~)
    uiresume(gcbf);
end

% Cancel callback
function doCancel(~, ~)
ud = get(gcbf, 'UserData');
ud.hBox.Value = 1;
uiresume(gcbf);
end