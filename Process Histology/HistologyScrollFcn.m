function HistologyScrollFcn(fig, evt)

ud = get(fig, 'UserData');
ud.key = 'scroll';


%modify based on scrolling
ud.contrast(ud.contrast_type) = ud.contrast(ud.contrast_type) + evt.VerticalScrollCount*.05;

% make sure within limit of 0 to 1
if ud.contrast(ud.contrast_type) < 0
    ud.contrast(ud.contrast_type) = 0;
elseif ud.contrast(ud.contrast_type) > 1
    ud.contrast(ud.contrast_type) = 1;
end


set(fig, 'UserData', ud);



