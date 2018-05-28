function SliceScrollFcn(fig, evt)

ud = get(fig, 'UserData');
ud.key = 'scroll';


%modify based on scrolling
ud.rotate_angle = ud.rotate_angle + evt.VerticalScrollCount*.333;



set(fig, 'UserData', ud);



