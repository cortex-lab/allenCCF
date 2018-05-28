
function FlipHotkeyFcn(fig, keydata)

ud = get(fig, 'UserData');


if strcmp(keydata.Key,'q') || strcmp(keydata.Key,'leftarrow')  || strcmp(keydata.Key,'rightarrow') % break
    ud.break = 1;    
elseif strcmp(keydata.Key,'f')  % flip
    if ud.flip == 1
        ud.flip = 0;
    else
        ud.flip = 1;  
    end
elseif strcmp(keydata.Key,'s')  % move forward    
    disp('')
end



ud.key = keydata.Key;
set(fig, 'UserData', ud);