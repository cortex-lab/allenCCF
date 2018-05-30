
function HistologyHotkeyFcn(fig, keydata)

ud = get(fig, 'UserData');

switch lower(keydata.Key)    
    case 'd' % lower min contrast
        if ud.contrast(1) >= .05
            ud.contrast(1) = ud.contrast(1) - .05; ud.contrast
        end
    case 'a' % raise min contrast
        if ud.contrast(1) <= .95
            ud.contrast(1) = ud.contrast(1) + .05; ud.contrast
        end
    case 's' % raise max contrast
        if ud.contrast(2) <= .95
            ud.contrast(2) = ud.contrast(2) + .05; ud.contrast
        end
    case 'w' % lower max contrast
        if ud.contrast(2) >= .05
            ud.contrast(2) = ud.contrast(2) - .05; ud.contrast
        end   
    case 'e' % show original
        ud.contrast
    case 'r' % return to original
        ud.contrast = [0 1];       
    case 'space'
        disp('switch contrast effect')
        if ud.contrast_type==2
            ud.contrast_type = 1;
        elseif ud.contrast_type==1
            ud.contrast_type = 2;
    end
    case 'q' % break
        disp('done with this channel')
        ud.break = 1;    
end

ud.key = keydata.Key;
set(fig, 'UserData', ud);