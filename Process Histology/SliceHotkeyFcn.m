
function SliceHotkeyFcn(fig, keydata)

ud = get(fig, 'UserData');


if strcmp(keydata.Key,'q') || strcmp(keydata.Key,'space') % break
    ud.break = 1;    
elseif strcmp(keydata.Key,'g') %|| strcmp(keydata.Key,'c')
    if (sum(ud.grid(:))==0 && strcmp(keydata.Key,'g')) %|| (sum(ud.grid(:))~=0 && strcmp(keydata.Key,'c'))
        ud.grid = zeros(ud.size,'uint8');
        ud.grid(1:50:end,:,:) = 150; ud.grid(:,1:50:end,:) = 150;
    else
        ud.grid = zeros(ud.size,'uint8');
    end
    
end

ud.key = keydata.Key;
set(fig, 'UserData', ud);