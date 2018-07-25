

function plotAsProbe(data, ann_percents, yc, cm, cm2, sqSizeY, active_site_start, rpl)
% data is vector, same lenght as xc and yc
% cm is colormap, size [n x 3]
%
% other usage: if data is empty and cm's first dimension is the same length
% as xc and yc, then use those colors literally
sqCoordsY = [sqSizeY/2 sqSizeY/2 -sqSizeY/2 -sqSizeY/2];

for q = 1:size(ann_percents,1)
    
    % get x-boundaries corresponding to percent certainty of region
    percent_region = [0 ann_percents(q,1)*100 ann_percents(q,1)*100 0];
    percent_second_region = [100-ann_percents(q,2)*100 100 ...
                                        100 100-ann_percents(q,2)*100];
    y = yc(q)+sqCoordsY;        

    if yc(q) < active_site_start  || yc(q) > rpl*10
        cur_alpha = .4;
    else
        cur_alpha = .8;
    end
    
    fill(percent_region,y,cm(q,:), 'EdgeAlpha', 0, 'FaceAlpha', cur_alpha);
    fill(percent_second_region,y,cm2(q,:), 'EdgeAlpha', 0, 'FaceAlpha', cur_alpha);
    
    if ann_percents(q,1) + ann_percents(q,2) < .995
        percent_in_between = [ann_percents(q,1)*100 100-ann_percents(q,2)*100 ...
                                        100-ann_percents(q,2)*100 ann_percents(q,1)*100];
        fill(percent_in_between,y,[.8 .8 .8], 'EdgeAlpha', 0, 'FaceAlpha', cur_alpha);
    end
    
    hold on;
end


