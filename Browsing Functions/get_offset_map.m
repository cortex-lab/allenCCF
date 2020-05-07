% ------------------------------------------------------------
% Generate offset map (for third dimension of a tilted slice)
% ------------------------------------------------------------

function offset_map = get_offset_map(slice_angle, ref_size)

  offset_map = zeros(ref_size);

  % get range of offsets along 1st angle dimension to query for offset in the third dimension
  if slice_angle(1)==0; offset_AP = 0;
  else; offset_AP = -slice_angle(1):sign(slice_angle(1)):slice_angle(1);
  end; start_index_AP = 1; 
 
  % loop through AP offsets
  for cur_offset_AP = offset_AP
    if cur_offset_AP == slice_angle(1); end_index_AP = 800;
    else; end_index_AP = start_index_AP + floor( 800 / length(offset_AP)) - 1;
    end
      
    % get range of offsets along 2nd angle dimension to query for offset in the third dimension
    if slice_angle(2)==0;  offset_ML = 0;
    else; offset_ML = -slice_angle(2):sign(slice_angle(2)):slice_angle(2);
    end; start_index_ML = 1;
    
    % nested: loop through ML offsets
    for cur_offset_ML = offset_ML
      if cur_offset_ML == slice_angle(2)
         end_index_ML = 1140;
      else
         end_index_ML = start_index_ML + floor( 1140 / length(offset_ML)) - 1;
      end

      % generate offset map
      offset_map(start_index_AP:end_index_AP, start_index_ML:end_index_ML) = cur_offset_AP + cur_offset_ML;

      start_index_ML = end_index_ML + 1;
    end
      start_index_AP = end_index_AP + 1;
  end  
  
end
