function y_mm = accf2pxs_mm(x_mm)
% convert Allen Common Coordinates into standard values (Paxinos and Franklin)
% x_mm   Allen Common Coordinates 
% y_mm Paxinos and Franklin coordinates

if size(x_mm,2) ==1 
    y_mm = 0.921 * x_mm - 0.834; % global
    % y = 0.908 * x - 0.957; % local

elseif size(x_mm,2) == 3
    % assume that the 2nd column is the DV
    y_mm = x_mm;
    y_mm(:,2) = 0.921 * x_mm(:,2) - 0.834; % global
    % y(:,2) = 0.908 * x(:,2) - 0.957;  % local

else
    error('unexpected size')
end