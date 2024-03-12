function y_mm = accf2pxs_mm(x_mm, plane, numtype)
% convert Allen Common Coordinates into standard values (Paxinos and Franklin, 2012, 4th edition)
%
% SYNTAX
% y_mm = accf2pxs_mm(x_mm, plane)
%
% INPUT ARGUMENTS
% x_mm        Allen Common Coordinates
%             if this is scaler, x_mm should be depth.
%             if this is n x 3 array, x_mm is an array of the same (or similar) size
%
% plane       'coronal' | 'sagittal' | 'transverse'
%
% numtype     'coordinate' (defalut) |  'distance'
%
% OUTPUT ARGUMENTS
% y_mm        Paxinos and Franklin coordinates
%
% Written by Kouichi C. Nakamura Ph.D.
% MRC Brain Network Dynamics Unit
% University of Oxford
% kouichi.c.nakamura@gmail.com
% 09-Aug-2023 11:06:59
%
% See also
% apdvml2info

arguments
 x_mm 
 plane (1,1) string {mustBeMember(plane, {'coronal','sagittal','transverse'})}
 numtype (1,1) string {mustBeMember(numtype, {'coordinate','distance'})} = "coordinate"

end


switch plane

    case "coronal"
        if size(x_mm,2) ==1
            switch numtype
                case "coordinate"
                    y_mm = 0.825 * x_mm - 0.131;
                case "distance"
                    y_mm = 0.825 * x_mm ; 
            end

        elseif size(x_mm,2) == 3

            % assume that the 2nd column is the DV
            y_mm = x_mm;
            switch numtype
                case "coordinate"
                    y_mm(:,2) = 0.825 * x_mm(:,2) - 0.131;
                case "distance"
                    y_mm = 0.825 * x_mm(:,2) ;
            end
        else
            error('unexpected size')
        end        
    case "sagittal"
        if size(x_mm,2) ==1
            switch numtype
                case "coordinate"
                    y_mm = 0.921 * x_mm - 0.834;
                case "distance"
                    y_mm = 0.921 * x_mm ;
            end
        elseif size(x_mm,2) == 3
            % assume that the 2nd column is the DV
            y_mm = x_mm;
            switch numtype
                case "coordinate"
                    y_mm(:,2) = 0.921 * x_mm(:,2) - 0.834;
                case "distance"
                    y_mm = 0.921 * x_mm(:,2) ;
            end
        else
            error('unexpected size')
        end
    case "transverse"
        error("not implemented yet")
    otherwise
        error("wrong input for 'plane'")
end