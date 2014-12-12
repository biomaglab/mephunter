% function to eliminate single electrodes far away from clusters
function mask_fin = clusterfilt(mask_ini)
mask_fin = mask_ini;
% first neighbors
[cord(:,1), cord(:,2)] = find(mask_ini);
mm = ones(3,3);
mm(2,2) = 0;
if size(cord, 1) > 1
    for i = 1:size(cord, 1)
        % first row
        if cord(i,1) == 1
            if cord(i,2) == 1
                aux  = mm(2:3, 2:3).*mask_ini(cord(i,1):cord(i,1)+1, cord(i,2):cord(i,2)+1);
            elseif cord(i,2) == 5
                aux  = mm(2:3, 1:2).*mask_ini(cord(i,1):cord(i,1)+1, cord(i,2)-1:cord(i,2));
            else
                aux  = mm(2:3, 1:3).*mask_ini(cord(i,1):cord(i,1)+1, cord(i,2)-1:cord(i,2)+1);
            end
            % last row
        elseif cord(i,1) == 13
            if cord(i,2) == 1
                aux  = mm(1:2, 2:3).*mask_ini(cord(i,1)-1:cord(i,1), cord(i,2):cord(i,2)+1);
            elseif cord(i,2) == 5
                aux  = mm(1:2, 1:2).*mask_ini(cord(i,1)-1:cord(i,1), cord(i,2)-1:cord(i,2));
            else
                aux  = mm(1:2, 1:3).*mask_ini(cord(i,1)-1:cord(i,1), cord(i,2)-1:cord(i,2)+1);
            end
            % middle rows
        else
            if cord(i,2) == 1
                aux  = mm(1:3, 2:3).*mask_ini(cord(i,1)-1:cord(i,1)+1, cord(i,2):cord(i,2)+1);
            elseif cord(i,2) == 5
                aux  = mm(1:3, 1:2).*mask_ini(cord(i,1)-1:cord(i,1)+1, cord(i,2)-1:cord(i,2));
            else
                aux  = mm.*mask_ini(cord(i,1)-1:cord(i,1)+1, cord(i,2)-1:cord(i,2)+1);
            end
        end
        
        if cumsum(cumsum(aux,1),2) == 0;
            mask_fin(cord(i,1), cord(i,2)) = 0;
        end
    end
end
end