%% Cluster test
clear all;
load('maps_rms.mat');
% coordinates of electrodes
[x, y] = meshgrid(1:5, 1:13);
% amplitudes for electrodes on cluster from watershed
d = masks.*data;
cls1 = zeros(13,5,8);
% cluster calculation
for i = 1:8
    d1 = d(:,:,i);
    % reshaped amplitude and coordinates map
    d1r = reshape(d1,65,1);
    xr = reshape(x,65,1);
    yr = reshape(y,65,1);
    % data set prepared for hierarchical clustering (X, Y, AMP)
    dataset = horzcat(xr, yr, d1r);
%     cls1r = clusterdata(dataset, 'maxclust', 2);
    cls1r = clusterdata(dataset, 'criterion', 'distance',...
        'distance', 'spearman', 'maxclust', 3);
    % find max peak-peak amplitude
    ind = find(d1r == max(d1r));
    clmax = (cls1r == cls1r(ind)).*cls1r;
    cls1(:,:,i) = reshape(clmax,13,5);
    % plot of clusters
    figure; imagesc(cls1(:,:,i));
    figure; imagesc(d1);
end

% m2 = cls1;
% % only execute cluster cleaning if the cluster has more than one electrode
% for j = 1:8
%     m1 = masks(:,:,j);
%     % first neighbors
%     [cord(:,1), cord(:,2)] = find(m1);
%     mm = ones(3,3);
%     mm(2,2) = 0;
%     if size(cord, 1) > 1
%         for i = 1:size(cord, 1)
%             % first row
%             if cord(i,1) == 1
%                 if cord(i,2) == 1
%                     aux  = mm(2:3, 2:3).*m1(cord(i,1):cord(i,1)+1, cord(i,2):cord(i,2)+1);
%                 elseif cord(i,2) == 5
%                     aux  = mm(2:3, 1:2).*m1(cord(i,1):cord(i,1)+1, cord(i,2)-1:cord(i,2));
%                 else
%                     aux  = mm(2:3, 1:3).*m1(cord(i,1):cord(i,1)+1, cord(i,2)-1:cord(i,2)+1);
%                 end
%                 % last row
%             elseif cord(i,1) == 13
%                 if cord(i,2) == 1
%                     aux  = mm(1:2, 2:3).*m1(cord(i,1)-1:cord(i,1), cord(i,2):cord(i,2)+1);
%                 elseif cord(i,2) == 5
%                     aux  = mm(1:2, 1:2).*m1(cord(i,1)-1:cord(i,1), cord(i,2)-1:cord(i,2));
%                 else
%                     aux  = mm(1:2, 1:3).*m1(cord(i,1)-1:cord(i,1), cord(i,2)-1:cord(i,2)+1);
%                 end
%                 % middle rows
%             else
%                 if cord(i,2) == 1
%                     aux  = mm(1:3, 2:3).*m1(cord(i,1)-1:cord(i,1)+1, cord(i,2):cord(i,2)+1);
%                 elseif cord(i,2) == 5
%                     aux  = mm(1:3, 1:2).*m1(cord(i,1)-1:cord(i,1)+1, cord(i,2)-1:cord(i,2));
%                 else
%                     aux  = mm.*m1(cord(i,1)-1:cord(i,1)+1, cord(i,2)-1:cord(i,2)+1);
%                 end
%             end
%             
%             if cumsum(cumsum(aux,1),2) == 0;
%                 m2(cord(i,1), cord(i,2), j) = 0;
%             end
%         end
%     end
%     clear cord
%     figure; imagesc(m2(:,:,j));
% end

m2 = cls1;
% only execute cluster cleaning if the cluster has more than one electrode
for j = 1:8
    m2(:,:,j) = clusterfilt(masks(:,:,j));
    figure; imagesc(m2(:,:,j));
end