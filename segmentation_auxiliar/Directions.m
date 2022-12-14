function [principalDirections] = Directions(vxLocation,idx, epsilon)
%% Directions
% Apply PCA to analise the direction of each cluster using their 1st
% eigenvector.
% dbscan to detect the 2 principal directions
% 
%--------------------------------------------------------------------------
% INPUTS:
% vxLocation: numeric Nx3. Location of a cloud
%
% idx: numeric Nx1. Index of the cluster of each point.
%
%--------------------------------------------------------------------------
% OUTPUTS:
% principalDirections: numeric Mx3. Principal directions of the analysed
%                      clusters sorted by its Z.
%
%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 21/04/21

%%
[num] = groupcounts(idx); 
[num, clusters] = sort(num,'descend');

directions = zeros(length(clusters),3);
for i = 1:numel(clusters)
    bar = vxLocation(idx == clusters(i),:);
    direction = pcacov(cov(bar));
    if all(size(direction) == [3,3])
        if direction(1,1) < 0
            rot180 = makehgtform('axisrotate',direction(:,3),pi);
            rot180 = rot180(1:3,1:3); % 3x3 matrix
            direction = rot180 * direction;    
        end
        directions(i,:) = direction(:,1)';
    end
end

% Remove directions = 0,0,0
num = num(~all(directions == [0,0,0],2),:);
directions = directions(~all(directions == [0,0,0],2),:);
% figure; pcshow(directions,'g', 'MarkerSize', 50);

% Grouping clusters with the same direction.
idx2 = dbscan(directions, epsilon,10); % antes 20. Mirar esto con el de Dom Luis
clusters_directions = unique(idx2);
clusters_directions = clusters_directions(clusters_directions ~= -1);

% Sum the number of points in clusters with the same directions.
num_direction_cluster = zeros(size(clusters_directions));
for i = 1:length(clusters_directions)
    num_direction_cluster(i) = sum(num(idx2 == clusters_directions(i)));
end

% Sorting clusters_directions by the total number of points of their clusters.
[~, order] = sort(num_direction_cluster,'descend');
clusters_directions = clusters_directions(order);

% i = 0;
% i = i+1;
% hold on; pcshow(directions(idx2 == clusters_directions(i),:),'y', 'MarkerSize', 50);

% Picking the 2 principal directions or NaN if there are not at least 2.
principalDirections = zeros(2,3);
if length(clusters_directions) < 2
    principalDirections = NaN;
    return
end
for i = 1:2
    principalDirections(i,:) = mean(directions(idx2 == clusters_directions(i),:));
end

%% Normalising
for i = 1: length(principalDirections(:,1))
    principalDirections(i,:) = principalDirections(i,:) / norm(principalDirections(i,:));
end

%% Sorting principalDirections using their Z
[~,order] = sort(principalDirections(:,3), 'descend');
principalDirections = principalDirections(order,:);

end

