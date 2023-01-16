% Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
% Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
% Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
% Under the terms of the GNU General Public License.
% See the LICENCE.md and COPYING.md file for more details.

function [idxOut]= UnirPalitroques(vxLocation, idxIn, width, directions, tolerance, minSeed)
%% UnirPalitroques
% Prolong the groups following a main direction and add to its
% group any other that falls into the prolongation and shares
% its direction.
%
% Generate AXIs for each direction. Necesary to analyse the clusters since
% the vector direction is not enough. The 1st vector is the direction, the
% 2nd is in the XY plane and the 3rd is orthonormal.
% Each cluster is assinged to a principal direction if its own direction 
% is similar to any principal direction.
% Starting with the bigger clusters, the clusters are merged. All the
% cluster are reorineted using the AXIs of the cluster under study. The
% groups with any voxel inside the prologation of the cluster under study,
% assigned to the same principal direction, and that are not for any group
% are merged. If some groups are marge, this process is repeated using the
% new group to compare with the others.
%
%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 21/04/21

%% AXIs for principal direction
barPca = cell(length(directions(:,1)),1);
for i = 1:numel(barPca)
    barPca{i} = zeros(3);
    barPca{i}(:,1) = directions(i,:);
    if all(barPca{i}(1:2,1) == [0;0])
        barPca{i}(:,2) = [-barPca{i}(3,1), 0, barPca{i}(1,1)] / norm([-barPca{i}(3,1), 0, barPca{i}(1,1)]); % orthonormal to pca(1) and in the xy plane
    else
        barPca{i}(:,2) = [-barPca{i}(2,1), barPca{i}(1,1), 0] / norm([-barPca{i}(2,1), barPca{i}(1,1), 0]); % orthonormal to pca(1) and in the xy plane
    end
    
    barPca{i}(:,3) = cross(barPca{i}(:,1), barPca{i}(:,2)); % orthonormal to mpca(1) and mpca(2)
end

%% Analysis of the clusters
[num] = groupcounts(idxIn); 
[num, clusters] = sort(num,'descend');

% clusters = clusters(num > minSeed); % removing the smallest seeds
clusters = clusters(num > mean(num)); % removing the smallest seeds
% Analysis of their direction
seedsDirection = zeros(length(clusters(:,1)),1);
for i = 1:numel(clusters)
    direction = pca( vxLocation(idxIn == clusters(i),:));
    if all(size(direction) == [3,3])
        for j = 1:length(barPca)
            if abs(dot(direction(:,1)',barPca{j}(:,1))) > tolerance
                seedsDirection(i) = j;
                break;
            end
        end
    end
end

% figure; pcshow(vxLocation, 'g', 'MarkerSize', 50);
% hold on; pcshow(vxLocation(ismember(idxIn, clusters),:), 'r', 'MarkerSize', 50);
% 
% i = 0;
% i = i+1
% hold on; pcshow(vxLocation(ismember(idxIn, clusters(i)),:), 'y', 'MarkerSize', 50);


%% Merging clusters
idxSeed = zeros(size(clusters)); % Detec
idxOut = zeros(length(idxIn),1);
i = 1;
while i <= length(clusters(:,1))
    clustersInRange = [];
    if seedsDirection(i) ~= 0 && any(idxIn == clusters(i)) % cluster with a correct direction and not assigned to any bar
        bar = vxLocation(idxIn == clusters(i),:); % Orienting the cloud
        vxOriented = vxLocation - mean(bar);
        vxOriented = vxOriented * barPca{seedsDirection(i)};
        
        vxInRange = sqrt(vxOriented(:,2).^2 + vxOriented(:,3).^2) < width; % voxels in range
        clustersInRange = unique(idxIn(vxInRange)); % clusters with any voxel in range
        clustersInRange = clusters(seedsDirection == seedsDirection(i) & idxSeed == 0 & ismember(clusters,clustersInRange)); % clusters in range with the same direction than cluster(i)
        idxSeed(ismember(clusters,clustersInRange)) = i; 
        
        % Change idxIn of the clustersInRange and wirthe in idxOut the
        % orientation
        if ~isempty(clustersInRange)
            idxIn(ismember(idxIn, clustersInRange)) = clusters(i);
            idxOut(ismember(idxIn, clusters(i))) = seedsDirection(i);
        end
    end
    
    if isempty(clustersInRange)
        i = i + 1;
    end
end

% figure; pcshow(vxOriented, 'w');
% hold on; pcshow(vxOriented(ismember(idxIn, clustersInRange), :), 'g');
% hold on; pcshow(vxOriented(idxIn == clusters(i), :), 'r');

idxOut = [idxIn, idxOut];

end

% figure; pcshow(vxOriented, 'g', 'MarkerSize', 50);
% hold on; pcshow(vxOriented(totalBar,:), 'y', 'MarkerSize', 50);
% hold on; pcshow(vxOriented(idx == clusters(i),:), 'r', 'MarkerSize', 50);
% 
% figure; pcshow(vxLocation, 'w', 'MarkerSize', 50);
% hold on; pcshow(vxLocation(idxOut(:,2) == 1,:), 'g', 'MarkerSize', 50);
% hold on; pcshow(vxLocation(idxOut(:,2) == 2,:), 'r', 'MarkerSize', 50);
