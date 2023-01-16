% Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
% Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
% Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
% Under the terms of the GNU General Public License.
% See the LICENCE.md and COPYING.md file for more details.

function [lateralDiagonals] = LateralDiagonals(varargin)
%% Lateral Diagonals:
% Select voxels using localPcaThreshold, specifing the minimum and maximun
% eigenvalues and eigenvectors.
% localPcaThreshold{1}{1} and localPcaThreshold{1}{1} contains the minimum
% and maximum eigenvectors, size(3,3) each. 
% localPcaThreshold{2}{1} and localPcaThreshold{2}{2} contains the minimum
% and maximum eigenvalues, size(1,3) each. 
% Group neighbouring voxels.
% Calculate the 2 predominant directions of the groups if directions are
% not specified.
% Prolong the groups following a main direction and add to its group any 
% other that falls into the prolongation and shares its direction.
% Erase small clusters.
% Erase bars that are outside of the face plane.
% Denoise bars.

%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 29/07/21

%% Checking inputs
parser = inputParser;
parser.CaseSensitive = true;

parser.addRequired('vx', @(x)validateattributes(x,{'Voxels'}, {}));
parser.addRequired('vxLocalPca', @(x)validateattributes(x,{'struct'}, {})); % Reslut of localPca();
parser.addRequired('minPointsSeed', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]}));
parser.addRequired('toleranceDirection', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]}));
parser.addRequired('barWidth', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]}));
parser.addRequired('minPointsBar', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]}));
parser.addRequired('localPcaThreshold', @(x)validateattributes(x,{'cell'}, {}));
parser.addOptional('directions', [], @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [2,3]}));
parser.addOptional('clean_y', true, @(x)validateattributes(x,{'logical'},{}));

parser.parse(varargin{:});

vx = parser.Results.vx;
vxLocalPca = parser.Results.vxLocalPca;         
minPointsSeed = parser.Results.minPointsSeed;  
toleranceDirection = parser.Results.toleranceDirection;
barWidth = parser.Results.barWidth;         
minPointsBar = parser.Results.minPointsBar; 
localPcaThreshold = parser.Results.localPcaThreshold;
directions = parser.Results.directions;
clean_y = parser.Results.clean_y;

%% Removing nodes
minEigenvector = localPcaThreshold{1}{1};
maxEigenvector = localPcaThreshold{1}{2};
minEigenvalue  = localPcaThreshold{2}{1};
maxEigenvalue  = localPcaThreshold{2}{2};

bar = find(~any(abs(vxLocalPca.eigenvectors) <= minEigenvector,2) & ~any(vxLocalPca.eigenvalues <= minEigenvalue,2) ...
                & ~any(abs(vxLocalPca.eigenvectors) >= maxEigenvector,2) & ~any(vxLocalPca.eigenvalues >= maxEigenvalue,2));

%     figure; pcshow(vx.Location,'g', 'MarkerSize', 50);
%     hold on; pcshow(vx.Location(bar,:), 'r','MarkerSize', 50);

%% Clustering             
vxBar = select(vx, bar);
idxSeed = dbscan(vxBar.Location, 1.3*vx.grid,1);

% [num] = groupcounts(idxSeed); 
% [~, clusters] = sort(num,'descend');
% figure; pcshow(vxBar.Location, 'w', 'MarkerSize', 50);
% j = 0;
% j = j +1
% hold on; pcshow(vxBar.Location(idxSeed == clusters(j),:), 'r', 'MarkerSize', 50);

%% Merging clusters by their directions
if isempty(directions) % If directions is not specified they are calculated.
    directions = Directions(vxBar.Location,idxSeed, 2*vxBar.grid);
end

if isnan(directions)
    lateralDiagonals{1} = [];
    lateralDiagonals{2} = [];
    return
end
idx = UnirPalitroques(vxBar.Location, idxSeed, barWidth*2, directions, toleranceDirection, minPointsSeed);

% figure; pcshow(vx.Location,'w', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(bar,:),'b', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(bar(idx(:,2) == 1),:), 'r','MarkerSize', 50);
% hold on; pcshow(vx.Location(bar(idx(:,2) == 2),:), 'y','MarkerSize', 50);

%% Deleting small clusters
num = groupcounts(idx(:,1)); 
[num, order] = sort(num,'descend');
clusters = unique(idx(:,1));
clusters = clusters(order);
num      = num((clusters ~=0) & clusters ~=-1);
clusters = clusters(clusters ~=0 & clusters ~=-1);
clusters = clusters(num > minPointsBar);

if isempty(clusters)
    lateralDiagonals{1} = [];
    lateralDiagonals{2} = [];
    return
end

%% Orienting the section with the selected clusters
vx.Location = vx.Location * PcaFlattering(vx.Location(bar(ismember(idx(:,1), clusters)),:));
vx.Location = vx.Location - mean(vx.Location(bar(ismember(idx(:,1), clusters)),:));

% figure; pcshow(vx.Location(bar(ismember(idx(:,1), clusters)),:), 'MarkerSize', 50);
% i = 0;
% i = i+1
% hold on; pcshow(vx.Location(bar(ismember(idx(:,1), clusters(33))),:), 'g', 'MarkerSize', 50);

%% Spliting the 2 types of bars and calculating their mean X
clusters = [clusters(:,1), zeros(length(clusters(:,1)),3)];
for j = 1:length(clusters(:,1))
    if all(idx(idx(:,1) == clusters(j),2) == 1)
        clusters(j,2) = 1;
    elseif all(idx(idx(:,1) == clusters(j),2) == 2)
        clusters(j,2) = 2;
    end

    clusters(j,3) = mean(vx.Location(bar(idx(:,1) == clusters(j)),1)); % mean X
    clusters(j,4) = abs(mean(vx.Location(bar(idx(:,1) == clusters(j)),2))); % mean Y
end

%% Removing false positives analysing their Y deviation from the centre
if clean_y
    clusters = clusters(clusters(:,4) < 1.33*barWidth,:);
end

lateralDiag1 = clusters(clusters(:,2) == 1,[1,3]);
lateralDiag2 = clusters(clusters(:,2) == 2,[1,3]);

% Sorting each type of bar in X
[~,order] = sort(lateralDiag1(:,2),'ascend');
lateralDiag1 = lateralDiag1(order,:);
[~,order] = sort(lateralDiag2(:,2),'ascend');
lateralDiag2 = lateralDiag2(order,:);

%% Saving
lateralDiagonals = cell(1,2);

for j = 1:length(lateralDiag1(:,1))
    thisBar = bar(idx(:,1) == lateralDiag1(j));
    [~,idxDenoised,~] = pcdenoise(pointCloud(vx.Location(thisBar,:)));
    lateralDiagonals{1}{j} = thisBar(idxDenoised);
end

for j = 1:length(lateralDiag2(:,1))
    thisBar = bar(idx(:,1) == lateralDiag2(j));
    [~,idxDenoised,~] = pcdenoise(pointCloud(vx.Location(thisBar,:)));
    lateralDiagonals{2}{j} = thisBar(idxDenoised);
end 

% figure; pcshow(vx.Location, 'g', 'MarkerSize', 50);
% j = 0;
% j = j +1
% hold on; pcshow(vx.Location(cat(1,lateralDiagonals{2}{:}),:), 'w', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(cat(1,lateralDiagonals{1}{:}),:), 'r', 'MarkerSize', 50);
end

