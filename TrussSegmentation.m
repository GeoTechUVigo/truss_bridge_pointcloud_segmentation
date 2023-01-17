%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
This file is part of the program Automated instance and semantic 
segmentation of point clouds of large metallic truss bridges with modelling
purposes.
The program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or any later version.
The program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details.
You should have received a copy of the GNU General Public License along 
with the program in COPYING. If not, see <https://www.gnu.org/licenses/>. 
%}

function [] = TrussSegmentation(varargin)
%% Summary
%
% Truss segmentation for box-like trusses.
% Load and voxelize the cloud.
% Denoise the cloud and reorient it.
% Truss is split into its lateral and horizontal faces.
% Analyze the lateral and horizontal faces in an isolated manner.
% With the information of the vertical members, the inner faces are located
% and analysed in an isolated manner.
% Lateral Faces:
    % LocalPCA
    % Lateral Diagonals:
        % Select voxels with linear dispersion that are neither too vertical nor too horizontal.
        % Group neighbouring voxels.
        % Calculate predominant directions of the groups.
        % Prolong the groups following a main direction and add to its
        % group any other that falls into the prolongation and shares
        % its direction.
        % Erase small clusters.
        % Erase bars that are outside of the face plane.
        % Denoise bars.
    % Vertical Members:
        % Select voxels with vertical dispersions that were discarded
        % perviously.
        % Generate a histogram in the X-axis
        % Obtain the max points in the histogram that are equidistant from
        % one another.
        % Select the voxels that are in the neighbourhood of said max points
        % using the bar width as guideline.
    % Chords:
        % Same procedure as with vertical members but searching for
        % horizontal dispersion and analyzing the histogram in the Z-axis.
        % Besides, the chords are denoised in Y-axis using trimHist().
% Horizontal Faces:
    % Rotate 90 degrees to match the lateral faces positioning.
    % Repeat the procedure for lateral faces but without the chord
    % analysis.
% Inner faces:
% Computed which vertical members are at the same level in the bridge.
% Calcalating its centre, the position of the inner faces is known.
% For each inner faces:
    % Rotate 90 degrees to match the lateral faces positioning.
    %  Repeat the procedure for lateral faces to segmented the diagonals,
    % specifying the direction by computing the angle, knowing that these
    % bars go through the entire face.
    % Horizontal member is segmented selecting voxels with horizontal
    % dispersion, selecting the voxels in the middel of the face and
    % filtering by Z and Y taking into account the width of the bar.
%
% Finally, the coordinates of each segmented bar are save in a .csv with
% the decimals specified. The header of each column specifies:
% Truss, face, type, id, dimension.
% In the cases of the inner faces, the id of their vertical members is also
% specified.
%
%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 02/08/2021

%% Checking inputs
parser = inputParser;
parser.CaseSensitive = true;

parser.addRequired('path_in', @(x)validateattributes(x,{'string', 'char'}, {}));
parser.addRequired('path_out', @(x)validateattributes(x,{'string', 'char'}, {}));
parser.addRequired('width_faces_v_h_i', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,3]}));
parser.addRequired('lateralBarWidth', @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addRequired('verticalBarWidth', @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addRequired('chordWidthZ', @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('grid', 0.05, @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('inside_board', false, @(x)validateattributes(x,{'logical'}, {}));
parser.addOptional('maxDistCentralFace', 0.5, @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('numHorizontalFaces', 2, @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('numVerticalFaces', 2, @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('x_limits', [-inf, inf],@(x)validateattributes(x,{'numeric'}, {}));
parser.addOptional('alignmentFile', [],@(x)validateattributes(x,{'string', 'char'}, {}));
parser.addOptional('plot_vx', false, @(x)validateattributes(x,{'logical'}, {}));
parser.addOptional('plot', false, @(x)validateattributes(x,{'numeric'}, {'real','nonnan'}));
parser.addOptional('inner_vertical', true, @(x)validateattributes(x,{'logical'}, {}));
parser.addOptional('inner_horizontal', true, @(x)validateattributes(x,{'logical'}, {}));
parser.addOptional('path_out_csv', [], @(x)validateattributes(x,{'string', 'char'}, {}));

parser.parse(varargin{:});

path_in = parser.Results.path_in;
path_out = parser.Results.path_out;         
width_faces_v_h_i = parser.Results.width_faces_v_h_i;
lateralBarWidth = parser.Results.lateralBarWidth;
verticalBarWidth = parser.Results.verticalBarWidth;
chordWidthZ = parser.Results.chordWidthZ;
grid = parser.Results.grid;
inside_board = parser.Results.inside_board;
maxDistCentralFace = parser.Results.maxDistCentralFace;  % max distance in X of 2 vertical members to be consider part of the same central face
numHorizontalFaces = parser.Results.numHorizontalFaces;
numVerticalFaces = parser.Results.numVerticalFaces;
x_limits = parser.Results.x_limits;
alignmentFile = parser.Results.alignmentFile;
plot_vx = parser.Results.plot_vx;
plot = parser.Results.plot;
inner_vertical = parser.Results.inner_vertical;
inner_horizontal = parser.Results.inner_horizontal;
path_out_csv = parser.Results.path_out_csv;

%% Add paths
addpath classes segmentation_auxiliar segmentation_element utils

%% Loading and voxeling cloud
cloudOrigin = LASread(path_in);
vxOrigin = Voxels(pointCloud_([cloudOrigin.record.x , cloudOrigin.record.y, cloudOrigin.record.z]),grid);
clear cloudOrigin;

% Denoising
[~,vxIdx,~] = pcdenoise(pointCloud(vxOrigin.Location));
% figure; pcshow(vxOrigin.Location, 'r', 'MarkerSize', 50);
% hold on; pcshow(vxOrigin.Location(vxIdx,:), 'g', 'MarkerSize', 50);
vxIdx = vxIdx';
vx = select(vxOrigin,vxIdx);

%% Reorienting the point cloud
if ~isempty(alignmentFile)

    % Loading alignment and selecting the alignment of this cloud
    alignmentOrigin = csvread(alignmentFile);
    
    alignmentIdx = all(alignmentOrigin(:,1:2) >= min(vx.Location(:,1:2)),2) & all(alignmentOrigin(:,1:2) < max(vx.Location(:,1:2)),2); % idx of alignment in the bounding box of vx
    if find(alignmentIdx,1,'first') > 1 % Adding one point before
        alignmentIdx(find(alignmentIdx,1,'first')-1) = true;
    end
    if find(alignmentIdx,1,'last') < length(alignmentIdx) % Adding one point after
        alignmentIdx(find(alignmentIdx,1,'last')+1) = true;
    end
    
    alignmentOrigin = alignmentOrigin(alignmentIdx,:);
        
%     figure; pcshow(vxOrigin.Location, 'g', 'MarkerSize', 50);
%     hold on; pcshow(alignmentOrigin, 'r', 'MarkerSize', 100);
    
    % Orienting
    meanAlignmentOrigin = mean(alignmentOrigin);
    alignment = alignmentOrigin - meanAlignmentOrigin;
    pcaAlignment = PcaFlattering(alignment);
    alignment    = alignment * pcaAlignment;
    
    % Checking if alignment is in the the correct order
    if alignment(end,1) < alignment(1,1)
        pcaAlignment = RotateAxes(pcaAlignment, 180, pcaAlignment(:,3));
        alignment    = alignmentOrigin - meanAlignmentOrigin;
        alignment    = alignment * pcaAlignment;
    end

else
    meanAlignmentOrigin = mean(vx.Location);
    pcaAlignment = PcaFlattering(vx.Location);
end

vx.Location = vx.Location - meanAlignmentOrigin;
vx.Location = vx.Location * pcaAlignment;

% figure; pcshow(vx.Location,'g','MarkerSize',50);
% hold on; pcshow(alignment,'w','MarkerSize',200);
% hold on; pcshow(alignment(1,:),'r','MarkerSize',200);

%% Looking for the faces
% figure; histogram(vx.Location(:,1),'BinWidth', lateralBarWidth*3); title('X')
% figure; histogram(vx.Location(:,2),'BinWidth', lateralBarWidth*3); title('Y')
% figure; histogram(vx.Location(:,3),'BinWidth', lateralBarWidth*2); title('Z')

% Histogram peaks analysis
[numY, edgesY] = HistogramPeaks(vx.Location(:,2),width_faces_v_h_i(1)/3); % Y
[numZ, edgesZ] = HistogramPeaks(vx.Location(:,3),width_faces_v_h_i(2)/3); % Z

%%
% Y limits
[~,order] = sort(numY,'descend');
edgesY = edgesY(order);

% Z limits
[~,order] = sort(numZ,'descend');
edgesZ = edgesZ(order);

facesLocationY = edgesY(1:numVerticalFaces);

if inside_board
    first=2;
else
    first=1;
end

facesLocationZ = edgesZ(first:first+numHorizontalFaces-1);
deck_z= edgesZ(1);

%% Sorting the faces
facesLocationY = sort(facesLocationY,'ascend');
facesLocationZ = sort(facesLocationZ,'ascend');

%% Analysing vertical faces
for i = 1:numVerticalFaces
    idxSec = find(vx.Location(:,2) > (facesLocationY(i) - 0.5*width_faces_v_h_i(1)) & vx.Location(:,2) < (facesLocationY(i) + 0.5*width_faces_v_h_i(1)) ...
                & vx.Location(:,3) > min(facesLocationZ - 0.5*width_faces_v_h_i(2)) & vx.Location(:,3) < max(facesLocationZ + 0.5*width_faces_v_h_i(1)) ...
                & vx.Location(:,1) > x_limits(1)                                  & vx.Location(:,1) < x_limits(2));
    vxSec = select(vx, idxSec);
   
%     figure;pcshow(vxSec.Location)
    
    %% Local PCA analysis
    
    vxLocalPca = PcaLocal(vxSec,6*grid,1);
    
    %% Lateral Diagonals
    minPointsSeed = 1.5/vxSec.grid;
    toleranceDirection = 0.97;
    minPointsBar = 0; %2.0/vx.grid;
    
    localPcaThreshold{1}{1} = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]; % min eigenvectors.
    localPcaThreshold{1}{2} = [0.9,NaN,0.9,NaN,NaN,NaN,NaN,NaN,NaN]; % max eigenvectors.
    localPcaThreshold{2}{1} = [0.7,NaN,NaN]; % min eigenvalues.
    localPcaThreshold{2}{2} = [NaN,NaN,NaN]; % max eigenvalues.

    components.verticalFaces{i}.lateralDiagonals = LateralDiagonals(vxSec, vxLocalPca, minPointsSeed, toleranceDirection, lateralBarWidth, minPointsBar, localPcaThreshold);
        
%     figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%     j = 0;
%     j = j+1
%     hold on; pcshow(vxSec.Location(components.verticalFaces{i}.lateralDiagonals{1}{j},:),'r','MarkerSize', 50); 

    %% Vertical members
    histWidth = 2*vxSec.grid;
    notAnalysed = true(length(vxSec.Location),1);
    if ~isempty(components.verticalFaces{i}.lateralDiagonals{1}) || ~isempty(components.verticalFaces{i}.lateralDiagonals{2})
        analysed = [cat(1,components.verticalFaces{i}.lateralDiagonals{1}{:});cat(1,components.verticalFaces{i}.lateralDiagonals{2}{:})];
        notAnalysed(analysed) = false;
    end
    
%     figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%     hold on; pcshow(vxSec.Location(notAnalysed,:),'g','MarkerSize', 50); 

    components.verticalFaces{i}.verticals = VerticalMembers(vxSec,notAnalysed, vxLocalPca, histWidth, verticalBarWidth);
    
%      figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%      hold on; pcshow(vxSec.Location(cat(1,components.verticalFaces{i}.verticals{:}),:),'r','MarkerSize', 50); 

   %% Chords
   histWidth = 2*vxSec.grid;
   if ~isempty(cat(1,components.verticalFaces{i}.verticals{:}))
        notAnalysed(cat(1,components.verticalFaces{i}.verticals{:})) = false;
   end
   
%    figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%    hold on; pcshow(vxSec.Location(notAnalysed,:),'g','MarkerSize', 50); 

   components.verticalFaces{i}.chords = ChordsMembers(vxSec,notAnalysed, vxLocalPca, histWidth, chordWidthZ);
   
%    figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%    hold on; pcshow(vxSec.Location(cat(1,components.verticalFaces{i}.chords{:}),:),'r','MarkerSize', 50); 

    %% Indexes in vxOring, not in vxSec
    % Diagonals
    for j = 1:length(components.verticalFaces{i}.lateralDiagonals)
        for k = 1:length(components.verticalFaces{i}.lateralDiagonals{j})
            components.verticalFaces{i}.lateralDiagonals{j}{k} = vxIdx(idxSec(components.verticalFaces{i}.lateralDiagonals{j}{k}));
        end
    end
    
    % Vertical members
    for j = 1:length(components.verticalFaces{i}.verticals)
        components.verticalFaces{i}.verticals{j} = vxIdx(idxSec(components.verticalFaces{i}.verticals{j}));
    end
    % Chords
    for j = 1:length(components.verticalFaces{i}.chords)
        components.verticalFaces{i}.chords{j} = vxIdx(idxSec(components.verticalFaces{i}.chords{j}));
    end
    
%     figure; pcshow(vxOrigin.Location, 'w', 'MarkerSize', 50);
%     figure; pcshow(vxOrigin.Location(vxIdx(idxSec),:), [0,0,0], 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.lateralDiagonals{1}{:}),:), 'r', 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.lateralDiagonals{2}{:}),:), 'g', 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.verticals{:}),:), 'b', 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.chords{:}),:), 'y', 'MarkerSize', 50);
%     WhitePcshow();
end

%% Horizontal sections
if ~inside_board
    % The top face is the board
    numHorizontalFaces = numHorizontalFaces - 1;
end
for i = 1:numHorizontalFaces
    
    idxSec = find(vx.Location(:,3) > (facesLocationZ(i) - 0.5*width_faces_v_h_i(2)) & vx.Location(:,3) < (facesLocationZ(i) + 0.5*width_faces_v_h_i(2)) ...
                & vx.Location(:,2) > min(facesLocationY - 0.5*width_faces_v_h_i(1)) & vx.Location(:,2) < max(facesLocationY + 0.5*width_faces_v_h_i(1)) ...
                & vx.Location(:,1) > x_limits(1)                                    & vx.Location(:,1) < x_limits(2));
            
    vxSec = select(vx, idxSec);
    
    % Rotate the section to work as it is a vertical section
    vxSec.Location = RotateAxes(vxSec.Location', 90, [1,0,0]')';
%     figure; pcshow(vxSec.Location, 'g', 'MarkerSize',50);

    %% Local PCA analysis
    
    vxLocalPca = PcaLocal(vxSec,6*grid,1);

    %% Lateral Diagonals
    minPointsSeed = 1/vxSec.grid;
    toleranceDirection = 0.97;
    minPointsBar = 0; %2.0/vx.grid;
    
    localPcaThreshold{1}{1} = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
    localPcaThreshold{1}{2} = [0.9,NaN,0.9,NaN,NaN,NaN,NaN,NaN,NaN];
    localPcaThreshold{2}{1} = [0.7,NaN,NaN];
    localPcaThreshold{2}{2} = [NaN,NaN,NaN];

    components.horizontalFaces{i}.lateralDiagonals = LateralDiagonals(vxSec, vxLocalPca, minPointsSeed, toleranceDirection, lateralBarWidth, minPointsBar, localPcaThreshold);
        
%     figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
% 
%     j = 0;
%     j = j+1
%     hold on; pcshow(vxSec.Location(components.horizontalFaces{i}.lateralDiagonals{1}{j},:),'r','MarkerSize', 50);

    %% Vertical members
    histWidth = 2*vxSec.grid;
    notAnalysed = true(length(vxSec.Location),1);
    analysed = [];
    if ~isempty(components.horizontalFaces{i}.lateralDiagonals{1}) 
        analysed = cat(1,components.horizontalFaces{i}.lateralDiagonals{1}{:});
    end
    if ~isempty(components.horizontalFaces{i}.lateralDiagonals{2})
        analysed = [analysed;cat(1,components.horizontalFaces{i}.lateralDiagonals{2}{:})];
    end
    notAnalysed(analysed) = false;
    
%     figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%     hold on; pcshow(vxSec.Location(notAnalysed,:),'g','MarkerSize', 50); 

    components.horizontalFaces{i}.verticals = VerticalMembers(vxSec,notAnalysed, vxLocalPca, histWidth, verticalBarWidth);
    
%      figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
%      hold on; pcshow(vxSec.Location(cat(1,components.horizontalFaces{i}.verticals{:}),:),'r','MarkerSize', 50); 

    %% Indexes in vxOring, not in vxSec
    % Diagonals
    for j = 1:length(components.horizontalFaces{i}.lateralDiagonals)
        for k = 1:length(components.horizontalFaces{i}.lateralDiagonals{j})
            components.horizontalFaces{i}.lateralDiagonals{j}{k} = vxIdx(idxSec(components.horizontalFaces{i}.lateralDiagonals{j}{k}));
        end
    end
    
    % Vertical members
    for j = 1:length(components.horizontalFaces{i}.verticals)
        components.horizontalFaces{i}.verticals{j} = vxIdx(idxSec(components.horizontalFaces{i}.verticals{j}));
    end
    
%     figure; pcshow(vxOrigin.Location, 'w', 'MarkerSize', 50);
%     figure; pcshow(vxOrigin.Location(vxIdx(idxSec),:), [0,0,0], 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.lateralDiagonals{1}{:}),:), 'r', 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.lateralDiagonals{2}{:}),:), 'g', 'MarkerSize', 50);
%     hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.verticals{:}),:), 'b', 'MarkerSize', 50);
%     WhitePcshow();

end

%% Inner Faces
%% Inner Faces locations
% Check if there is any central face without vertical member is all the
% vertical faces.
% Detect the face using the vertical members. The mean X of each vertical
% face is the X coordinate of the plane of each Central face.
facesLocationXAux = [];
for i = 1:numel(components.verticalFaces)
    facesLocationXAux{i} = zeros(numel(components.verticalFaces{1}.verticals),1);
    for j = 1:numel(components.verticalFaces{i}.verticals)
        facesLocationXAux{i}(j) = mean(vx.Location(ismember(vxIdx,components.verticalFaces{i}.verticals{j}),1));
    end
end

% Analyse only central faces with vertical members in all the faces.
% this vector has the idx of vertical faces pairs in each face.
innerRelVertical = zeros(numel(facesLocationXAux{i}), numel(facesLocationXAux));
innerRelVertical(:,1) = 1:numel(innerRelVertical(:,1)); % first face
% Compare all the elements the first face with the others
for i = 1:numel(facesLocationXAux{1})
    % All the other sections
    for j = 2:numel(facesLocationXAux)
        % All the elements in this section
        for k = 1:numel(facesLocationXAux{j})
            if abs(facesLocationXAux{1}(i) - facesLocationXAux{j}(k)) < maxDistCentralFace
                innerRelVertical(i,j) = k;
                break;
            end
        end
    end
end

% Only faces with all its vertical members.
innerRelVertical = innerRelVertical(all(innerRelVertical ~=0, 2),:);

% Calculate facesLocationX as the middle of all its vertical members.
facesLocationX = zeros(size(innerRelVertical,1),1);
for i = 1:size(innerRelVertical,1)
    allVertical = [];
    for j = 1:size(innerRelVertical,2)
        allVertical = [allVertical, facesLocationXAux{j}(innerRelVertical(i,j))];
    end
    facesLocationX(i) = mean(allVertical);
end

%% Inner faces analysis. They are analysed all together.
idxInnerFaces = [];
for i = 1:numel(facesLocationX)
       
    idxSec = find(vx.Location(:,1) > (facesLocationX(i)- 0.5*width_faces_v_h_i(3))    & vx.Location(:,1) < (facesLocationX(i) + 0.5*width_faces_v_h_i(3)) ...
                & vx.Location(:,2) > (min(facesLocationY) + 0.5*width_faces_v_h_i(1)) & vx.Location(:,2) < (max(facesLocationY) - 0.5*width_faces_v_h_i(1)) ...
                & vx.Location(:,3) > min(facesLocationZ)                              & vx.Location(:,3) < deck_z);

    idxInnerFaces{i} = idxSec;

end

idxInnerFacesAll = cat(1,idxInnerFaces{:});
 vxSec = select(vx, idxInnerFacesAll);
 vxSec.Location = RotateAxes(vxSec.Location', 90, [0,0,1]')';
%  figure; pcshow(vxSec.Location, 'g', 'MarkerSize',50);


%% Local PCA analysis
vxLocalPca = PcaLocal(vxSec,6*grid,1);

%% Lateral Diagonals
minPointsSeed = 1.5/vxSec.grid;
toleranceDirection = 0.97;
minPointsBar = 0; %2.0/vx.grid;

localPcaThreshold{1}{1} = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]; % min eigenvectors.
localPcaThreshold{1}{2} = [0.9,NaN,0.9,NaN,NaN,NaN,NaN,NaN,NaN]; % max eigenvectors.
localPcaThreshold{2}{1} = [0.7,NaN,NaN]; % min eigenvalues.
localPcaThreshold{2}{2} = [NaN,NaN,NaN]; % max eigenvalues.

innerFaces.lateralDiagonals = LateralDiagonals(vxSec, vxLocalPca, minPointsSeed, toleranceDirection, lateralBarWidth, minPointsBar, localPcaThreshold,'clean_y',false);
    
% figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(cat(1,innerFaces.lateralDiagonals{1}{:}),:),'r','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(cat(1,innerFaces.lateralDiagonals{2}{:}),:),'g','MarkerSize', 50); 

%% Vertical members
if inner_vertical
histWidth = 2*vxSec.grid;
notAnalysed = true(length(vxSec.Location),1);
if ~isempty(innerFaces.lateralDiagonals{1}) || ~isempty(innerFaces.lateralDiagonals{2})
    analysed = [cat(1,innerFaces.lateralDiagonals{1}{:});cat(1,innerFaces.lateralDiagonals{2}{:})];
    notAnalysed(analysed) = false;
end

% figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(notAnalysed,:),'g','MarkerSize', 50); 

innerFaces.verticals = VerticalMembers(vxSec,notAnalysed, vxLocalPca, histWidth, verticalBarWidth);

% figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(cat(1,innerFaces.verticals{:}),:),'r','MarkerSize', 50);
else
    innerFaces.verticals = {};
end

%% Horizontal member in the middle
if inner_horizontal
% Remove analysed voxels
notAnalysed = true(length(vxSec.Location),1);
analysed = [];
if ~isempty(innerFaces.lateralDiagonals{1}) 
    analysed = cat(1,innerFaces.lateralDiagonals{1}{:});
end
if ~isempty(innerFaces.lateralDiagonals{2})
    analysed = [analysed;cat(1,innerFaces.lateralDiagonals{2}{:})];
end
notAnalysed(analysed) = false;

% Horizontal voxels
minEigenvectorMastNoBracket = [0.9,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
minEigenvalueMastNoBracket  = [NaN,NaN,NaN];
maxEigenvectorMastNoBracket = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
maxEigenvalueMastNoBracket  = [NaN,NaN,NaN];

horizontals = (~any(abs(vxLocalPca.eigenvectors) <= minEigenvectorMastNoBracket,2) & ~any(vxLocalPca.eigenvalues <= minEigenvalueMastNoBracket,2) ...
            & ~any(abs(vxLocalPca.eigenvectors) >= maxEigenvectorMastNoBracket,2) & ~any(vxLocalPca.eigenvalues >= maxEigenvalueMastNoBracket,2));

horizontals = horizontals & notAnalysed;

% Calculate angle
ranges = range(vxSec.Location);

% Select the middle bar in Z
middleZ = min(vxSec.Location(:,3)) + ranges(3)/2;
middlePoints = false(size(horizontals));
middlePoints(abs(vxSec.Location(:,3) - middleZ) < lateralBarWidth*10 & horizontals) = 1;
middlePoints = find(middlePoints);

% Centre of the bar
% Z
histWidth = 6*vxSec.grid;
%     figure; histogram(vxSec.Location(middlePoints,3),'BinWidth', histWidth); title('Z')
[num, edges] = histcounts(vxSec.Location(middlePoints,3),'BinWidth',histWidth); 
[~, maxNum] = max(num);
centre = edges(maxNum) + histWidth/2;
middlePoints = middlePoints(abs(vxSec.Location(middlePoints, 3) - centre) < histWidth);
% Y
% histWidth = lateralBarWidth;
% %     figure; histogram(vxSec.Location(middlePoints,2),'BinWidth', histWidth); title('Y')
% [num, edges] = histcounts(vxSec.Location(middlePoints,2),'BinWidth',histWidth); 
% [~, maxNum] = max(num);
% centre = edges(maxNum) + histWidth/2;
% middlePoints = middlePoints(abs(vxSec.Location(middlePoints, 2) - centre) < histWidth);

% figure; pcshow(vxSec.Location,'w','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(horizontals,:),'g','MarkerSize', 50);
% hold on; pcshow(vxSec.Location(middlePoints,:),'r','MarkerSize', 50);

innerFaces.horizontal{1} = middlePoints;
else
    innerFaces.horizontal = {};
end

%% Split faces and write indexes in vxOrigin
% TODO: seguir aqui. Separar los miembros en cada cara y ver como segmento
% la horizontal del medio.
components.innerFaces = [];
for i = 1:numel(facesLocationX)

    % Diagonals
    components.innerFaces{i}.lateralDiagonals{1}=[];
    for j=1:numel(innerFaces.lateralDiagonals{1})
        element = idxInnerFacesAll(innerFaces.lateralDiagonals{1}{j});

        a = ismember(idxInnerFaces{i}, element);

        if any(a)
            components.innerFaces{i}.lateralDiagonals{1}{length(components.innerFaces{i}.lateralDiagonals{1})+1} = vxIdx(idxInnerFaces{i}(a));
        end

%         figure; pcshow(vxOrigin.Location, 'g');
%         hold on; pcshow(vxOrigin.Location(components.innerFaces{i}.lateralDiagonals{1}{end}, :), 'r', 'MarkerSize',100);
    end

    components.innerFaces{i}.lateralDiagonals{2}=[];
    for j=1:numel(innerFaces.lateralDiagonals{2})
        element = idxInnerFacesAll(innerFaces.lateralDiagonals{2}{j});

        a = ismember(idxInnerFaces{i}, element);

        if any(a)
            components.innerFaces{i}.lateralDiagonals{2}{length(components.innerFaces{i}.lateralDiagonals{2})+1} = vxIdx(idxInnerFaces{i}(a));
        end
    end

    % Vertical members
    components.innerFaces{i}.verticals=[];
    for j=1:numel(innerFaces.verticals)
        element = idxInnerFacesAll(innerFaces.verticals{j});

        a = ismember(idxInnerFaces{i}, element);

        if any(a)
            components.innerFaces{i}.verticals{length(components.innerFaces{i}.verticals)+1} = vxIdx(idxInnerFaces{i}(a));
        end
    end

    % Horizontal members
    components.innerFaces{i}.horizontal=[];
    for j=1:numel(innerFaces.horizontal)
        element = idxInnerFacesAll(innerFaces.horizontal{j});

        a = ismember(idxInnerFaces{i}, element);

        if any(a)
            components.innerFaces{i}.horizontal{length(components.innerFaces{i}.horizontal)+1} = vxIdx(idxInnerFaces{i}(a));
        end
    end
end

%% Plot
if plot_vx
markerSize = 100;
figure; pcshow(vxOrigin.Location, [0.5, 0.5, 0.5], 'MarkerSize', markerSize);
for i = 1:numel(components.verticalFaces)
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.lateralDiagonals{1}{:}),:), 'g', 'MarkerSize', markerSize);
        hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.lateralDiagonals{2}{:}),:), 'g', 'MarkerSize', markerSize);
        hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.verticals{:}),:), 'b', 'MarkerSize', markerSize);
        hold on; pcshow(vxOrigin.Location(cat(1,components.verticalFaces{i}.chords{:}),:),  [0.7,0.7,0.2], 'MarkerSize', markerSize);
    end
end
for i = 1:numel(components.horizontalFaces)
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.lateralDiagonals{1}{:}),:), [0,0,0], 'MarkerSize', markerSize);
        hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.lateralDiagonals{2}{:}),:), [0,0,0], 'MarkerSize', markerSize);
        hold on; pcshow(vxOrigin.Location(cat(1,components.horizontalFaces{i}.verticals{:}),:), [0.5,0,0.9], 'MarkerSize', markerSize);
    end
end
for i = 1:numel(components.innerFaces)
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.innerFaces{i}.lateralDiagonals{1}{:}),:), 'r', 'MarkerSize', markerSize);
    end
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.innerFaces{i}.lateralDiagonals{2}{:}),:), 'r', 'MarkerSize', markerSize);
    end
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.innerFaces{i}.verticals{:}),:), 'b', 'MarkerSize', markerSize);
    end
    try
        hold on; pcshow(vxOrigin.Location(cat(1,components.innerFaces{i}.horizontal{:}),:), [0.7,0.5,0], 'MarkerSize', markerSize);
    end
end
WhitePcshow();
end

%% Plot in the not voxelised cloud
if plot
markerSize = 100;
figure; pcshow(vxOrigin.parent_cloud(), [0.5,0.5,0.5]);
for i = 1:numel(components.verticalFaces)
    try
        element = vxOrigin.parent_idx((cat(1,components.verticalFaces{i}.lateralDiagonals{1}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), 'g', 'MarkerSize', markerSize);

        element = vxOrigin.parent_idx((cat(1,components.verticalFaces{i}.lateralDiagonals{2}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), 'g', 'MarkerSize', markerSize);


        element = vxOrigin.parent_idx((cat(1,components.verticalFaces{i}.verticals{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), 'b', 'MarkerSize', markerSize);

        element = vxOrigin.parent_idx((cat(1,components.verticalFaces{i}.chords{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0.7,0.7,0], 'MarkerSize', markerSize);
    end
end
for i = 1:numel(components.horizontalFaces)
    try
        element = vxOrigin.parent_idx((cat(1,components.horizontalFaces{i}.lateralDiagonals{1}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0,0,0], 'MarkerSize', markerSize);

        element = vxOrigin.parent_idx((cat(1,components.horizontalFaces{i}.lateralDiagonals{2}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0,0,0], 'MarkerSize', markerSize);
        
        element = vxOrigin.parent_idx((cat(1,components.horizontalFaces{i}.verticals{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0.5,0,0.9], 'MarkerSize', markerSize);
    end
end
for i = 1:numel(components.innerFaces)
    try
        element = vxOrigin.parent_idx((cat(1,components.innerFaces{i}.lateralDiagonals{1}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), 'r', 'MarkerSize', markerSize);
      
        element = vxOrigin.parent_idx((cat(1,components.innerFaces{i}.lateralDiagonals{2}{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), 'r', 'MarkerSize', markerSize);
    end
    try
        element = vxOrigin.parent_idx((cat(1,components.innerFaces{i}.vertical{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0.7,0.5,0], 'MarkerSize', markerSize);
    end
    try
        element = vxOrigin.parent_idx((cat(1,components.innerFaces{i}.horizontal{:})));
        element = cat(1,element{:});
        hold on; pcshow(vxOrigin.parent_cloud(element,:), [0.7,0.5,0], 'MarkerSize', markerSize);
    end
end
WhitePcshow();
end

%% Save original point cloud
SaveLas(vxOrigin.parent_idx, components, path_in, path_out);

%% Save in csv
if ~isempty(path_out_csv)
TrussCsv(components, vxOrigin, path_out_csv, innerRelVertical, 3);
end

end

