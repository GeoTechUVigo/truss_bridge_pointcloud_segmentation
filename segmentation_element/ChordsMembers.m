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

function [horizontal] = ChordsMembers(vx,notAnalysed, vxLocalPca, histWidth, chordWidth)
%% Chords Members:
% Select voxels with horizontal dispersions that were discarded perviously.
% Generate a histogram in the Z-axis.
% Obtain the max points in the histogram that are equidistant.
% Select the voxels that are in the neighbourhood of said max points
% using the bar width as guideline.
% Clusters are denoised in Y-axis using trimHist().
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
% figure; pcshow(vx.Location,'w', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(notAnalysed,:),'g', 'MarkerSize', 50);

minEigenvectorMastNoBracket = [0.9,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
minEigenvalueMastNoBracket  = [NaN,NaN,NaN];
maxEigenvectorMastNoBracket = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
maxEigenvalueMastNoBracket  = [NaN,NaN,NaN];

horizontals = (~any(abs(vxLocalPca.eigenvectors) <= minEigenvectorMastNoBracket,2) & ~any(vxLocalPca.eigenvalues <= minEigenvalueMastNoBracket,2) ...
            & ~any(abs(vxLocalPca.eigenvectors) >= maxEigenvectorMastNoBracket,2) & ~any(vxLocalPca.eigenvalues >= maxEigenvalueMastNoBracket,2));

horizontals = horizontals & notAnalysed;

% figure; pcshow(vx.Location,'w', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(notAnalysed,:),'b', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(horizontals,:),'r', 'MarkerSize', 50);

% figure; histogram(vx.Location(horizontals,3),'BinWidth', histWidth); title('X', histWidth);

[num, edges] = histcounts(vx.Location(horizontals,3),'BinWidth', histWidth);
[peaksEdges,~] = EquidistantHistPeaks(num, edges, vx.grid/10);

horizontal = cell(1,length(peaksEdges));
for j = 1:length(peaksEdges)
    %% Z limits of this element
    element = find(vx.Location(:,3) > (peaksEdges(j)+0.5*histWidth) - 0.5*chordWidth ...
                       & vx.Location(:,3) < (peaksEdges(j)+0.5*histWidth) + 0.5*chordWidth ...
                       & horizontals);
    %% denoising 
    [~,aux,~] = pcdenoise(pointCloud(vx.Location(element,:)));
%     figure; pcshow(vx.Location(element,:), 'r', 'MarkerSize',50);
%     hold on; pcshow(vx.Location(element(aux),:), 'g', 'MarkerSize',50);
    
    %% Orienting
    element = element(aux);
    elementLocation = vx.Location(element,:) * PcaFlattering(vx.Location(element,:));
    elementLocation = elementLocation - mean(elementLocation);
    
%     figure; pcshow(elementLocation, 'g', 'MarkerSize',50);
%     figure; histogram(elementLocation(:,2),'BinWidth', vx.grid);               
     
    %% Y limits
    [num, edges] = histcounts(elementLocation(:,2),'BinWidth', vx.grid);
    trimmedEdges = TrimHist(num, edges, 0.05);
    
    
    %% Saving
    horizontal{j} = element(elementLocation(:,2) > trimmedEdges(1) ...
                          & elementLocation(:,2) < trimmedEdges(2));
   
%     hold on; pcshow(vx.Location(horizontal{j},:),'g', 'MarkerSize', 50);
     
end

% hold on; pcshow(vx.Location(cat(1,horizontal{:}),:),'r', 'MarkerSize', 50);

end