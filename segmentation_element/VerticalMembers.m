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

function [vertical] = VerticalMembers(vx,notAnalysed, vxLocalPca, histWidth, verticalBarWidth)
%% Vertical Members:
% Select voxels with vertical dispersions that were discarded
% perviously.
% Generate a histogram in the X-axis
% Obtain the max points in the histogram that are equidistant from
% one another.
% Select the voxels that are in the neighbourhood of said max points
% using the bar width as guideline.
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

minEigenvector = [NaN,NaN,0.9,NaN,NaN,NaN,NaN,NaN,NaN];
minEigenvalue  = [NaN,NaN,NaN];
maxEigenvector = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
maxEigenvalue  = [NaN,NaN,NaN];

verticals = (~any(abs(vxLocalPca.eigenvectors) <= minEigenvector,2) & ~any(vxLocalPca.eigenvalues <= minEigenvalue,2) ...
            & ~any(abs(vxLocalPca.eigenvectors) >= maxEigenvector,2) & ~any(vxLocalPca.eigenvalues >= maxEigenvalue,2));

verticals = verticals & notAnalysed;

% figure; pcshow(vx.Location,'w', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(notAnalysed,:),'b', 'MarkerSize', 50);
% hold on; pcshow(vx.Location(verticals,:),'r', 'MarkerSize', 50);

% figure; histogram(vx.Location(verticals,1),'BinWidth', histWidth); title('X', histWidth);

[numX, edgesX] = histcounts(vx.Location(verticals,1),'BinWidth', histWidth);
[peaksEdgesX,~] = EquidistantHistPeaks(numX, edgesX, vx.grid/10);

vertical = cell(1,length(peaksEdgesX));
for j = 1:length(peaksEdgesX)
    vertical{j} = find(vx.Location(:,1) > (peaksEdgesX(j)+0.5*histWidth) - 0.5*verticalBarWidth ...
                     & vx.Location(:,1) < (peaksEdgesX(j)+0.5*histWidth) + 0.5*verticalBarWidth ...
                     & verticals);
end

% hold on; pcshow(vx.Location(cat(1,vertical{:}),:),'r', 'MarkerSize', 50);

end

