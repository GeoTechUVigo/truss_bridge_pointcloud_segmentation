%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
Copyright (C) 2023 Andrés Justo Domínguez andres.justo.dominguez@uvigo.gal
This file is part of the program Automated instance and semantic 
segmentation of point clouds of large metallic truss bridges with modelling
purposes.
The program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) 
any later version.
The program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
more details.
You should have received a copy of the GNU General Public License along 
with the program in COPYING. If not, see <https://www.gnu.org/licenses/>. 
%}

function [num,edges] = HistogramPeaks(data, BinWidth)
%% HistogramPeaks:
% Return the num and edges of the constructed histogram with data and
% BinWidth, but only those that are local maximums.
% Edge values are not the ones returned by histcounts, but the midpoint of
% the bar (edge(i) + 0.5*BinWidth).

%--------------------------------------------------------------------------
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 16/09/2022

% Histogram
[num, edges] = histcounts(data,'BinWidth',BinWidth);

% Add values to detect peaks in first and last bars
num = [0,num,0];

% Detect peaks
[num, locs] = findpeaks(num, 'MinPeakHeight', 1);
locs = locs-1;

edges = edges(locs) + 0.5*BinWidth; 

end

