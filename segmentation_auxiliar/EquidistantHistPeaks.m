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

function [peaksEdges,distances] = EquidistantHistPeaks(varargin)
%% EquidistanHistPeaks
%
% This method looks for the maximum number of equidistance peaks on a
% histogram by analysing the variance of the distances between them.
%
%--------------------------------------------------------------------------
% INPUTS:
% num: numeric Nx1. First output of histcounts.
%
% edges: numeric Nx1. Second output of histcounts.
% 
% maxVariance: numeric 1x1. Maximum variance for iteration to stop.
%
%--------------------------------------------------------------------------
% OUTPUT:
% peaksEdges: numeric Mx1. Edges of the selected peaks.
%
% distances: numeric (M-1)x1. Distance between peaksEdges.
%
%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 21/04/21

%% Checking inputs
parser = inputParser;
parser.CaseSensitive = true; 
parser.addRequired('num', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,NaN]}));
parser.addRequired('edges', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,length(varargin{1})+1]}));
parser.addRequired('maxVariance', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]})); 
parser.parse(varargin{:});

num = parser.Results.num;
edges = parser.Results.edges;         
maxVariance = parser.Results.maxVariance;

%%
% Solving findpeaks problems with first and last values
num = [0,num,0];
% Local max
[peaks,locs] = findpeaks(num, 'MinPeakHeight', 1);
locs = locs-1; % Solving findpeaks problems with first and last values

continue_loop=true;
while continue_loop
    continue_loop=false;

    peaksEdges = edges(locs);
    
    % variance of the distancees between each peak and their next peak 
    if var(peaksEdges(2:end) - peaksEdges(1:end-1)) > maxVariance
        continue_loop=true;
    end

    % Updating values removing the min peak
    locs = locs(peaks > min(peaks));
    peaks = peaks(peaks > min(peaks));
    
end

distances = peaksEdges(2:end) - peaksEdges(1:end-1);
end

