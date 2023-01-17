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

function [trimmedEdges] = TrimHist(varargin)
% TrimHist
% Method to trim a histogram based on the percentual increase between
% iterations. the histogram increases from its maximum to both sizes until
% in one iteration the percentual increase is lower than threshold or all
% the histrogram have been added. If the limit of a lateral is reached, 
% it only grews in the other direction.
%
%--------------------------------------------------------------------------
% INPUTS:
% num: numeric Nx1. First output of histcounts.
%
% edges: numeric Nx1. Second output of histcounts.
%
% threshold:  numeric 1x1. Minimum increased between iteration to stop.
% 
%--------------------------------------------------------------------------
% OUTPUTS:
% trimmedEdges: numerci 2x1. First and last edge.
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
parser.addRequired('threshold', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]})); 
parser.parse(varargin{:});

num = parser.Results.num;
edges = parser.Results.edges;         
threshold = parser.Results.threshold;         

%% Maximum in the histogram
[~,maxLoc] = max(num);
i = 1;
j = 1;
endLeft = false;
endRight = false;

lastOcupation = sum(num(maxLoc))/sum(num);
while true
    
    %% Updating the ocupation
    currentOcupation = sum(num(maxLoc-i:maxLoc+j))/sum(num);
    
    %% Checking if the increase of the ocupation is enought
    if currentOcupation - lastOcupation < threshold 
        break
    end
    lastOcupation = currentOcupation;
    
    %% Checking laterals of the histogram to determine its grew
    if maxLoc-i > 1
        i = i+1;
    else
        endLeft = true;
    end
    if maxLoc+j < length(num)
        j = j+1;
    else
        endRight = true;
    end
    
    %% Checking if all the histogram has been added
    if endRight && endLeft
        break
    end
end

%% Deleting the last bars if the while has end because of the threshold
if currentOcupation - lastOcupation < threshold
    if ~endRight
        j = j-1;
    end
    if ~endLeft
        i = i-1;
    end    
end

trimmedEdges = [edges(maxLoc-i),edges(maxLoc+j+1)];

end


