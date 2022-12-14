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

