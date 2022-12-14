function [location,idx] = loadCsv(varargin)
% Load the file.csv with assigning an index for each member

%--------------------------------------------------------------------------
% INPUTS:
%
% file : char or string. File of the .csv
%
%--------------------------------------------------------------------------
% OUTPUTS
%
% location : Nx3 numeric. Matrix with all the coordinates XYZ of all the
%            poitns in the file.
%
% idx : Nx1 numeric. Array with the member index of each point in location.
%
%--------------------------------------------------------------------------
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 05/05/2021

%% Checking inputs
parser = inputParser;
parser.CaseSensitive = true; 
parser.addRequired('file',@(x)validateattributes(x,{'string', 'char'},{'nonempty'}));
parser.parse(varargin{:});

file = parser.Results.file;

%% Loading .csv
file = table2array(readtable(file));

%% matrix where the data will be reshaped
location = zeros(numel(file)/3,3);
idx = zeros(numel(file)/3,1);

%% Concatening the table saving an index for each element
rowsForElement = size(file,1);
for i = 1:size(file,2)/3
    location(rowsForElement*(i-1)+1:rowsForElement*i,:) =  file(:,(i-1)*3+1:(i-1)*3+1+2);
    idx(rowsForElement*(i-1)+1:rowsForElement*i) = i;
end

idx = idx(~any(isnan(location),2));
location = location(~any(isnan(location),2),:);


% i = 0;
% figure;
% 
% i = i+1;
% if mod(i,2) ==0
%     color = 'r';
% else
%     color = 'g';
% end
% hold on; pcshow(location(idx == i,:), color,'MarkerSize', 100);

