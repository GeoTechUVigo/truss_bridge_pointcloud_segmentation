function [vectors] = RotateAxes(varargin)
% Function to rotate column vector around the axi rotationAxis.
%
%--------------------------------------------------------------------------
% INPUTS:
%
% vectors : numeric. 3D column vectors
%
% degrees : numeric. Rotation angle in degrees
% 
% rotationAxis : numeric. 3D column vector. The rotation vector
%
%--------------------------------------------------------------------------
% OUTPUTS:
%
% VECTORS : numeric. 3D column vectors rotated
%
%--------------------------------------------------------------------------
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 30/04/2021

%% Checking inputs
parser = inputParser;
parser.CaseSensitive = true; 
parser.addRequired('vectors', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [3,NaN]}));
parser.addRequired('degrees', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [1,1]}));
parser.addRequired('rotationAxis', @(x)validateattributes(x,{'numeric'}, {'real','nonnan','size', [3,1]})); 
parser.parse(varargin{:});

vectors = parser.Results.vectors;
degrees = parser.Results.degrees;         
rotationAxis = parser.Results.rotationAxis;   

%% 
rotationMatrix = makehgtform('axisrotate',rotationAxis,degrees*pi/180); % tranform matrix 4x4. Rotation degrees in the axis rotationAxis
rotationMatrix = rotationMatrix(1:3,1:3); % rotation matrix 3x3
vectors        = rotationMatrix * vectors; % Rotation operation of column vectors

end

