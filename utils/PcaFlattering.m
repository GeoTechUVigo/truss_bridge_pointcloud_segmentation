function [mpca] = PcaFlattering(location)
% This function modifies the result of pca(location).
% The first eigenvector is not modified. The second is orthonormal to the
% first but in XY plane. The thrid is orthonormal to both and pointing
% upwards (its Z component is positive)
%
%--------------------------------------------------------------------------
% INPUTS:
%
% location : Nx3 numeric.
%                
% -------------------------------------------------------------------------
% OUTPUTS:
% 
% mpca : 3x3 numeric.
%                           
% -------------------------------------------------------------------------
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 23/12/2020

%%
[cloudPca,~,~] = pca(location, 'Economy', false); % normal pca
mpca = zeros(size(cloudPca));
mpca(:,1) = cloudPca(:,1); % principal eigenvector
mpca(:,2) = [-cloudPca(2,1), cloudPca(1,1), 0] / norm([-cloudPca(2,1), cloudPca(1,1), 0]); % orthonormal to pca(1) and in the xy plane
mpca(:,3) = cross(mpca(:,1), mpca(:,2)); % orthonormal to mpca(1) and mpca(2)

rot90Pca2= makehgtform('axisrotate',mpca(:,2),90*pi/180); % rotate 90 degrees in pca(2) 4x4 matrix
rot90Pca2 = rot90Pca2(1:3,1:3); % 3x3 matrix
mpca(:,3) = mpca(:,1)' * rot90Pca2; % orthonormal to mpca(1) and mpca(2)

if mpca(3,3) < 0 % if mpca(:,3) is not pointing upwards (its Z component is not positive)
    mpca = RotateAxes(mpca, 180, mpca(:,1));
end
end

