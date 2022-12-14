function [directions] = PcaLocal(vx,maxDistance,mode)
% This fuction calculates the principal components of each voxels and its
% closest ones. For that, pca() is applied to each voxel and the voxels 
% closer than maxDistance.
%
% In mode == 2 it is the same but not only considering distance but also 
% they have to be a neighbourhood propagation.
% 
%--------------------------------------------------------------------------
% INPUTS:
%
% vx : Voxels. Cloud 
% 
% maxDistance : numeric. Max distance between a points and the points that
%               are considered to do th pca() of this voxel.
% 
% mode : numeric. If 1 --> considering distances.
%                 2 --> considering distances voxels's neighborhood.
%                
% -------------------------------------------------------------------------
% OUTPUTS:
% 
% directions.eigenvectors : Nx9 numeric. 1:3 is the first autovector of the
%                          covariance matrix of the voxel N and its
%                          neighborhood. 4:6 the second and 7:9 the third.
%
%
% directions.eigenvalues : Nx3 numeric. 1 is the first eigenvalues of the
%                         covariance matrix of the voxel N and its
%                         neighborhood. 2 the second and 3 the third.
%
% -------------------------------------------------------------------------
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 23/12/2020

if mode == 1
%% With distance

directions.eigenvectors = zeros(size(vx.Location,1),9);
directions.eigenvalues = zeros(size(vx.Location,1),3);

idx = rangesearch(vx.Location,vx.Location,maxDistance);

% figure; pcshow(vx.Location);

for i = 1:length(vx.Location)
%     hold on; pcshow(vx.Location(idx{i},:),'g','MarkerSize', 100);

    % Cálculo de autovectores y autovalores en la matriz de covarianza. 
    % Es más rapido que hacer pca()
    [coeff,latent] = pcacov(cov(vx.Location(idx{i},:)));
%     [coeff,~,latent] = pca(vx.Location(idx{i},:));
    
    if ~isempty(latent)
        
        directions.eigenvectors(i,1:3) = coeff(:,1);
        directions.eigenvalues(i,1)    = latent(1) /sum(latent);

        if numel(latent) >= 2
            
            directions.eigenvectors(i,4:6) = coeff(:,2);
            directions.eigenvalues(i,2)    = latent(2) /sum(latent);
        
            if numel(latent) == 3
                
                directions.eigenvectors(i,7:9) = coeff(:,3);
                directions.eigenvalues(i,3)    = latent(3) /sum(latent);
        
            end
        
        end
        
    end
   
end

elseif mode == 2
%% With neighbours and distance

firstEigenvectors = zeros(length(vx.vx.Location),3);
firstEigenvalues = zeros(length(vx.vx.Location),1);

secondEigenvectors = zeros(length(vx.vx.Location),3);
secondEigenvalues = zeros(length(vx.vx.Location),1);

thirdEigenvectors = zeros(length(vx.vx.Location),3);
thirdEigenvalues = zeros(length(vx.vx.Location),1);

idx = false(length(vx.neighbours_rows),1);

% figure; pcshow(VX.vx.Location);

for i = 1:length(vx.neighbours_rows)
    distance = pdist2(vx.Location(:,:), vx.Location(i,:)); % Distance to the voxel
    
    % Cleaning the cluster
    idx(:) = false;
    
    % The voxel is part of the cluster
    idx(i)= true;
    
    
    for j = 1:ceil((maxDistance / vx.grid))
        neighbors_propagation = vx.neighbours_rows(idx,:);
        neighbors_propagation = cat(1,neighbors_propagation);
        neighbors_propagation = neighbors_propagation(~isnan(neighbors_propagation));
        neighbors_propagation = unique(neighbors_propagation); % neighbours of the new voxels
        
        idx(neighbors_propagation) = true;
    end
    
    idx(distance > maxDistance) = false;

%     hold on; pcshow(VX.vx.Location(idx,:),'g','MarkerSize', 1000);

    
    [coeff,~,latent] = pca(vx.Location(idx,:));
    
    if ~isempty(latent)
        firstEigenvectors(i,:) = coeff(:,1);
        firstEigenvalues(i) = latent(1) /sum(latent);

        if numel(latent) >= 2
            secondEigenvectors(i,:) = coeff(:,2);
            secondEigenvalues(i) = latent(2) /sum(latent);
        
            if numel(latent) == 3
                thirdEigenvectors(i,:) = coeff(:,3);
                thirdEigenvalues(i) = latent(3) /sum(latent);
            end
        
        end
        
    end
   
end

directions.first.vector = firstEigenvectors;
directions.first.value = firstEigenvalues;
    
directions.second.vector = secondEigenvectors;
directions.second.value = secondEigenvalues;

directions.third.vector = thirdEigenvectors;
directions.third.value = thirdEigenvalues;
    
end
end

