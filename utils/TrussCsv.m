function [] = TrussCsv(components, vx, path, innerRelVertical, decimals)
% Function to save the coordinates of each element in a .csv
%
% The headers of each element is : TrussTag_face_type_id_dim
% In the case of inner faces also the vertical member that connect them:
% TrussTag_face_type_id_dim_idVerticalMembers
%
% First, The headers are created.
% Then, a NaN matrix is created. The number of columns is the number of
% headers. The number of row is the number of points of the element with
% hte most number of points. Once the matrix is generated, all the elements
% are looped in the same order as generating the header, writing its
% coordenates in the matrix.
% Lastly, a table with the matrix and the headers is saved as .csv
%
%--------------------------------------------------------------------------
% INPUTS :
% components : struct with the indexes of each element
%
% vx : point cloud
%
% path : string. Path to save the .csv
%
% innerRelVertical: matrix Nº inner faces x Nº vertical faces. Contains the
%                  index of the vertical members of each vertical face
%                  conected to each inner face.
%
% decimals: number of decimals.
%--------------------------------------------------------------------------
% OUTPUTS :
% A .csv is save in path with the coordinates in vx.Location of each
% element in components
%
%--------------------------------------------------------------------------
% Andrés Justo Domínguez.
% Daniel Lamas Novoa.
% Enxeñaría dos materiais, mecánica aplicada e construción.
% Escola de enxeñería industrial
% Grupo de xeotecnoloxía aplicada.
% Universidade de Vigo.
% 21/04/21

%% Count the elements to generate the output matrix knowing its dimensions
headers = string([]);
dimension = 0;

cont = 0;

for i = 1:length(components(:,1))
    TrussTag = "TRUSS" + string(i); % Truss
    
    %% Face VerticalFaces
    for j = 1:length(components.verticalFaces)
        if j == 1
            face = "R";
        else
            face = "L";
        end
        
        %% Lateral diagonals
        for k = 1:length(components.verticalFaces{j}.lateralDiagonals)
            if k == 1
                type = "LU"; % Up
            else
                type = "LD"; % Down
            end
            for q = 1:length(components.verticalFaces{j}.lateralDiagonals{k})
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z';
                
                dimension= max(dimension,length(vx.Location(components.verticalFaces{j}.lateralDiagonals{k}{q},1)));
            end       
        end
        
        %% Verticals
        type = "VM";
        for q = 1:length(components.verticalFaces{j}.verticals)
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z';
                
                dimension= max(dimension,length(vx.Location(components.verticalFaces{j}.verticals{q},1)));
                
                % Rewrite id of this member
                innerRelVertical(q, j) = cont;
        end
        
        %% Chords
        type = "CH";
        for q = 1:length(components.verticalFaces{j}.chords)
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z';
                
                dimension= max(dimension,length(vx.Location(components.verticalFaces{j}.chords{q},1)));
        end
    end 
    %% Horizontals
     for j = 1:length(components.horizontalFaces)
         if j == 1
             face = "D";
         else
             face = "U";
         end
         
        %% Horizontal diagonals
        for k = 1:length(components.horizontalFaces{j}.lateralDiagonals)
            if k == 1
                type = "HU"; % Up
            else
                type = "HD"; % Down
            end
            for q = 1:length(components.horizontalFaces{j}.lateralDiagonals{k})
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z';
                
                dimension= max(dimension,length(vx.Location(components.horizontalFaces{j}.lateralDiagonals{k}{q},1)));
            end       
        end
       
        %% Horizontal members      
        type = "HM";
        for q = 1:length(components.horizontalFaces{j}.verticals)
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y';
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z';
                
                dimension= max(dimension,length(vx.Location(components.horizontalFaces{j}.verticals{q},1)));
        end
     end
    %% San Andres
    for j = 1:length(components.innerFaces)
        
        face = "I" + string(j);

        %% Inner face
        for k = 1:length(components.innerFaces{j}.lateralDiagonals)
            if k == 1
                type = "IU"; % Up
            else
                type = "ID"; % Down
            end
            for q = 1:length(components.innerFaces{j}.lateralDiagonals{k})
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                
                dimension= max(dimension,length(vx.Location(components.innerFaces{j}.lateralDiagonals{k}{q},1)));
            end       
        end
       
        %% San Andres horizontal     
        type = "IH";
        for q = 1:length(components.innerFaces{j}.horizontal)
                cont = cont+1;
                id = type + string(cont); % number of the element
                
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'X' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Y' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                headers(numel(headers)+1) = TrussTag + '_' + face + '_' + type + '_' + id + '_' + 'Z' + '_' + 'VM' + string(innerRelVertical(j,1)) + '_' + 'VM' + string(innerRelVertical(j,2));
                
                dimension= max(dimension,length(vx.Location(components.innerFaces{j}.horizontal{q},1)));
        end
     end
end

%% Fill the output matrix
output = NaN(dimension,numel(headers));
col = 1;
for i = 1:length(components(:,1))
    
    %% Face VerticalFaces
    for j = 1:length(components.verticalFaces)
        
        %% Lateral diagonals
        for k = 1:length(components.verticalFaces{j}.lateralDiagonals)
            for q = 1:length(components.verticalFaces{j}.lateralDiagonals{k})               
                output(1:length(components.verticalFaces{j}.lateralDiagonals{k}{q}),col:col+2) = vx.Location(components.verticalFaces{j}.lateralDiagonals{k}{q},:);
                col = col+3;
            end       
        end
        
        %% Verticals
        for q = 1:length(components.verticalFaces{j}.verticals)
            output(1:length(components.verticalFaces{j}.verticals{q}),col:col+2) = vx.Location(components.verticalFaces{j}.verticals{q},:);
            col = col+3;
        end
        
        %% Chords
        for q = 1:length(components.verticalFaces{j}.chords)
            output(1:length(components.verticalFaces{j}.chords{q}),col:col+2) = vx.Location(components.verticalFaces{j}.chords{q},:);
            col = col+3;
        end
    end 
    %% Horizontals
     for j = 1:length(components.horizontalFaces)
         
        %% Horizontal diagonals
        for k = 1:length(components.horizontalFaces{j}.lateralDiagonals)
            for q = 1:length(components.horizontalFaces{j}.lateralDiagonals{k})
                output(1:length(components.horizontalFaces{j}.lateralDiagonals{k}{q}),col:col+2) = vx.Location(components.horizontalFaces{j}.lateralDiagonals{k}{q},:);
                col = col+3;
            end       
        end
       
        %% Horizontal members      
        for q = 1:length(components.horizontalFaces{j}.verticals)
            output(1:length(components.horizontalFaces{j}.verticals{q}),col:col+2) = vx.Location(components.horizontalFaces{j}.verticals{q},:);
            col = col+3;
        end
     end
    %% San Andres
     for j = 1:length(components.innerFaces)
         
        %% San Andres cross
        for k = 1:length(components.innerFaces{j}.lateralDiagonals)
            for q = 1:length(components.innerFaces{j}.lateralDiagonals{k})
                output(1:length(components.innerFaces{j}.lateralDiagonals{k}{q}),col:col+2) = vx.Location(components.innerFaces{j}.lateralDiagonals{k}{q},:);
                col = col+3;
            end       
        end
       
        %% San Andres horizontal     
        for q = 1:length(components.innerFaces{j}.horizontal)
            output(1:length(components.innerFaces{j}.horizontal{q}),col:col+2) = vx.Location(components.innerFaces{j}.horizontal{q},:);
            col = col+3;
        end
    end
end

%% Save the matrix with the headres as .csv
output = round(output, decimals);
T = array2table(output);
T.Properties.VariableNames = headers;
writetable(T,path, 'Delimiter', ',');
