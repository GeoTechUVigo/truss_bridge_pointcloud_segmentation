
function [] = SaveLas(parent_idx, components, path_in, path_out)


cloud = LASread(path_in);

%% Initialising
cloud.record.classification(:) = 0;
cloud.record.point_source_id(:) = 0;
cloud.record.user_data(:) = 0;
group = 1;

%% Classification
% Vertical face Lateral -> 1
% Vertical face Vertical -> 2
% Chord -> 3
% Horizontal face Lateral -> 4
% Horizontal face Vertical -> 5
% Inner face Lateral -> 6
% Inner face Horizontal -> 7

%% Vertical Faces
for i = 1:numel(components.verticalFaces)
    % Lateral Diagonals
    for j = 1:numel(components.verticalFaces{i}.lateralDiagonals)
        for k =1:numel(components.verticalFaces{i}.lateralDiagonals{j})
            % Indexes
            element = parent_idx(components.verticalFaces{i}.lateralDiagonals{j}{k});
            element = cat(1,element{:});

            % Classification
            cloud.record.classification(element)= 1;

            % Group
            cloud.record.point_source_id(element)= group;
            group = group + 1;
        end
    end

    % Vertical members
    for j = 1:numel(components.verticalFaces{i}.verticals)
        % Indexes
        element = parent_idx(components.verticalFaces{i}.verticals{j});
        element = cat(1,element{:});

        % Classification
        cloud.record.classification(element)= 2;
        % Group
        cloud.record.point_source_id(element)= group;
        group = group + 1;
    end
    % Chords
    for j = 1:numel(components.verticalFaces{i}.chords)
        % Indexes
        element = parent_idx(components.verticalFaces{i}.chords{j});
        element = cat(1,element{:});

        % Classification
        cloud.record.classification(element)= 3;
        % Group
        cloud.record.point_source_id(element)= group;
        group = group + 1;
    end
end

%% Horizontal sections
for i = 1:numel(components.horizontalFaces)
    % Lateral Diagonals
    for j = 1:numel(components.horizontalFaces{i}.lateralDiagonals)
        for k = 1:numel(components.horizontalFaces{i}.lateralDiagonals{j})
            % Indexes
            element = parent_idx(components.horizontalFaces{i}.lateralDiagonals{j}{k});
            element = cat(1,element{:});

            % Classification
            cloud.record.classification(element)= 4;

            % Group
            cloud.record.point_source_id(element)= group;
            group = group + 1;
        end
    end
    % Vertical members
    for j = 1:numel(components.horizontalFaces{i}.verticals)
        % Indexes
        element = parent_idx(components.horizontalFaces{i}.verticals{j});
        element = cat(1,element{:});
    
        % Classification
        cloud.record.classification(element)= 5;
        % Group
        cloud.record.point_source_id(element)= group;
        group = group + 1;
    end
end
%% Inner Faces
for i = 1:numel(components.innerFaces)
    % Lateral Diagonals
    for j = 1:numel(components.innerFaces{i}.lateralDiagonals)
        for k = 1:numel(components.innerFaces{i}.lateralDiagonals{j})
            % Indexes
            element = parent_idx(components.innerFaces{i}.lateralDiagonals{j}{k});
            element = cat(1,element{:});

            % Classification
            cloud.record.classification(element)= 6;

            % Group
            cloud.record.point_source_id(element)= group;
            group = group + 1;
        end
    end
    % Horizontal member in the middle
    for j = 1:numel(components.innerFaces{i}.horizontal)
        % Indexes
        element = parent_idx(components.innerFaces{i}.horizontal{j});
        element = cat(1,element{:});
    
        % Classification
        cloud.record.classification(element)= 7;
        % Group
        cloud.record.point_source_id(element)= group;
        group = group + 1;
    end
end

%% Save
LASwrite(cloud,char(path_out));

