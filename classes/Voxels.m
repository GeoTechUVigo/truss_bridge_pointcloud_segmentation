%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
Copyright (C) 2023 Daniel Lamas Novoa daniel.lamas.novoa@uvigo.gal
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

classdef Voxels
    %Añadir las propiedades que se especifican en el ejercicio. No es
    %necesario definir el tipo de dato de la propiedad - puede hacerse si
    %se considera oportuno -. Para definir una propiedad sin tipo de dato,
    %simplemente se debe escribir su nombre, por ejemplo:
    %properties
    %  prop1; 
    %  prop2; 
    %end
    
    properties 
        
        parent_cloud; % Matriz Nx3 con los N puntos de la nube a partir de la cual se voxeliza (copia de la nube de puntos original).
        parent_idx; % Cell array Mx1 donde M es el número total de vóxeles ocupados en la estructura. El elemento parent_idx{i} contendrá los índices de los puntos en el vóxel de índices índices(i).
        index; % Vector Mx1 con los índices de los vóxeles que contengan al menos un punto de la nube.
        Location; % Matriz Mx3 con las coordenadas de los centroides de cada vóxel. (El elemento Location(i,:) contendrá las coordenadas del vóxel de índice índices(i), en el cual se encuentran los puntos de índice parent_idx{i} para la nube en parent_cloud.
        mean_z; % Vector Mx1 con el valor medio de las alturas de los puntos en cada vóxel. mean_z(i) indica la media de las alturas de los puntos en el vóxel de índice índices(i).
        variance_z, % Vector Mx1 con la varianza de las alturas de los puntos en cada vóxel. variance_z(i) indica la varianza de las alturas de los puntos en el vóxel de índice índices(i).
        range_z;
        max_z;
        min_z;
        neighbours; % Cell array Mx1 con los índices de los vóxeles vecinos que contengan puntos respecto a cada vóxel. neighbors(i) contendrá los índices de los vóxeles vecinos al vóxel de índice índices(i).
        neighbours_rows; % Igual que neighbours pero en vez de indicar los índices de los vecinos, indica las filas en .Location de los vecinos
        grid_size; % Valor utilizado como tamaño de rejilla (input para el constructor de clase)
        xlimit; % Número de vóxeles en la dirección X
        ylimit; % Número de vóxeles en la dirección Y
        zlimit; % Número de vóxeles en la dirección Z
        intensity;
        timeStamp;
        Color;
        angle;
        Normal;
        grid;
       
    end
    %Constructor de la clase. Es el único método que es necesario
    %definir para este ejercicio. En él debe realizarse la
    %rasterización, y deben crearse las imágenes que se piden en el
    %enunciado. 
    methods
        %Constructor
        function obj = Voxels(cloud, grid)
        %Inputs: 
        %- cloud deberá ser un objeto de tipo pointCloud_ (no es
        %necesario realizar comprobaciones. Puede hacerse si se
        %considera oportuno.)
        % - grid será un número (double) que defina el tamaño
        % de celda. 
        % Recordar que, para definir una propiedad, y dado que nuestro
        % output se llama 'obj', nos referiremos a ella como
        % obj.nombrePropiedad. Por ejemplo, podemos definir
        % inmediatamente la propiedad parent_cloud:
            obj.parent_cloud = cloud.Location;
         %%%% Generación estructura ráster
        % Cálculo de mínimos y máximos de las coordenadas x,y,z de la
        % nube de puntos (~6 líneas)
            min_x = min(cloud.Location(:,1));
            max_x = max(cloud.Location(:,1));
            min_y = min(cloud.Location(:,2));
            max_y = max(cloud.Location(:,2));
            min_z = min(cloud.Location(:,3));
            max_z = max(cloud.Location(:,3));

        
        % Cálculo del número de vóxeles en x,y,z (~3 lineas)
            n_x = ceil((max_x - min_x) / grid);
            n_y = ceil((max_y - min_y) / grid);
            n_z = ceil((max_z - min_z) / grid);
        %Cálculo del tamaño real de vóxel tras ajuste a los puntos. (~3
        %líneas).

            paso_x =(max_x - min_x) / n_x;
            paso_y = (max_y - min_y) / n_y;
            paso_z =(max_z - min_z) / n_z;
            
        %Cálculo de las coordenadas x e y de cada punto (~ 2 líneas)
            dim_x = uint64(ceil((obj.parent_cloud(:,1) - min_x) / paso_x));
            dim_y = uint64(ceil((obj.parent_cloud(:,2) - min_y) / paso_y));
            dim_z = uint64(ceil((obj.parent_cloud(:,3) - min_z) / paso_z));

        %Mover el caso dim_x = 0 y dim_y = 0 a la celda de coordenada 1
            dim_x(dim_x==0) = 1;
            dim_y(dim_y==0) = 1;
            dim_z(dim_z==0) = 1;


        %Definición de un índice líneal para cada triplete de coordenadas.
        %(~1 línea)
        idx = dim_x + (dim_y - 1) * uint64(n_x) + (dim_z - 1) * uint64(n_x) *  uint64(n_y);
        
        %%%%% Estructuración de la información: Búsqueda, de la forma más
        %eficiente posible, del número de vóxeles ocupados y del número de
        %puntos en cada vóxel. 
        
        %Crear una variable de nombre data, Nx4, que contenga los puntos de
        %la nube y su posición de fila (vector columna [1,2,3...N].(1-2 líneas). 
        
        data = [obj.parent_cloud, [1:length(obj.parent_cloud)]'];
        %Ordenar el vector de índices de menor a mayor, obteniendo los
        %índices que describan el cambio de filas que sufre el vector 
        %original (~ 1 línea. Función sortrows). 
        [idx, shift_index] = sortrows(idx);
        
        % Utilizar la función unique para obtener un vector que contenga
        % unicamente los vóxeles ocupados, así como un vector de índices
        % correspondientes al primer vóxel de cada valor único. (~ 1 linea)
        
        [occ_voxels, first_idx] = unique(idx);
        
        % Obtener el número de puntos en cada vóxel utilizando la
        % información en first_idx. El número de puntos con indice
        % idx(1) será la diferencia entre first_idx(2) y first_idx(1).
        % Crear una segunda columna en occ_pixels de tal forma que
        % occ_pixels(:,1) tenga los índices ocupados, y occ_pixels(:,2)
        % el número de puntos con dicho índice. (~ 2 lineas).
        n_points = [first_idx(2:length(first_idx)) - first_idx(1:length(first_idx)-1) ; length(idx) + 1 - first_idx(length(first_idx))];
        occ_voxels=[occ_voxels, n_points];
        
        %%%% Cálculo de características de cada vóxel.
        
        %Con los datos que tenemos, podemos definir el número de vóxeles
        %ocupados (~ 1 línea) y el cell array que contendrá los índices de
        %los puntos en cada voxel (~ 1 linea): 
        nVoxels = length(occ_voxels);
        parent_idx = cell(nVoxels,1);
        %Y definir también el tamaño de los vectores donde almacenaremos
        %las propiedades (utilizar función zeros) (~ 3 lineas):
        centroids   = zeros(nVoxels,3);
        mean_z      = zeros(nVoxels,1);
        variance_z  = zeros(nVoxels,1);
        range_z     = zeros(nVoxels,1);
        max_z       = zeros(nVoxels,1);
        min_z       = zeros(nVoxels,1);
        
        %Juntar los puntos de 'data' y sus índices idx en una misma
        %variable de nombre 'data_sort'. Para eso, reordenarlos según el vector
        %shift_index para que los puntos se muestren de forma ordenada en
        %función del vóxel al que pertenecen.(~ 1 línea)
        
        data_sort = data(shift_index,:);
        
        has_intensity = false;
        has_timeStamp = false;
        has_Color = false;
        has_angle = false;
        has_Normal = false;
        
        if ~isempty(cloud.intensity)
            intensity_sort = [cloud.intensity(shift_index,:),idx];
            intensity   = zeros(nVoxels,1);
            has_intensity = true;
        end
        
        if ~isempty(cloud.timeStamp)
            timeStamp_sort = [cloud.timeStamp(shift_index,:),idx];
            timeStamp   = zeros(nVoxels,1);
            has_timeStamp = true;
        end
        
        if ~isempty(cloud.Color)
            Color_sort = [cloud.Color(shift_index,:),idx];
            Color       = zeros(nVoxels,3);
            has_Color = true;
        end
        
        if ~isempty(cloud.angle)
            angle_sort = [cloud.angle(shift_index,:),idx];
            angle       = zeros(nVoxels,1);
            has_angle = true;
        end
        
        if ~isempty(cloud.Normal)
            Normal_sort = [cloud.Normal(shift_index,:),idx];
            Normal      = zeros(nVoxels,3);
            has_Normal = true;
        end
        
        
        
        %% Cálculo de vóxeles vecinos. Conociendo la estructura de
        %voxelización y la forma de asignar índices, obtener los índices de
        %los vóxeles vecinos a cada vóxel ocupado.
%         neighbors = cell(nVoxels,1);
        
        %Para el voxel (dim_x, dim_y, dim_z) de indice idx, le
        %corresponderán 26 posibles vecinos (9 abajo, 9 arriba, 8 en su
        %plano). Calcular en neighbors{i} los índices vecinos a idx(i)
        %utilizando el propio idx o sus correspondientes coordenadas
        %(dim_x, dim_y, dim_z). 
        
        
        neighbor_offset =[(- n_x*n_y - n_x - 1), (- n_x*n_y - n_x), (- n_x*n_y - n_x + 1),...
            (- n_x*n_y - 1), (- n_x*n_y), (- n_x*n_y + 1),...
            (- n_x*n_y + n_x -1), (- n_x*n_y + n_x), (- n_x*n_y + n_x + 1),...
            (- n_x - 1), (-n_x), (-n_x + 1),...
            (- 1), (+ 1),...
            (+ n_x - 1), (+ n_x), (+ n_x +1),...
            (+ n_x*n_y - n_x - 1), (+ n_x*n_y - n_x), (+ n_x*n_y - n_x + 1),...
            (+ n_x*n_y - 1), (+ n_x*n_y), (+ n_x*n_y + 1),...
            (+ n_x*n_y + n_x - 1), (+ n_x*n_y + n_x), (+ n_x*n_y + n_x + 1)];

        idx_neighbor = double(occ_voxels(:,1)) + neighbor_offset; % matriz con los posibles vecinos de cada voxel

        % Pongo como NaN los que no son vecinos
        idx_neighbor(rem(occ_voxels(:,1), n_x) == 0,[3 6 9 12 14 17 20 23 26]) = NaN; % Límite superior X
        idx_neighbor(rem(occ_voxels(:,1), n_x) == 1, [1 4 7 10 13 15 18 21 24]) = NaN; % Límte inferior X
        idx_neighbor(rem(ceil(occ_voxels(:,1) / n_x), n_y) == 0,[7 8 9 15 16 17 24 25 26]) = NaN; % Límite superior Y
        idx_neighbor(rem(ceil(occ_voxels(:,1) / n_x), n_y) == 1, [1 2 3 10 11 12 18 19 20]) = NaN; % Límte inferior X
        idx_neighbor(idx_neighbor > n_x * n_y * n_z) = NaN; % Límite superior Z
        idx_neighbor(idx_neighbor < 0) = NaN; % Límite inferior Z

        %%
        % Elimino los vecinos que no hay porque el voxel no tiene puntos de
        % la nube
        [a,rows_neighbor]=ismember(idx_neighbor,occ_voxels(:,1));
        idx_neighbor(~a) = NaN;
        rows_neighbor(~a) = NaN;

        % Relleno la lista de celdas
%         for i = 1: nVoxels
%             neighbors{i} = idx_neighbor_(i,:);
%         end
%       XY_not_ground_neighbors = VX.neighbours(not_ground(1:3),10:17);
     
        %% Programar un bucle for para calcular las propiedades. Puede ser
        % necesario definir variables auxiliares antes de comenzar el
        % bucle. Esquemáticamente: 
        % 1. Recorrer todos los vóxeles ocupados.
        % 2. Conocer el número de puntos en cada vóxel.
        % 3. Extraer la información (en 'data_sort') de forma
        % secuencial (si en el primer vóxel hay 10 puntos y en el
        % segundo 12, las diez primeras filas de data componen el
        % vóxel idx(1)  y las filas 11-22 el voxel idx(2).. y así sucesivamente)
        % 4. Para la primera iteración, estaremos rellenando la
        % información del voxel en centroids(i,:), mean_z(i)...
        % 5. Para rellenar parent_idx{i}, recordar que los índices de
        % cada punto deberían estar en la última columna de data_sort.
        % (~15-20 líneas).
        
            lastPoint = 1;
            for i = 1:nVoxels
                
                nextPoint = lastPoint + occ_voxels(i,2);
                points = lastPoint : nextPoint - 1;
                
                centroids(i,:)  = [mean(data_sort(points , 1),1),mean(data_sort(points , 2),1),mean(data_sort(points , 3),1)];
                parent_idx{i}   = data_sort(points,4);
                mean_z(i)       = mean(data_sort(points,3),1);
                variance_z(i)   = var(data_sort(points,3),0,1);
                range_z(i)      = range(data_sort(points,3),1);
                max_z(i)        = max(data_sort(points,3),[],1);
                min_z(i)        = min(data_sort(points,3),[],1);
                                
                if has_intensity
                    intensity(i)    = mean(intensity_sort(points,1),1);
                end
                if has_timeStamp
                    timeStamp(i)    = mean(timeStamp_sort(points,1),1);
                end
                if has_Color
                    Color(i,:)      = [mean(Color_sort(points , 1),1),mean(Color_sort(points , 2),1),mean(Color_sort(points , 3),1)];
                end
                if has_angle
                    angle(i)        = mean(angle_sort(points,1),1);
                end
                if has_Normal
                   Normal(i,:)     = [mean(Normal_sort(points , 1),1),mean(Normal_sort(points , 2),1),mean(Normal_sort(points , 3),1)];
                end
                
                lastPoint = nextPoint;
            end
           
            
     
        %% Asignar los resultados a la salida del constructor (obj).
            obj.Location    = centroids;
            obj.parent_idx  = parent_idx;
            obj.mean_z = mean_z;
            obj.variance_z = variance_z;
            obj.range_z = range_z;
            obj.max_z = max_z;
            obj.min_z = min_z;
            obj.neighbours = idx_neighbor;
            obj.neighbours_rows = rows_neighbor;
            obj.index = occ_voxels(:,1);
            obj.grid_size = n_x*n_y*n_z;
            obj.xlimit = n_x;
            obj.ylimit = n_y;
            obj.zlimit = n_z;
            obj.grid = grid;
            
            if ~isempty(cloud.intensity)
            obj.intensity = intensity;
            end
            if ~isempty(cloud.timeStamp)
            obj.timeStamp = timeStamp;
            end
            if ~isempty(cloud.Color)
            obj.Color = Color;
            end
            if ~isempty(cloud.angle)
            obj.angle = angle;
            end
            if ~isempty(cloud.Normal)
            obj.Normal = Normal;
            end
        end
        
        function voxelOut = select(obj, varargin)
            %select Select points specified by index.
            %
            %  voxelOut = select(voxel, indices) returns a Voxelized
            %  object that contains the points selected using linear
            %  indices.
            %
            %  voxelOut = select(voxel, row, column) returns a Voxelized
            %  object that contains the points selected using row and
            %  column subscripts. This syntax applies only to organized
            %  point cloud (M-by-N-by-3).

            narginchk(2, 3);
            
            if nargin == 3
                if (strcmp(varargin{2},'index'))
                    rows = ismember(obj.index,varargin{1});
                end

            elseif nargin == 2
                rows = varargin{1};
%                  validateattributes(rows, {'numeric'}, ...
%                      {'real','nonsparse', 'vector', 'integer'}, mfilename, 'index');
            else
                % Subscript syntax is only for organized point cloud
                if ndims(this.Location) ~= 3
                    error(message('vision:pointcloud:organizedPtCloudOnly'));
                end
                
                row = varargin{1};
                column = varargin{2};
                validateattributes(row, {'numeric'}, ...
                    {'real','nonsparse', 'vector', 'integer'}, mfilename, 'row');
                validateattributes(column, {'numeric'}, ...
                    {'real','nonsparse', 'vector', 'integer'}, mfilename, 'column');
                rows = sub2ind([size(this.Location,1), size(this.Location,2)], row, column);
            end
            
            
            
            obj.Location = obj.Location(rows,:);
            obj.index = obj.index(rows);
            obj.parent_idx  = obj.parent_idx(rows);
            obj.mean_z = obj.mean_z(rows);
            obj.variance_z = obj.variance_z(rows);
            obj.range_z = obj.range_z(rows);
            obj.max_z = obj.max_z(rows);
            obj.min_z = obj.min_z(rows);
            obj.grid = obj.grid;
            
            
            if ~isempty(obj.intensity)
            obj.intensity = obj.intensity(rows);
            end
            if ~isempty(obj.timeStamp)
            obj.timeStamp = obj.timeStamp(rows);
            end
            if ~isempty(obj.Color)
            obj.Color = obj.Color(rows);
            end
            if ~isempty(obj.angle)
            obj.angle = obj.angle(rows);
            end
            if ~isempty(obj.Normal)
            obj.Normal = obj.Normal(idx);
            end
            
            new_neighbours = obj.neighbours(rows,:);
            [a,rows_neighbor]=ismember(new_neighbours,obj.index);
            new_neighbours(~a) = NaN;
            rows_neighbor(~a) = NaN; 
            

            obj.neighbours = new_neighbours;            
            obj.neighbours_rows = rows_neighbor;
            
            voxelOut = obj;
        end        
    end
end
