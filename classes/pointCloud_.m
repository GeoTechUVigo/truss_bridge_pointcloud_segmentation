%{
Copyright (C) 2023 GeoTECH Group geotech@uvigo.gal
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


classdef pointCloud_ < matlab.mixin.Copyable & vision.internal.EnforceScalarHandle
%Clase modificada a partir de la clase pointCloud_ de matlab. Se añaden las
%propiedades Intensity, timeStamp, angle, y sensor al objeto para representar mejor
%la nube de puntos. 
    
% pointCloud Object for storing a 3-D point cloud.
%   ptCloud = pointCloud(xyzPoints) creates a point cloud object whose
%   coordinates are specified by an M-by-3 or M-by-N-by-3 matrix xyzPoints.
%
%   ptCloud = pointCloud(xyzPoints, Name, Value) specifies additional
%   name-value pair arguments described below:
%
%   'Color'             A uint8 matrix C specifying the color of each point
%                       in the point cloud. If the coordinates are an
%                       M-by-3 matrix, C must be an M-by-3 matrix, where
%                       each row contains a corresponding RGB value. If the
%                       coordinates are an M-by-N-by-3 matrix xyzPoints, C
%                       must be an M-by-N-by-3 matrix containing
%                       1-by-1-by-3 RGB values for each point.
%
%                       Default: []
%
%   'Normal'            A matrix NV specifying the normal vector at
%                       each point. If the coordinates are an M-by-3
%                       matrix, NV must be an M-by-3 matrix, where each row
%                       contains a corresponding normal vector. If the
%                       coordinates are an M-by-N-by-3 matrix xyzPoints, NV
%                       must be an M-by-N-by-3 matrix containing
%                       1-by-1-by-3 normal vector for each point.
%
%                       Default: []
%
%   pointCloud properties :
%   Location        - Matrix of [X, Y, Z] point coordinates (read only)
%   Color           - Matrix of point RGB colors
%   Intensity       - Vector of point intensities
%   timeStamp       - Vector of time stamps
%   angle           - Vector of angles
%   sensor          - Vector with the index of the sensor´s point
%   Normal          - Matrix of [NX, NY, NZ] point normal directions 
%   Count           - Number of points (read only)
%   XLimits         - The range of coordinates along X-axis (read only)
%   YLimits         - The range of coordinates along Y-axis (read only)
%   ZLimits         - The range of coordinates along Z-axis (read only)
%
%   pointCloud methods:
%   findNearestNeighbors   - Retrieve K nearest neighbors of a point
%   findNeighborsInRadius  - Retrieve neighbors of a point within a specified radius
%   findPointsInROI        - Retrieve points within a region of interest
%   select                 - Select points specified by index
%   removeInvalidPoints    - Remove invalid points
%
%   Notes 
%   ----- 
%   - To reduce memory use, point locations are suggested to be stored as single.
%
%   Example : Find the shortest distance between two point clouds
%   -------------------------------------------------------------
%   ptCloud1 = pointCloud(rand(100, 3, 'single'));
%        
%   ptCloud2 = pointCloud(1+rand(100, 3, 'single'));
%
%   minDist = inf;
%
%   % Find the nearest neighbor in ptCloud2 for each point in ptCloud1
%   for i = 1 : ptCloud1.Count
%       point = ptCloud1.Location(i,:);
%       [~, dist] = findNearestNeighbors(ptCloud2, point, 1);
%       if dist < minDist 
%           minDist = dist;
%       end
%   end
%
% See also pcshow, pcfromkinect, reconstructScene, pcread, pcwrite, 
%          pctransform, pcdenoise, pcdownsample, pcregrigid, pcmerge
 
%  Copyright 2014 The MathWorks, Inc.
%
% References
% ----------
% Marius Muja and David G. Lowe, "Fast Approximate Nearest Neighbors with
% Automatic Algorithm Configuration", in International Conference on
% Computer Vision Theory and Applications (VISAPP'09), 2009

    properties (GetAccess = public, SetAccess = public)
        % Location is an M-by-3 or M-by-N-by-3 matrix. Each entry specifies
        % the x, y, z coordinates of a point.
        Location = single([]);
        intensity = double([]);
        timeStamp = double([]);
        sensor = uint8([]);
    end
    
    properties (Access = public)
        % Color is an M-by-3 or M-by-N-by-3 uint8 matrix. Each entry
        % specifies the RGB color of a point.
        Color = uint8([]);

        angle = double([]);
        % Normal is an M-by-3 or M-by-N-by-3 matrix. Each entry
        % specifies the x, y, z component of a normal vector.
        Normal = single([]);                
    end
    
    properties(Access = protected, Transient)
        Kdtree = [];
    end
    
    properties(Access = protected, Hidden)
        IsOrganized;
        KdtreeIndexState
    end
       
    properties(Dependent)
        % Count specifies the number of points in the point cloud.
        Count;
        % XLimits is a 1-by-2 vector that specifies the range of point
        % locations along X axis.
        XLimits;
        % YLimits is a 1-by-2 vector that specifies the range of point
        % locations along Y axis.
        YLimits;
        % ZLimits is a 1-by-2 vector that specifies the range of point
        % locations along Z axis.
        ZLimits;
    end
    
    methods
        %==================================================================
        % Constructor
        %==================================================================
        % pcd = pointCloud(xyzPoints);
        % pcd = pointCloud(..., 'Color', C);
        % pcd = pointCloud(..., 'Normal', nv);
        function this = pointCloud_(varargin)
            narginchk(1, 13);
            
            [xyzPoints, C, nv, intensity, timeStamp, angle, sensor] = validateAndParseInputs(varargin{:}); 
            
            this.Location = double(xyzPoints);
            this.Color = C;
            this.Normal = nv;
            this.intensity = intensity;
            this.timeStamp = timeStamp;
            this.angle = angle;
            this.sensor = sensor;
            
            this.IsOrganized = ~ismatrix(this.Location); % M-by-N-by-3
        end
        
        %==================================================================
        % K nearest neighbor search
        %==================================================================
        function [indices, dists] = findNearestNeighbors(this, point, K, varargin)
            %findNearestNeighbors Find the nearest neighbors of a point.
            %
            %  [indices, dists] = findNearestNeighbors(ptCloud, point, K)
            %  returns the K nearest neighbors of a query point. The point
            %  is an [X, Y, Z] vector. indices is a column vector that 
            %  contains K linear indices to the stored points in the point 
            %  cloud. dists is a column vector that contains K Euclidean 
            %  distances to the point. 
            %
            %  [...] = findNearestNeighbors(... , Name, Value)
            %  specifies additional name-value pairs described below:
            %  
            %  'Sort'          True or false. When set to true,
            %                  the returned indices are sorted in the
            %                  ascending order based on the distance
            %                  from a query point.
            %
            %                  Default: false
            %
            %  'MaxLeafChecks' An integer specifying the number of leaf
            %                  nodes that are searched in Kdtree. If the tree
            %                  is searched partially, the returned result may
            %                  not be exact. A larger number may increase the
            %                  search accuracy but reduce efficiency. When this
            %                  number is set to inf, the whole tree will be
            %                  searched to return exact search result. 
            %
            %                  Default: inf
            %
            %  Example : Find K-nearest neighbors
            %  ---------
            %  % Create a point cloud object with randomly generated points
            %  ptCloud = pointCloud(1000*rand(100,3,'single'));
            %
            %  % Define a query point and number of neighbors
            %  point = [50, 50, 50];
            %  K = 10;
            %
            %  % Get the indices and distances of 10 nearest points
            %  [indices, dists] = findNearestNeighbors(ptCloud, point, K);

            narginchk(3, 7);

            validateattributes(point, {'single', 'double'}, ...
                {'real', 'nonsparse', 'finite', 'size', [1, 3]}, mfilename, 'point');
            
            if isa(this.Location, 'single')
                point = single(point);
            else
                point = double(point);
            end

            validateattributes(K, {'single', 'double'}, ...
                {'nonsparse', 'scalar', 'positive', 'integer'}, mfilename, 'K');

            K = min(double(K), numel(this.Location)/3);

            [doSort, MaxLeafChecks] = validateAndParseSearchOption(varargin{:});

            % Use bruteforce if there are fewer than 500 points
            if numel(this.Location)/3 < 500
                if ~this.IsOrganized
                    allDists = visionSSDMetric(point', this.Location');
                else
                    allDists = visionSSDMetric(point', reshape(this.Location, [], 3)');
                end
                % This function will ensure returning actual number of neighbors
                % The result is already sorted
                [dists, indices] = vision.internal.partialSort(allDists, K);
                tf = isfinite(dists);
                indices = indices(tf);
                dists = dists(tf);
            else
                this.buildKdtree();

                searchOpts.checks = MaxLeafChecks;
                searchOpts.eps = 0;
                [indices, dists, valid] = this.Kdtree.knnSearch(point, K, searchOpts);
                
                % This step will ensure returning actual number of neighbors
                indices = indices(1:valid);
                dists = dists(1:valid);

                % Sort the result if specified
                if doSort
                    [dists, IND] = sort(dists);
                    indices = indices(IND);
                end
            end

            if nargout > 1
                dists = sqrt(dists);
            end
        end

        %==================================================================
        % Radius search
        %==================================================================
        function [indices, dists] = findNeighborsInRadius(this, point, radius, varargin)
            %findNeighborsInRadius Find the neighbors within a radius.
            %
            %  [indices, dists] = findNeighborsInRadius(ptCloud, point, radius) 
            %  returns the neighbors within a radius of a query point.
            %  The query point is an [X, Y, Z] vector. indices is a column
            %  vector that contains the linear indices to the stored points in
            %  the point cloud object ptCloud. dists is a column vector that
            %  contains the Euclidean distances to the query point.
            %
            %  [...] = findNeighborsInRadius(... , Name, Value) specifies
            %  specifies additional name-value pairs described below:
            %  
            %  'Sort'          True or false. When set to true,
            %                  the returned indices are sorted in the
            %                  ascending order based on the distance
            %                  from a query point.
            %
            %                  Default: false
            %
            %  'MaxLeafChecks' An integer specifying the number of leaf
            %                  nodes that are searched in Kdtree. If the tree
            %                  is searched partially, the returned result may
            %                  not be exact. A larger number may increase the
            %                  search accuracy but reduce efficiency. When this
            %                  number is set to inf, the whole tree will be
            %                  searched to return exact search result. 
            %
            %                  Default: inf
            %
            %   Example : Find neighbors within a given radius using Kdtree
            %   -----------------------------------------------------------
            %   % Create a point cloud object with randomly generated points
            %   ptCloud = pointCloud(100*rand(1000, 3, 'single'));
            %
            %   % Define a query point and search radius
            %   point = [50, 50, 50];
            %   radius =  5;
            %
            %   % Get all the points within a radius
            %   [indices, dists] = findNeighborsInRadius(ptCloud, point, radius);
        
            narginchk(3, 7);

            validateattributes(point, {'single', 'double'}, ...
                {'real', 'nonsparse', 'finite', 'size', [1, 3]}, mfilename, 'point');

            if isa(this.Location, 'single')
                point = single(point);
            else
                point = double(point);
            end

            validateattributes(radius, {'single', 'double'}, ...
                {'nonsparse', 'scalar', 'nonnegative', 'finite'}, mfilename, 'radius');
            
            radius = double(radius);
            
            [doSort, MaxLeafChecks] = validateAndParseSearchOption(varargin{:});
            
            % Use bruteforce if there are less than 500 points
            if numel(this.Location)/3 < 500
                if ~this.IsOrganized
                    allDists = visionSSDMetric(point', this.Location');
                else
                    allDists = visionSSDMetric(point', reshape(this.Location, [], 3)');
                end
                indices = uint32(find(allDists <= radius^2))';
                dists = allDists(indices)';
                tf = isfinite(dists);
                indices = indices(tf);
                dists = dists(tf);
            else
                this.buildKdtree();

                searchOpts.checks = MaxLeafChecks;
                searchOpts.eps = 0;
                [indices, dists] = this.Kdtree.radiusSearch(point, radius, searchOpts);
            end
            
            % Sort the result if specified
            if doSort
                [dists, IND] = sort(dists);
                indices = indices(IND);
            end
            
            if nargout > 1
                dists = sqrt(dists);
            end
        end       
        
        %==================================================================
        % Box search
        %==================================================================
        function indices = findPointsInROI(this, roi)
            %findPointsInROI Find points within a region of interest.
            %
            %  indices = findPointsInROI(ptCloud, roi) returns the points
            %  within a region of interest, roi. The roi is a cuboid
            %  specified as a 1-by-6 vector in the format of [xmin, xmax,
            %  ymin, ymax, zmin, zmax]. indices is a column vector that
            %  contains the linear indices to the stored points in the
            %  point cloud object, ptCloud.
            %
            %   Example : Find points within a given cuboid
            %   -------------------------------------------
            %   % Create a point cloud object with randomly generated points
            %   ptCloudA = pointCloud(100*rand(1000, 3, 'single'));
            %
            %   % Define a cuboid
            %   roi = [0, 50, 0, inf, 0, inf];
            %
            %   % Get all the points within the cuboid
            %   indices = findPointsInROI(ptCloudA, roi);
            %   ptCloudB = select(ptCloudA, indices);
            %
            %   pcshow(ptCloudA.Location, 'r');
            %   hold on;
            %   pcshow(ptCloudB.Location, 'g');
            %   hold off;
            
            narginchk(2, 2);

            validateattributes(roi, {'single', 'double'}, ...
                {'real', 'nonsparse', 'numel', 6}, mfilename, 'roi');
            
            if isvector(roi)
                roi = reshape(roi, [2, 3])';
            end
            
            if any(roi(:, 1) > roi(:, 2))
                error(message('vision:pointcloud:invalidROI'));
            end

            roi = double(roi);
            
            % Use bruteforce if there are less than 500 points
            if numel(this.Location)/3 < 500
                if ~this.IsOrganized
                    tf = this.Location(:,1)>=roi(1)&this.Location(:,1)<=roi(4) ...
                        &this.Location(:,2)>=roi(2)&this.Location(:,2)<=roi(5) ...
                        &this.Location(:,3)>=roi(3)&this.Location(:,3)<=roi(6);
                    indices = uint32(find(tf));
                else
                    tf = this.Location(:,:,1)>=roi(1)&this.Location(:,:,1)<=roi(4) ...
                        &this.Location(:,:,2)>=roi(2)&this.Location(:,:,2)<=roi(5) ...
                        &this.Location(:,:,3)>=roi(3)&this.Location(:,:,3)<=roi(6);
                    indices = uint32(find(tf));
                end
            else
                this.buildKdtree();

                indices = this.Kdtree.boxSearch(roi);
            end                        
        end
        
        %==================================================================
        % Obtain a subset of this point cloud object
        %==================================================================
        function ptCloudOut = select(this, varargin)
            %select Select points specified by index.
            %
            %  ptCloudOut = select(ptCloud, indices) returns a pointCloud
            %  object that contains the points selected using linear
            %  indices.
            %
            %  ptCloudOut = select(ptCloud, row, column) returns a pointCloud
            %  object that contains the points selected using row and
            %  column subscripts. This syntax applies only to organized
            %  point cloud (M-by-N-by-3).
            %
            %  Example : Downsample a point cloud with fixed step
            %  -------------------------------------------------
            %   ptCloud = pcread('teapot.ply');
            %
            %   % Downsample a point cloud with fixed step size 
            %   stepSize = 10;
            %   indices = 1:stepSize:ptCloud.Count;
            %
            %   ptCloudOut = select(ptCloud, indices);
            %
            %   pcshow(ptCloudOut);

            narginchk(2, 3);
            
            if nargin == 2
                indices = varargin{1};
                validateattributes(indices, {'numeric'}, ...
                    {'real','nonsparse', 'vector', 'integer'}, mfilename, 'indices');
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
                indices = sub2ind([size(this.Location,1), size(this.Location,2)], row, column);
            end
            
            % Obtain the subset for every property
            [loc, c, nv] = this.subsetImpl(indices);
            angle = this.angle(indices);
            intensity = this.intensity(indices);
            timeStamp = this.timeStamp(indices);
            if ~isempty(this.sensor)
                sensor = this.sensor(indices);
            else
                sensor = uint8.empty;
            end
            
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv, 'intensity', intensity, 'angle', angle, 'timeStamp', timeStamp, 'sensor', sensor);
        end
        
        %==================================================================
        % Obtain a rotated version of this point cloud object
        %==================================================================
        function ptCloudOut = rotate_cloud(this, varargin)
            %select Select points specified by index.
            %
            %  ptCloudOut = select(ptCloud, indices) returns a pointCloud
            %  object that rotates the points of the original cloud

            narginchk(2, 3);
            
            if nargin == 2
                angle_rot = varargin{1};
                validateattributes(angle_rot, {'numeric','double'}, ...
                    {'nonsparse','finite'}, mfilename, 'angle_rot');
            else
                % Subscript syntax is only for organized point cloud
                if ndims(this.Location) ~= 3
                    error(message('vision:pointcloud:organizedPtCloudOnly'));
                end
                    
            end
            
            % Obtain the subset for every property
            loc     = (rotz(angle_rot) * (this.Location'))';
            ang   = this.angle;
            int = this.intensity;
            time = this.timeStamp;
            c = this.Color;
            nv = this.Normal;
            sen = this.sensor;
            
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv, 'intensity', int, 'angle', ang, 'timeStamp', time, 'sensor', sen);
        end
        
        function ptCloudOut = rotate_cloud_x(this, varargin)
            %select Select points specified by index.
            %
            %  ptCloudOut = select(ptCloud, indices) returns a pointCloud
            %  object that rotates the points of the original cloud

            narginchk(2, 3);
            
            if nargin == 2
                angle_rot = varargin{1};
                validateattributes(angle_rot, {'numeric','double'}, ...
                    {'nonsparse','finite'}, mfilename, 'angle_rot');
            else
                % Subscript syntax is only for organized point cloud
                if ndims(this.Location) ~= 3
                    error(message('vision:pointcloud:organizedPtCloudOnly'));
                end
                    
            end
            % Center cloud
            center  =   mean(this.Location);
            loc     =   this.Location - center;
            % Obtain the subset for every property
            loc     =   (rotx(angle_rot) * (loc'))';
            loc     =   loc + center;
            ang   = this.angle;
            int = this.intensity;
            time = this.timeStamp;
            c = this.Color;
            nv = this.Normal;
            sen = this.sensor;
            
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv, 'intensity', int, 'angle', ang, 'timeStamp', time, 'sensor', sen);
        end
        
        function ptCloudOut = translate_cloud(this, varargin)
            %select Select points specified by index.
            %
            %  ptCloudOut = select(ptCloud, indices) returns a pointCloud
            %  object that rotates the points of the original cloud

            narginchk(2, 3);
            translate_point = varargin{1};
%             if nargin == 2
%                 translate_point = varargin{1};
%                 validateattributes(angle_rot, {'numeric','double'}, ...
%                     {'nonsparse','finite'}, mfilename, 'angle_rot');
%             else
%                 % Subscript syntax is only for organized point cloud
%                 if ndims(this.Location) ~= 3
%                     error(message('vision:pointcloud:organizedPtCloudOnly'));
%                 end
%                     
%             end
            
            % Obtain the subset for every property
            loc     =   this.Location - translate_point;
            ang     =   this.angle;
            int     =   this.intensity;
            time    =   this.timeStamp;
            c       =   this.Color;
            nv      =   this.Normal;
            sen     =   this.sensor;
            
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv, 'intensity', int, 'angle', ang, 'timeStamp', time, 'sensor', sen);
        end
        
        function ptCloudOut = align_cloud(this, varargin)


            narginchk(2, 3);
            PCA = varargin{1};
%             if nargin == 2
%                 translate_point = varargin{1};
%                 validateattributes(angle_rot, {'numeric','double'}, ...
%                     {'nonsparse','finite'}, mfilename, 'angle_rot');
%             else
%                 % Subscript syntax is only for organized point cloud
%                 if ndims(this.Location) ~= 3
%                     error(message('vision:pointcloud:organizedPtCloudOnly'));
%                 end
%                     
%             end
            
            % Obtain the subset for every property
            loc     =   this.Location * PCA;
            ang     =   this.angle;
            int     =   this.intensity;
            time    =   this.timeStamp;
            c       =   this.Color;
            nv      =   this.Normal;
            sen     =   this.sensor;
            
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv, 'intensity', int, 'angle', ang, 'timeStamp', time, 'sensor', sen);
        end
        
        function ptCloudOut = clean_cloud(this)

%             narginchk(2, 3);
%             if nargin == 2
%                 translate_point = varargin{1};
%                 validateattributes(angle_rot, {'numeric','double'}, ...
%                     {'nonsparse','finite'}, mfilename, 'angle_rot');
%             else
%                 % Subscript syntax is only for organized point cloud
%                 if ndims(this.Location) ~= 3
%                     error(message('vision:pointcloud:organizedPtCloudOnly'));
%                 end
%                     
%             end
            
            % Obtain the subset for every property
            [~, id, ~]  = unique(this.Location,'rows');
            
            ptCloudOut = select(this, id);
        end
        
        function vis(varargin)
            if (numel(varargin) == 1)
                this = varargin{1};
                pcshow(this.Location);
            elseif (numel(varargin) == 2)
                this = varargin{1};
                colors = varargin{2}; 
                pcshow(this.Location, colors)
            end
        end
        %==================================================================
        % Remove invalid points from this point cloud object
        %==================================================================
        function [ptCloudOut, indices] = removeInvalidPoints(this)
            %removeInvalidPoints Remove invalid points.
            %
            %  [ptCloudOut, indices] = removeInvalidPoints(ptCloud) removes
            %  points whose coordinates contain Inf or NaN. The second
            %  output, indices, is a vector of linear indices indicating
            %  locations of valid points in the point cloud.
            %
            %  Note :
            %  ------
            %  An organized point cloud (M-by-N-by-3) will become
            %  unorganized (X-by-3) after calling this function.
            %               
            %  Example : Remove NaN valued points from a point cloud
            %  ---------------------------------------------------------
            %  ptCloud = pointCloud(nan(100,3))
            %
            %  ptCloud = removeInvalidPoints(ptCloud)

            % Find all valid points
            tf = isfinite(this.Location);
            if ~this.IsOrganized
                indices = (sum(tf, 2) == 3);
            else
                indices = (sum(reshape(tf, [], 3), 2) == 3);
            end

            [loc, c, nv] = this.subsetImpl(indices);
            ptCloudOut = pointCloud_(loc, 'Color', c, 'Normal', nv);
            if nargout > 1
                indices = find(indices);
            end
        end
    end
    
    methods
        %==================================================================
        % Writable Property
        %==================================================================
        function set.Color(this, value) 
            validateattributes(value,{'uint8'}, {'real','nonsparse'});
            if ~isempty(value) && ~isequal(size(value), size(this.Location)) %#ok<MCSUP>
                error(message('vision:pointcloud:unmatchedXYZColor'));
            end
            this.Color = value;
        end
        
        function set.Normal(this, value)
            validateattributes(value,{'single', 'double'}, {'real','nonsparse'});
            if ~isempty(value) && ~isequal(size(value), size(this.Location)) %#ok<MCSUP>
                error(message('vision:pointcloud:unmatchedXYZNormal'));
            end
            if isa(this.Location,'double') %#ok<MCSUP>
                value = double(value);
            else
                value = single(value);
            end
            this.Normal = value;
        end
        
%         function set.Location(this, value)
%             validateattributes(value,{'double'}, {'real','nonsparse'});
%             if ~isempty(value) && ~isequal(size(value), size(this.Location)) 
%                 error(message('vision:pointcloud:unmatchedXYZNormal'));
%             end
%             this.Location = value;
%         end   
        
        %==================================================================
        % Dependent Property
        %==================================================================
        function xlim = get.XLimits(this)
            tf = ~isnan(this.Location);
            if ~this.IsOrganized
                tf = (sum(tf, 2)==3);
                xlim = [min(this.Location(tf, 1)), max(this.Location(tf, 1))];
            else
                tf = (sum(tf, 3)==3);
                X = this.Location(:, :, 1);
                xlim = [min(X(tf)), max(X(tf))];
            end                
        end  
        %==================================================================
        function ylim = get.YLimits(this)
            tf = ~isnan(this.Location);
            if ~this.IsOrganized
                tf = (sum(tf, 2)==3);
                ylim = [min(this.Location(tf, 2)), max(this.Location(tf, 2))];
            else
                tf = (sum(tf, 3)==3);
                Y = this.Location(:, :, 2);
                ylim = [min(Y(tf)), max(Y(tf))];
            end                
        end  
        %==================================================================
        function zlim = get.ZLimits(this)
            tf = ~isnan(this.Location);
            if ~this.IsOrganized
                tf = (sum(tf, 2)==3);
                zlim = [min(this.Location(tf, 3)), max(this.Location(tf, 3))];
            else
                tf = (sum(tf, 3)==3);
                Z = this.Location(:, :, 3);
                zlim = [min(Z(tf)), max(Z(tf))];
            end                
        end
        %==================================================================
        function count = get.Count(this)
            if ~this.IsOrganized
                count = size(this.Location, 1);
            else
                count = size(this.Location, 1)*size(this.Location, 2);
            end                
        end  
    end
    
    methods (Access = public, Hidden)
        %==================================================================
        % helper function to get subset for each property
        %==================================================================
        function [loc, c, nv] = subsetImpl(this, indices)
            if ~isempty(this.Location)
                if ~this.IsOrganized
                    loc = this.Location(indices, :);
                else
                    loc = reshape(this.Location, [], 3);
                    loc = loc(indices, :);
                end
            else
                loc = zeros(0, 3, 'single');
            end
            
            if nargout > 1
                if ~isempty(this.Color)
                    if ~this.IsOrganized
                        c = this.Color(indices, :);
                    else
                        c = reshape(this.Color, [], 3);
                        c = c(indices, :);
                    end
                else
                    c = uint8.empty;
                end
            end
            
            if nargout > 2
                if ~isempty(this.Normal)
                    if ~this.IsOrganized
                        nv = this.Normal(indices, :);
                    else
                        nv = reshape(this.Normal, [], 3);
                        nv = nv(indices, :);
                    end
                else
                     if isa(loc, 'single')
                        nv = single.empty;
                     else
                        nv = double.empty;
                     end
                end
            end
        end
        
        %==================================================================
        % helper function to support multiple queries in KNN search
        % indices, dists: K-by-numQueries
        % valid: 1-by-K
        % Note, the algorithm may return less than K results for each
        % query. Therefore, only 1:valid(n) in n-th column of indices and
        % dists are valid results. Invalid indices are all zeros.
        %==================================================================
        function [indices, dists, valid] = multiQueryKNNSearchImpl(this, points, K)
            % Validate the inputs
            validateattributes(points, {'single', 'double'}, ...
                {'real', 'nonsparse', 'size', [NaN, 3]}, mfilename, 'points');

            if isa(this.Location, 'single')
                points = single(points);
            else
                points = double(points);
            end

            validateattributes(K, {'single', 'double'}, ...
                {'nonsparse', 'scalar', 'positive', 'integer'}, mfilename, 'K');

            K = min(double(K), this.Count);
            
            this.buildKdtree();

            % Use exact search in Kdtree
            searchOpts.checks = 0;
            searchOpts.eps = 0;
            [indices, dists, valid] = this.Kdtree.knnSearch(points, K, searchOpts);
        end
        
        %==================================================================
        % helper function to compute normals
        % normals: the same size of the Location matrix
        %
        % Note, the algorithm uses PCA to fit local planes around a point,
        % and chooses the normal direction (inward/outward) arbitrarily.
        %==================================================================
        function normals = surfaceNormalImpl(this, K)            
            % Reset K if there are not enough points
            K = min(double(K), this.Count);
            
            if this.Count <= 2
                normals = NaN(size(this.Location), 'like', this.Location);
                return;
            end
            
            this.buildKdtree();

            if ~this.IsOrganized
                loc = this.Location;
            else
                loc = reshape(this.Location, [], 3);
            end
            
            % Use exact search in Kdtree
            searchOpts.checks = 0;
            searchOpts.eps = 0;

            % Find K nearest neighbors for each point
            [indices, ~, valid] = this.Kdtree.knnSearch(loc, K, searchOpts);

            % Find normal vectors for each point
            normals = visionPCANormal(loc, indices, valid);
            
            if this.IsOrganized
                normals = reshape(normals, size(this.Location));
            end
        end
    end

    methods (Access = protected)        
        %==================================================================
        % helper function to index data
        %==================================================================
        function buildKdtree(this)
            if isempty(this.Kdtree)                             
                % Build a Kdtree to index the data
                this.Kdtree = vision.internal.Kdtree();                
                createIndex = true;                                        
            elseif this.Kdtree.needsReindex(this.Location)                 
                createIndex = true;                          
            else
                createIndex = false;
            end
            
            if createIndex
                % store rand state prior to indexing
                this.KdtreeIndexState = rng;
                this.Kdtree.index(this.Location);
            end
        end
    end    
    
    methods(Static, Access=private)
        %==================================================================
        % load object 
        %==================================================================
        function this = loadobj(s)
            this = pointCloud_(s.Location,...
                'Color', s.Color,...
                'Normal', s.Normal,...
                'intensity', s.intensity,...
                'timeStamp', s.timeStamp,...
                'angle', s.angle,...
                'sensor', s.sensor);
            
            if isfield(s, {'KdtreeIndexState'}) % added in R2015b
                % set the saved state prior to indexing to recreate the
                % same KD-Tree.
                                
                this.KdtreeIndexState = s.KdtreeIndexState;
                
                % only create index if it was created previously
                if ~isempty(this.KdtreeIndexState)                    
                    sprev = rng(this.KdtreeIndexState);
                    
                    buildKdtree(this);
                    
                    rng(sprev); % restore previous state
                end                         
            end       
        end  
    end
    
    methods(Access=private)
        %==================================================================
        % save object 
        %==================================================================
        function s = saveobj(this)
            % save properties into struct
            s.Location      = this.Location;
            s.Color         = this.Color;
            s.Normal        = this.Normal;
            s.intensity     = this.intensity;
            s.timeStamp     = this.timeStamp;
            s.angle         = this.angle;
            s.sensor        = this.sensor;
            % save the index state to enable recreate duplicate on load.
            s.KdtreeIndexState = this.KdtreeIndexState;
        end
    end
    
    methods(Access=protected)
        %==================================================================
        % copy object 
        %==================================================================
        % Override copyElement method:
        function cpObj = copyElement(obj)
            % Make a copy except the internal Kdtree
            cpObj = pointCloud_(obj.Location, 'Color', obj.Color, 'Normal', obj.Normal,...
                                              'intensity', obj.intensity, 'timeStamp', obj.timeStamp, 'sensor', sen);
        end
    end    
end

%==================================================================
% parameter validation
%==================================================================
function [xyzPoints, C, nv, intensity, timeStamp, angle, sensor] = validateAndParseInputs(varargin)
    % Validate and parse inputs
    narginchk(1, 13);

    parser = inputParser;
    parser.CaseSensitive = false;
    % Parse the arguments according to the format of the first argument

    if ismatrix(varargin{1})
        dims = [NaN, 3];
    else
        dims = [NaN, NaN, NaN, NaN, NaN, 3];
    end

    parser.addRequired('xyzPoints', @(x)validateattributes(x,{'single', 'double'}, {'real','nonsparse','size', dims}));            
    parser.addParameter('Color', uint8([]), @(x)validateattributes(x,{'uint8'}, {'real','nonsparse'}));
    parser.addParameter('Normal', single([]),  @(x)validateattributes(x,{'single', 'double'}, {'real','nonsparse'}));
    parser.addParameter('intensity', double([]), @(x)validateattributes(x,{'single', 'double', 'uint8','uint16'}, {'real','nonsparse'}));
    parser.addParameter('timeStamp', double([]), @(x)validateattributes(x,{'single', 'double'}, {'real','nonsparse'}));
    parser.addParameter('angle', double([]), @(x)validateattributes(x,{'single', 'double', 'int8'}, {'real','nonsparse'}));
    parser.addParameter('sensor', uint8([]), @(x)validateattributes(x,{'uint8'}, {'real','nonsparse'}));
    parser.parse(varargin{:});

    xyzPoints = parser.Results.xyzPoints;
    C = parser.Results.Color;
    nv = parser.Results.Normal; 
    intensity = parser.Results.intensity;
    timeStamp = parser.Results.timeStamp;
    angle = parser.Results.angle;
    sensor = parser.Results.sensor;
    if isa(xyzPoints, 'single')
        nv = single(nv);
    else
        nv = double(nv);
    end
end                                
%==================================================================
% parameter validation for search
%==================================================================
function [doSort, maxLeafChecks] = validateAndParseSearchOption(varargin)
persistent p;
if isempty(p)
    % Validate and parse search options
    p = inputParser;
    p.CaseSensitive = false;
    p.addParameter('Sort', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
    p.addParameter('MaxLeafChecks', inf, @validateMaxLeafChecks);
    parser = p;
else
    parser = p;
end

parser.parse(varargin{:});

doSort = parser.Results.Sort;
maxLeafChecks = parser.Results.MaxLeafChecks;
if isinf(maxLeafChecks)
    % 0 indicates infinite search in internal function
    maxLeafChecks = 0;
end

end
%==================================================================
function validateMaxLeafChecks(value)
    % Validate MaxLeafChecks
    if any(isinf(value))
        validateattributes(value,{'double'}, {'real','nonsparse','scalar','positive'});
    else
        validateattributes(value,{'double'}, {'real','nonsparse','scalar','integer','positive'});
    end
end
