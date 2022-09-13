% Container for cylindrical quantitative structure model block information.
% A subclass of the abstract class QSMB.


% This file is part of QSM-Blocks.
% 
% QSM-Blocks is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QSM-Blocks is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QSM-Blocks.  If not, see <http://www.gnu.org/licenses/>.

classdef QSMBCylindrical < QSMB

    properties (SetAccess=protected)

        % Cylinder-level properties.
        cylinder_start_point = [];
        cylinder_axis        = [];
        cylinder_length      = [];
        cylinder_radius      = [];
        cylinder_parent      = [];
        cylinder_extension   = [];

        cylinder_branch_index    = [];
        cylinder_branch_order    = [];
        cylinder_index_in_branch = [];

        cylinder_end_point = [];
        cylinder_mid_point = [];
        cylinder_is_last   = [];

        cylinder_added = [];

        % Number of branches.
        branch_count = 0;

        % Branch-level properties.
        branch_order  = [];
        branch_parent = [];
        branch_volume = [];
        branch_length = [];
        branch_angle  = [];
        branch_height = [];

        % Triangulated stem.
        has_triangle_stem = false;

        % Vertices of the triangulated stem.
        stem_vertices = [];
        % Faces of the triangulated stem.
        stem_faces = [];

        % Indices of triangulated cylinders.
        triangle_cylinders = [];

    end
    
    properties(Access=private)
        
        % Coordinate system change matrix for each block.
        % Used for optimization in intersection checking.
        cylinder_coordinate_system;
        
    end
    
    properties(Constant)
        
        % Debug flag for easier development.
        debug = false;

        % Default color matrix to use while plotting.
        color_matrix = [ ...
            181,  66,  81; ...
            150, 180,   0; ...
              0, 160, 200; ...
             78,   0, 141; ...
            220, 220,   0; ...
            220,   0, 220; ...
              0, 220, 220; ...
            100, 100,   0; ...
            100,   0, 100; ...
              0, 100, 100 ...
        ]./255;
    end


    methods

        function ob = QSMBCylindrical(varargin)
        % Constructor of class object.
        %
        % QSM = QSMBCylindrical(Sta,Axe,Len,Rad,Par,BI,Add);
        % QSM = QSMBCylindrical(Sta,Axe,Len,Rad,Par,BI,Add,Tri);

            if nargin < 1
                error('Not enough input arguments.');
            end

            if nargin >= 6

                % Start point.
                assert( ...
                    isnumeric(varargin{1}) && ...
                    size(varargin{1},2) == 3, ...
                    'Start point must be a N x 3 matrix.' ...
                );

                ob.block_count = size(varargin{1},1);

                % Axis direction.
                assert( ...
                    isnumeric(varargin{2}) && ...
                    size(varargin{2},1) == ob.block_count && ...
                    size(varargin{2},2) == 3, ...
                    ['Axis direction must be a N x 3 matrix, where N ' ...
                    'matches the row count of the start point ' ...
                    'argument.'] ...
                );

                % Length.
                assert( ...
                    isnumeric(varargin{3}) && ...
                    size(varargin{3},1) == ob.block_count && ...
                    size(varargin{3},2) == 1, ...
                    ['Length must be a N x 1 vector, where N ' ...
                    'matches the row count of the start point ' ...
                    'argument.'] ...
                );

                % Radius.
                assert( ...
                    isnumeric(varargin{4}) && ...
                    size(varargin{4},1) == ob.block_count && ...
                    size(varargin{4},2) == 1, ...
                    ['Radius must be a N x 1 vector, where N ' ...
                    'matches the row count of the start point ' ...
                    'argument.'] ...
                );

                % Cylinder parent.
                assert( ...
                    isnumeric(varargin{5}) && ...
                    size(varargin{5},1) == ob.block_count && ...
                    size(varargin{5},2) == 1, ...
                    ['Parent must be a N x 1 vector, where N ' ...
                    'matches the row count of the start point ' ...
                    'argument.'] ...
                );

                % Branch index.
                assert( ...
                    isnumeric(varargin{6}) && ...
                    size(varargin{6},1) == ob.block_count && ...
                    size(varargin{6},2) == 1, ...
                    ['Branch index must be a N x 1 vector, where N ' ...
                    'matches the row count of the start point ' ...
                    'argument.'] ...
                );

                ob.cylinder_start_point = varargin{1};
                ob.cylinder_axis = varargin{2};
                ob.cylinder_length = varargin{3};
                ob.cylinder_radius = varargin{4};
                ob.cylinder_parent = varargin{5};
                ob.cylinder_branch_index = varargin{6};

                if nargin > 6
                    % Added.
                    Added = varargin{7};

                    if isempty(Added)
                        Added = false(ob.block_count,1);
                    end

                    assert( ...
                        islogical(Added) && ...
                        size(Added,1) == ob.block_count && ...
                        size(Added,2) == 1, ...
                        ['Added must be a N x 1 logical vector, where ' ...
                        'N matches the row count of the start point ' ...
                        'argument.'] ...
                    );

                    ob.cylinder_added = Added;
                else
                    ob.cylinder_added = false(ob.block_count,1);
                end

                if nargin > 7
                    % Triangle stem.
                    assert( ...
                        isstruct(varargin{8}) && ...
                        isfield(varargin{8},'vertices') && ...
                        isfield(varargin{8},'faces') && ...
                        isfield(varargin{8},'cylinders'), ...
                        [ ...
                            'Triangles must be a struct with the ' ...
                            'following fields: ''vertices'', ' ...
                            '''faces'', ''cylinders''.' ...
                        ] ...
                    );

                    ob.stem_vertices = varargin{8}.vertices;
                    ob.stem_faces = varargin{8}.faces;
                    ob.triangle_cylinders = varargin{8}.cylinders;

                    ob.has_triangle_stem = true;
                else
                    ob.has_triangle_stem = false;
                end

            % Old cell-based model format.
            elseif nargin == 1 && iscell(varargin{1})

                ModelData = varargin{1};

                assert( ...
                    isnumeric(ModelData{1}) && ...
                    size(ModelData{1},2) >= 11, ...
                    ['First element of input cell-array should have ' ...
                    'at least 11 columns.'] ...
                );


                % Set number of blocks.
                ob.block_count = size(ModelData{1},1);

                % Set cylinder-level properties.
                ob.cylinder_start_point = ModelData{1}(:,3:5);
                ob.cylinder_axis        = ModelData{1}(:,6:8);
                ob.cylinder_length      = ModelData{1}(:,2);
                ob.cylinder_radius      = ModelData{1}(:,1);
                ob.cylinder_parent      = ModelData{1}(:,9);

                ob.cylinder_branch_index    = ModelData{1}(:,11);

                if size(ModelData{1},2) > 13
                    ob.cylinder_added = ModelData{1}(:,14);
                else
                    ob.cylinder_added = false(ob.block_count,1);
                end

                % Add custom tree properties.
                if numel(ModelData) > 2
                    TreeData = ModelData{3};

                    if length(TreeData) > 10
                        ob.set_property('DBH CYL', TreeData(11));
                    end
                    if length(TreeData) > 11
                        ob.set_property('DBH TRI', TreeData(12));
                    end
                end

            % New struct-based model format.
            elseif nargin == 1 && isstruct(varargin{1})

                ModelData = varargin{1};

                assert( ...
                    isfield(ModelData,'cylinder') && ...
                    isfield(ModelData.cylinder,'start') && ...
                    isfield(ModelData.cylinder,'axis') && ...
                    isfield(ModelData.cylinder,'length') && ...
                    isfield(ModelData.cylinder,'radius') && ...
                    isfield(ModelData.cylinder,'parent') && ...
                    isfield(ModelData.cylinder,'branch'), ...
                    ['Struct type input argument should have struct ' ...
                    'field ''cylinder'', which in turn sould have ' ...
                    'fields: ''start'', ''axis'', ''length'', ' ...
                    '''radius'', ''parent'' and ''branch''.'] ...
                );

                % Set number of blocks.
                ob.block_count = length(ModelData.cylinder.radius);

                % Set cylinder-level properties.
                ob.cylinder_start_point = ModelData.cylinder.start;
                ob.cylinder_axis        = ModelData.cylinder.axis;
                ob.cylinder_length      = ModelData.cylinder.length;
                ob.cylinder_radius      = ModelData.cylinder.radius;
                ob.cylinder_parent      = ModelData.cylinder.parent;

                ob.cylinder_branch_index = ModelData.cylinder.branch;

                if isfield(ModelData.cylinder,'added')
                    ob.cylinder_added = ModelData.cylinder.added;
                else
                    ob.cylinder_added = false(ob.block_count,1);
                end

                % Add custom tree properties.
                if isfield(ModelData,'treedata')

                    % Tree-level data.
                    TreeData = ModelData.treedata;

                    % Store different DBH estimates as custom properties.
                    if isfield(TreeData, 'DBHcyl')
                        ob.set_property('DBH CYL', TreeData.DBHcyl);
                    end
                    if isfield(TreeData, 'DBHtri')
                        ob.set_property('DBH TRI', TreeData.DBHtri);
                    end
                end

                % Check if trianglulated stem data is present.
                if isfield(ModelData,'triangulation') && ...
                    isfield(ModelData.triangulation,'vert') && ...
                    not(isempty(ModelData.triangulation.vert))

                    % Triangle vertices.
                    ob.stem_vertices = ModelData.triangulation.vert;
                    % Triangle faces.
                    ob.stem_faces = ModelData.triangulation.facet;

                    % Indices of cylinders that are covered by the 
                    % triangulated stem.
                    ob.triangle_cylinders = ...
                        1:ModelData.triangulation.cylid;
                    %-

                    % Check triangulation flag to true.
                    ob.has_triangle_stem = true;

                end

            elseif nargin == 1 && ischar(varargin{1}) && ...
                strcmp(varargin{1}, 'example')

                ob.cylinder_start_point = [ ...
                     0.000  0.000  0.000; ...
                     0.000  0.000  1.000; ...
                     1.000  0.000  2.000; ...
                     0.000  0.000  1.000; ...
                     0.000  1.000  2.000; ...
                     0.000  0.000  1.000; ...
                     0.000 -1.000  2.000; ...
                     0.000  0.000  1.000; ...
                    -1.000  0.000  2.000; ...
                ];

                ob.cylinder_axis = [ ...
                     0.000  0.000  1.000; ...
                     0.707  0.000  0.707; ...
                     1.000  0.000  0.000; ...
                     0.000  0.707  0.707; ...
                     0.000  1.000  0.000; ...
                     0.000 -0.707  0.707; ...
                     0.000 -1.000  0.000; ...
                    -0.707  0.000  0.707; ...
                    -1.000  0.000  0.000; ...
                ];

                ob.cylinder_length = [ ...
                    1.000; ...
                    1.414; ...
                    1.000; ...
                    1.414; ...
                    1.000; ...
                    1.414; ...
                    1.000; ...
                    1.414; ...
                    1.000; ...
                ];

                ob.cylinder_radius = [ ...
                    0.300; ...
                    0.200; ...
                    0.100; ...
                    0.200; ...
                    0.100; ...
                    0.200; ...
                    0.100; ...
                    0.200; ...
                    0.100; ...
                ];

                ob.cylinder_parent = [ ...
                    0; ...
                    1; ...
                    2; ...
                    1; ...
                    4; ...
                    1; ...
                    6; ...
                    1; ...
                    8; ...
                ];

                ob.cylinder_branch_index = [ ...
                    1; ...
                    1; ...
                    1; ...
                    2; ...
                    2; ...
                    3; ...
                    3; ...
                    4; ...
                    4; ...
                ];

                ob.block_count = length(ob.cylinder_parent);

                ob.cylinder_added = false(ob.block_count,1);

            else
                error('Unknown input argument type.');
                

            end
            % End of various input types.
            
            % Set default values for twig distribution.
            ob.fun_twig_distribution = @default_twig_param_dist;
            ob.twig_length_limits = [0.02, 0.05];
            
            % Initialize cylinder coordinate system variable. Used during
            % triangle intersection detection. If matrix contains NaNs, the
            % respective matrix is computed during execution.
            ob.cylinder_coordinate_system = nan(3,3,ob.block_count);

            % Compute derived values and set datatypes.
            ob = ob.uncompress();

        end

        function Origin = origin(ob)
        % Return the origin point of the tree model. If the models has 
        % blocks (cylinders), then the origin is the start point of the 
        % first cylinder. Otherwise, the Cartesian origin (0,0,0) is
        % returned.

            if ob.block_count > 0
                % Origin of the tree is the start point of the first
                % cylinder.
                Origin = ob.cylinder_start_point(1,:);
            else
                Origin = zeros(1,3);
            end

        end

        function ob = translate(ob, Delta)
        % Move cylinder start points by given delta-vector.
        
            assert( ...
                (size(Delta,1) == 1 && size(Delta,2) == 3) || ...
                (size(Delta,1) == 3 && size(Delta,2) == 1), ...
                'Translate vector must be a three-element vector.' ...
            );

            ob.cylinder_start_point(1:ob.block_count,:) = bsxfun( ...
                @plus, ...
                ob.cylinder_start_point(1:ob.block_count,:), ...
                Delta ...
            );

        end

        function ob = rotate(ob, ang, ax, mode)
        % Rotate model around given axis by given angle. 
        % 
        % Inputs:
        %
        % ang       Angle of rotation in degrees.
        %
        % ax        Axis vector around which to rotate. z-axis by default.
        %           When empty, default is used.
        %
        % mode      Rotation coordinate system: 
        %
        %           'local' [default]:  Rotate around base of the tree.
        %           'global':           Rotate around world origin.
        %

            % Set default axis as z-axis.
            if nargin < 3 || isempty(ax)
                ax = [0, 0, 1];
            end

            if nargin < 4
                mode = 'local';
            end

            % Convert degrees to radians.
            ang = deg2rad(ang);

            % Get rotation matrix. Angle in radians.
            R = rotation_matrix(ax, ang);

            switch mode
                case 'local'

                    % Temporary copy of start points.
                    SP = ob.cylinder_start_point(1:ob.block_count,:);

                    % Origin of tree, i.e., first start point.
                    sp1 = SP(1,:);

                    % Ensure points start from origin.
                    SP = bsxfun(@minus, SP, sp1);

                    % Scale start points by factor.
                    SP = SP * R;

                    % Move back to previous origin.
                    SP = bsxfun(@plus, SP, sp1);

                    % Replace data in object.
                    ob.cylinder_start_point(1:ob.block_count,:) = SP;
                case 'global'
                    % Transform start points.
                    ob.cylinder_start_point(1:ob.block_count,:) = ...
                        ob.cylinder_start_point(1:ob.block_count,:) * R;
                    %-
                otherwise
                    error(['Unknown rotate mode: "' mode '"']);
            end

            % Transform axes.
            ob.cylinder_axis(1:ob.block_count,:) = ...
                ob.cylinder_axis(1:ob.block_count,:) * R;
            %-

        end

        function ob = scale(ob, s)
        % Scale model by given factor.
        %
        % Inputs:
        %
        % s         Scaling factor.

            % Transform start points.
            ob.cylinder_radius(1:ob.block_count,:) = ...
                ob.cylinder_radius(1:ob.block_count,:) * s;
            %-

            % Transform start points.
            ob.cylinder_length(1:ob.block_count,:) = ...
                ob.cylinder_length(1:ob.block_count,:) * s;
            %-

            % Temporary copy of start points.
            SP = ob.cylinder_start_point(1:ob.block_count,:);

            % Origin of tree, i.e., first start point.
            sp1 = SP(1,:);

            % Ensure points start from origin.
            SP = bsxfun(@minus, SP, sp1);

            % Scale start points by factor.
            SP = SP.*s;

            % Move back to previous origin.
            SP = bsxfun(@plus, SP, sp1);

            % Replace data in object.
            ob.cylinder_start_point(1:ob.block_count,:) = SP;

        end

        function Val = volume(ob, varargin)
        % Compute total or sub-volumes of the tree model. Without optional 
        % arguments the method returns total tree volume. The optional
        % argument can either be a string ('total', 'stem', 'branch', etc.)
        % or a custom logical vector with the same length as the QSM block
        % count. 
        %
        % Examples:
        % 
        % % Compute total volume.
        % Volume = QSM.volume();
        % % OR
        % Volume = QSM.volume('total');
        %
        % % Compute stem volume.
        % Volume = QSM.volume('stem');
        % % OR
        % Volume = QSM.volume('trunk');
        %
        % % Compute branch volume.
        % Volume = QSM.volume('branch');
        % % OR
        % Volume = QSM.volume('branches');
        %
        % % Compute crown volume.
        % Volume = QSM.volume('crown');
        %
        % % Custom volume.
        % IStemBranch = QSM.cylinder_branch_order == 1;
        % Volume = QSM.volume(IStemBranch);

            Val = ob.compute_property('volume', varargin{:});

        end

        function Val = area(ob, varargin)
        % Compute total or sub-areas of the tree model. Without optional 
        % arguments the method returns total tree area. The optional
        % argument can either be a string ('total', 'stem', 'branch', etc.)
        % or a custom logical vector with the same length as the QSM block
        % count. 
        %
        % Examples:
        % 
        % % Compute total area.
        % Area = QSM.area();
        % % OR
        % Area = QSM.area('total');
        %
        % % Compute stem area.
        % Area = QSM.area('stem');
        % % OR
        % Area = QSM.area('trunk');
        %
        % % Compute branch area.
        % Area = QSM.area('branch');
        % % OR
        % Area = QSM.area('branches');
        %
        % % Compute crown area.
        % Area = QSM.area('crown');
        %
        % % Custom area.
        % IStemBranch = QSM.cylinder_branch_order == 1;
        % Area = QSM.area(IStemBranch);

            Val = ob.compute_property('area', varargin{:});

        end

        function Val = length(ob, varargin)
        % Compute total or sub-lengths of the tree model. Without optional 
        % arguments the method returns total tree length. The optional
        % argument can either be a string ('total', 'stem', 'branch', etc.)
        % or a custom logical vector with the same length as the QSM block
        % count. 
        %
        % Examples:
        % 
        % % Compute total length.
        % Length = QSM.length();
        % % OR
        % Length = QSM.length('total');
        %
        % % Compute stem length.
        % Length = QSM.length('stem');
        % % OR
        % Length = QSM.length('trunk');
        %
        % % Compute branch length.
        % Length = QSM.length('branch');
        % % OR
        % Length = QSM.length('branches');
        %
        % % Compute crown length.
        % Length = QSM.length('crown');
        %
        % % Custom length.
        % IStemBranch = QSM.cylinder_branch_order == 1;
        % Length = QSM.length(IStemBranch);

            Val = ob.compute_property('length', varargin{:});

        end

        function Val = height(ob, varargin)
        % Compute total or sub-heights of the tree model. Without optional 
        % arguments the method returns total tree height. The optional
        % argument can either be a string ('total', 'stem', 'branch', etc.)
        % or a custom logical vector with the same height as the QSM block
        % count. 
        %
        % Examples:
        % 
        % % Compute total height.
        % Height = QSM.height();
        % % OR
        % Height = QSM.height('total');
        %
        % % Compute stem height.
        % Height = QSM.height('stem');
        % % OR
        % Height = QSM.height('trunk');
        %
        % % Compute branch height.
        % Height = QSM.height('branch');
        % % OR
        % Height = QSM.height('branches');
        %
        % % Compute crown height.
        % Height = QSM.height('crown');
        %
        % % Custom height.
        % IStemBranch = QSM.cylinder_branch_order == 1;
        % Height = QSM.height(IStemBranch);

            Val = ob.compute_property('height', varargin{:});

        end

        function Val = count(ob, varargin)
        % Compute total or sub-counts of the cylinders in the tree model. 
        % Without optional arguments the method returns total cylinder 
        % count. The optional argument can either be a string 
        % ('total', 'stem', 'branch', etc.) or a custom logical vector with 
        % the same count as the QSM block count. 
        %
        % Examples:
        % 
        % % Compute total count.
        % Count = QSM.count();
        % % OR
        % Count = QSM.count('total');
        %
        % % Compute stem count.
        % Count = QSM.count('stem');
        % % OR
        % Count = QSM.count('trunk');
        %
        % % Compute branch count.
        % Count = QSM.count('branch');
        % % OR
        % Count = QSM.count('branches');
        %
        % % Compute crown count.
        % Count = QSM.count('crown');
        %
        % % Custom count.
        % IStemBranch = QSM.cylinder_branch_order == 1;
        % Count = QSM.count(IStemBranch);

            Val = ob.compute_property('count', varargin{:});

        end

        function ISub = get_cylinder_set(ob, str)
        % Get a logical vector indicating which of the models cylinders
        % are part of a set, defined by a case-insensitive string
        % descriptor. Available sets are:
        %
        % 'Total'                   Select all cylinders.
        %
        % 'Trunk' OR 'Stem'         Selected as branch order == 0.
        %
        % 'Branch' OR 'Branches'    Selected as branch order > 0.
        %
        % 'Crown'                   Select branches that are part of the
        %                           tree crown but try to avoid including
        %                           short branch stubs lower in the stem.

            switch lower(str)

                % All cylinders.
                case 'total'
                    ISub = true(ob.block_count,1);

                % Select stem by branch order.
                case {'trunk', 'stem'}
                    ISub = ob.cylinder_branch_order == 0;
                    
                % Select branches by branch order.
                case {'branch', 'branches'}
                    ISub = ob.cylinder_branch_order > 0;

                % Get crown as a logical vector.
                case 'crown'
                    ISub = ob.get_crown();

                % Other string identifiers are unknown.
                otherwise
                    error(['Unknown cylinder set: ''' str '''']);
            end

        end

        function Val = compute_property(ob, Property, varargin)
        % Compute summative tree properties either using all cylinders 
        % or a subset of cylinders. 
        % 
        % Inputs:
        %
        % Property       Name of the property to compute. Possible values:
        %                'volume', 'area', 'length', 'height', 'count'.
        %                The argument is case-insensitive.
        %
        % [Optional]     String or logical vector to select subset of 
        %                cylinders. Possible string values:
        %                'total', 'stem', 'trunk', 'branch', 'branches'
        %                'crown'.
        %                Logical vector length should match block count.

            assert(ischar(Property), 'Property should be a string.');

            % If no optional argument given, return total property value.
            if nargin == 2
                varargin = {'total'};
            end

            % Get first (and only) argument.
            Arg = varargin{1};

            % If argument is string.
            if ischar(Arg)

                % Get subset of cylinders from string descriptor.
                ISub = ob.get_cylinder_set(Arg);

            % If argument is a logical (vector).
            elseif islogical(Arg)

                % Check that vector length matches cylinder/block count.
                assert(...
                    length(Arg(:)) == ob.block_count, ...
                    'Logical vector length should match block count.' ...
                );

                % Use extra argument as subset selector.
                ISub = Arg(:);

            else
                error([...
                    'Unknown optional argument. ' ...
                    'Argument should be a string or a logical vector.'...
                ]);
            end

            % Cylinders are selected. By default all cylinders are
            % selected.
            if any(ISub)

                switch lower(Property)
                    case 'volume'
                        % Cylinder volume. [Pi*R^2*L]
                        Val = pi*sum(ob.cylinder_radius(ISub).^2 ...
                                    .*ob.cylinder_length(ISub));
                        %-

                    case 'area'
                        % Cylinder area. [2*Pi*R*L]
                        Val = 2*pi*sum(ob.cylinder_radius(ISub) ...
                                        .*ob.cylinder_length(ISub));
                        %-

                    case 'length'
                        % Cylinder length.
                        Val = sum(ob.cylinder_length(ISub));

                    case 'count'
                        % Number of selected cylinders.
                        Val = nnz(ISub);

                    case 'height'
                        % Compute extreme heights from start and end
                        % points.
                        Heights = [...
                            ob.cylinder_start_point(ISub,3); ...
                            ob.cylinder_end_point(ISub,3) ...
                        ];

                        Val = max(Heights) - min(Heights);

                    otherwise
                        error([...
                            'Unknown property ''' ...
                            Property ...
                            ''' to compute.' ...
                        ]);
                end

            else
                % If no cylinder are selected, set value to zero.
                Val = 0;
            end

        end

        function DBH = dbh(ob)
        % Compute tree DBH from cylinder model. Here DBH is defined as
        % the diameter of the last cylinder below 1.3 meters. Length 
        % is computed along the tree stem.
        %
        % Examples:
        %
        % DBH = QSM.dbh();

            % Get stem cylinders.
            IStem = ob.cylinder_branch_order == 0;

            % Check that at least one stem cylinder exists.
            assert(any(IStem), ...
                'Model has no stem cylinders. Unable to compute DBH.');
            %-

            % Compute cumulative cylinder length along the stem.
            LenSum = cumsum(ob.cylinder_length(IStem));

            % "Shift" cumulative sum one index towards the top.
            % This achieves the cumulative length at the cylinder
            % bottom point.
            LenSum = [0; LenSum(1:end-1)];

            % Find last cylinder below 1.3 meters. As at least one
            % stem cylinder exists, FIND should always return an
            % index.
            JHeight = find(LenSum < 1.3, 1, 'last');

            % Diameter is twice the radius of the selected cylinder.
            DBH = ob.cylinder_radius(JHeight)*2;

        end

        function ICrown = get_crown(ob)
        % Compute the indices of QSM cylinders that are part
        % of the tree crown.
        %
        % Outputs:
        %
        % ICrown        Logical vector of length 
        %               [ob.block_count x 1] where TRUE means
        %               cylinder is part of the tree crown.

            % Logical index for stem cylinders.
            IStem = ob.cylinder_branch_order == 0;

            % Form initial crown set.
            for iOrder = 3:-1:0
                ICrown = ob.cylinder_branch_order > iOrder;

                if any(ICrown)
                    break;
                end
            end

            % Indices of crown extensions.
            JCrownExt = ICrown;

            % Search new cylinders in the parent direction (backwards).
            while true

                JCrownExt = unique(ob.cylinder_parent(JCrownExt));

                % Check if root cylinder is included (no parent).
                if JCrownExt(1) == 0

                    % If new cylinders other than the root was found
                    % remove root and continue.
                    if length(JCrownExt) > 1
                        JCrownExt = JCrownExt(2:end);
                    % Otherwise stop this phase.
                    else
                        break;
                    end
                end

                % New cylinders are not part of the crown and not 
                % part of the stem.
                INew = not(ICrown(JCrownExt)) & not(IStem(JCrownExt));

                % If no new cylinders found, stop adding.
                if not(any(INew))
                    break;
                end
                
                % Set new cylinders as part of crown.
                ICrown(JCrownExt(INew)) = true;

            end

            % Flag to see if new cylinders were still added.
            FNewFound = true;

            % Indices of new cylinders.
            JNew = find(ICrown);

            % Search new cylinders in the child direction (forward).
            while FNewFound
                FNewFound = false;

                % New cylinders have parents in the current crown set,
                % and are not part of the crown yeat.
                INew = ismember(ob.cylinder_parent, JNew) & not(ICrown);

                % New cylinders found.
                if any(INew)

                    % Update current crown set.
                    ICrown(INew) = true;
                    % Get indices of new cylinders.
                    JNew = find(INew);
                    % Flag that new cylinders found, thus continue.
                    FNewFound = true;
                end
            end

            % Find minimum crown attachment to stem.
            [CrownStartHeight,~] = min(...
                ob.cylinder_start_point(ICrown & ...
                    ob.cylinder_branch_order == 1,3)...
            );

            % Add all branch cylinders above lowest attachment point to
            % crown.
            ICrown(not(IStem) & ...
                ob.cylinder_start_point(:,3) >= CrownStartHeight) = true;

        end

        function Props = get_block_properties(ob,PropNames,vargin)
        % Function to read (or compute, in some cases) properties of
        % cylindrical blocks. Property names are given in a cell array of
        % strings. It is also possible to restrict the property
        % computations to certain cylinders with the variable <indices>.
        % The property values are returned as fields of a struct.
        %
        % Examples:
        %
        % props = ob.get_block_properties(properties)
        % props = ob.get_block_properties(properties, indices)

            % Check if filtering by cylinder ID is required.
            if nargin > 2
                JCyl = vargin{1};
            else
                JCyl = 1:ob.block_count;
            end

            % Number of properties to compute.
            NProp = numel(PropNames);

            % Properties to return.
            Props = [];

            for iProp = 1:NProp

                % Current property label.
                Label = lower(PropNames{iProp});

                switch Label

                    % Basic properties that are already stored in the
                    % object. Return appropriate matrices.
                    case {'start_point', 'axis', 'length',...
                          'radius', 'parent', 'extension',...
                          'branch_index','branch_order',...
                          'index_in_branch','is_last'}
                        %-

                        Props.(Label) = ob.(['cylinder_' Label])(JCyl,:);

                    % Relative height of cylinder, normalized by tree
                    % extreme points.
                    case 'relative_height'
                        
                        CylMeanHeight = ...
                                  mean([ob.cylinder_start_point(JCyl,3),...
                                        ob.cylinder_end_point(JCyl,3)],2);
                        %-
            
                        Props.(Label) = ...
                               (CylMeanHeight - ob.tree_limits(1,3)) / ...
                               (ob.tree_limits(2,3) - ob.tree_limits(1,3));
                        %-

                    % Relative position of cylinder along the branch.
                    % Computed from the end point of the cylinder. Thus,
                    % the relative position of the last cylinder in a
                    % branch is always one.
                    case 'relative_position'

                        BranchIndex   = ob.cylinder_branch_index(JCyl);
                        IndexInBranch = ob.cylinder_index_in_branch(JCyl);
                        Length        = ob.cylinder_length(JCyl);


                        RelativePosition = zeros(length(JCyl),1);

                        % Compute cylinder position on branch.
                        % Iterate over all branch numbers.
                        for i = 1:max(BranchIndex)
                            
                            % Indices of cylinder in given branch number.
                            IBranch = (BranchIndex == i);

                            if not(any(IBranch))
                                continue;
                            end
                            
                            % Sort from start to finnish.
                            [~,JCylinder] = sort(IndexInBranch(IBranch));
                            
                            JBranch = find(IBranch);
                            
                            % Cumulative position of end points.
                            CumulativeLength = ...
                                        cumsum(Length(JBranch(JCylinder)));

                            % Total length.
                            TotalBranchLength = CumulativeLength(end);

                            % Relative position of cylinders.
                            RelativePosition(JBranch) = CumulativeLength...
                                                      / TotalBranchLength;
                        end

                        Props.(Label) = RelativePosition;

                end

            end

        end
        
        
        function [TwigStart, TwigEnd] = generate_twigs(ob, LeafParent)
        % Generate twigs connecting leaves to the cylinders. The input
        % <LeafParent> is a sorted vector of cylinder indices of parents 
        % of the leaves. A single cylinder can occur multiple times. The 
        % function returns the start and end points of the twigs.
            
            % Number of twigs to genereate.
            NTwig = length(LeafParent);
            
            % Parameters of cylinders with twigs. Same parent can be
            % repeated multiple times.
            sp      = ob.cylinder_start_point(LeafParent,:);
            ax      = ob.cylinder_axis(LeafParent,:);
            h       = ob.cylinder_length(LeafParent);
            r       = ob.cylinder_radius(LeafParent);
            is_last = ob.cylinder_is_last(LeafParent);
            
            % Parameters of the twigs. Row = single twig.
            % Columns:
            %   1: position on axis on the inverval [0,1].
            %   2: position on radial axis on the inverval [0,1]. 
            %      Equals 1 if on side, less than 1 when on head.
            %   3: rotation around axis. 0 upwards, Pi downwards.
            %   4: twig elevation. -pi backwards, Pi towards tip.
            %   5: twig azimuth. -Pi/2 left, Pi/2 right.
            %   6: twig length.
            TwigParam = ob.fun_twig_distribution(sp, ax, h, r, is_last,...
                                                 ob.twig_length_limits);
            %-
            
            % Index of the last parent cylinder, used to check if cylinder
            % coordinate system needs to be updated.
            LastParent = [];
            
            % Initialize return values.
            TwigStart = zeros(NTwig,3);
            TwigEnd   = zeros(NTwig,3);
            
            % Iterate over twigs.
            for iTwig = 1:NTwig
                
                % Always compute on first run, and when parent index
                % changes.
                if iTwig == 1 || LeafParent(iTwig) ~= LastParent
                    
                    % Update parent index.
                    LastParent = LeafParent(iTwig);
                
                    % Axis pointing to the side, defining radial direction.
                    if all(ax(iTwig,:) == [0 0 1])
                        TwigSideBase = [0 1 0];
                        TwigUpBase = [1 0 0];
                    else
                        TwigSideBase = cross([0 0 1],ax(iTwig,:));
                        TwigSideBase = TwigSideBase/norm(TwigSideBase);

                        TwigUpBase = cross(ax(iTwig,:),TwigSideBase);
                    end
                    
                end
                
                % Debug: plot cylinder coordinate system axis vectors.
                if ob.debug
                    figure(2);
                    conf = [1 4];
                    TwigFront = ax(iTwig,:);
                    TwigUp = TwigUpBase;
                    TwigSide = TwigSideBase;
                    plot_axes(conf,1,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Rotation matrix to rotate around cylinder axis.
                Rax = rotation_matrix(ax(iTwig,:),TwigParam(iTwig,3));

                % Rotate required coordinate axes around cylinder axis.
                TwigUp   = (Rax*TwigUpBase')';
                TwigSide = (Rax*TwigSideBase')';

                % Compute twig start point.
                TwigStart(iTwig,:) = sp(iTwig,:) ...
                                   + ax(iTwig,:)*h(iTwig)...
                                                *TwigParam(iTwig,1) ...
                                   + TwigUp     *r(iTwig)...
                                                *TwigParam(iTwig,2);
                %-
                
                if ob.debug
                    plot_axes(conf,2,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Elevation parameter has to be inverted if the twig is
                % connected to the tip of the cylinder.
                if TwigParam(iTwig,2) == 1
                    ElParam = -TwigParam(iTwig,4);
                else
                    ElParam = TwigParam(iTwig,4);
                end
                
                % Rotation matrix to rotate twig either towards
                % cylinder axis.
                Rel = rotation_matrix(TwigSide,ElParam);

                % Rotate required coordinate axes around side axis.
                TwigUp    = (Rel*TwigUp')';
                TwigFront = (Rel*ax(iTwig,:)')';
                
                % Debug: plot updated coordinate system.
                if ob.debug
                    plot_axes(conf,3,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
                % Debug: plot twig direction unit vector.
                if ob.debug
                    
                    colors = eye(3);
                    E = [TwigUp; TwigSide; TwigFront];
                    
                    figure(2);
                    for i = 1:3
                        plott([TwigStart(iTwig,:); ...
                               TwigStart(iTwig,:)+E(i,:)],...
                              '-','Color',colors(i,:));
                        %-
                        hold on;
                    end
                    hold off;
                    axis equal;
                    
                end

                % Twig is attached to the envelope.
                if TwigParam(iTwig,2) == 1
                                        
                    % Rotation matrix to rotate to wanted azimuth, i.e.,
                    % turn the twig axis off cylinder center.
                    Raz = rotation_matrix(TwigFront,TwigParam(iTwig,5));
                    
                    % Rotate required coordinate axes around 
                    TwigUp = (Raz*TwigUp')';
                    
                    % Debug: update for later plot.
                    if ob.debug
                        TwigSide = (Raz*TwigSide')';
                    end
                    
                    % Compute twig end point.
                    TwigEnd(iTwig,:) = TwigStart(iTwig,:) ...
                                     + TwigUp*TwigParam(iTwig,6);
                    %-
                    
                    
                % Twig is attached to the tip of the cylinder.
                else
                    
                    % Rotation matrix to rotate to wanted azimuth, i.e.,
                    % turn the twig axis off cylinder center.
                    Raz = rotation_matrix(TwigUp,-TwigParam(iTwig,5));
                    
                    % Rotate required coordinate axes around.
                    TwigFront = (Raz*TwigFront')';
                    
                    % Debug: update for later plot.
                    if ob.debug
                        TwigSide = (Raz*TwigSide')';
                    end
                    
                    % Compute twig end point.
                    TwigEnd(iTwig,:) = TwigStart(iTwig,:) ...
                                     + TwigFront*TwigParam(iTwig,6);
                    %-
                    
                end
                
                % Plot updated coordinate system.
                if ob.debug
                    plot_axes(conf,4,TwigStart(iTwig,:),...
                              TwigUp,TwigSide,TwigFront);
                    %-
                end
                
            end
            
        end


        function fIntersect = block_triangle_intersection(ob, JBlock, Tris)
        % Check if triangles intersect with any of the blocks with given
        % indices. Returns True if at least one of the given triangles
        % intersects any of the given blocks.

            % Number of cylinders.
            NCyl = length(JBlock);
            
            % Number of triangles.
            NTri = size(Tris,1);

            % Flag for intersection.
            fIntersect = false;

            % Iterate over each cylinder-triangle pair.
            for iCyl = 1:NCyl

                for iTri = 1:NTri

                    % Check intersection of each pair.
                    fIntersect = ob.cylinder_triangle_intersect(...
                                               JBlock(iCyl),...
                                               Tris(iTri,:),...
                                               false);
                    %-
                    
                    % If an intersection occurs, skip the rest and return.
                    if fIntersect
                        return;
                    end
                end


            end

        end
        
        % Function to check the intersection of a single cylinder and as
        % single triangle.
        intersect = cylinder_triangle_intersect(ob, jBlock, Tri,...
                                                edgehit, debug)
        %-


        % toVoxels: Return center coordinates of voxels that contain parts
        % of the blocks.
        function BlockVoxelization = toVoxels(ob, edge, minp, maxp)

            % Initialize voxelized space object.
            BlockVoxelization = CubeVoxelization(edge, minp, maxp);
            
            % Iterate over blocks in the model to populate voxelization.
            for iCyl = 1:ob.block_count

                % Get cube coordinates of voxels occupied by cylinder.
                cc = ob.CylinderToVoxels(iCyl, edge, minp);
                
                cc = min(cc,BlockVoxelization.size);
                cc = max(cc,[1 1 1]);

                % Add cylinder index to occupied voxels.
                BlockVoxelization.add_object_by_cc(cc,iCyl);

            end

        end

        % External implementation.
        export(ob, Format, File, varargin)


        function export_blender(ob, File, Prec, Origin, varargin)
        % Print cylinder model parameters to file, e.g., for exporting 
        % to Blender, using the Blender QSM import addon.

            % Flag if extra parameters are given to print to the file.
            ExtraData = [];

            if nargin > 4
                for iArgin = 1:numel(varargin)
                    ExtraData = cat(2, ExtraData, varargin{iArgin});
                end
            end

            % Set origin override as zeros, if not given.
            if nargin < 4 || isempty(Origin)
                Origin = zeros(1,3);
            end

            ob.export(...
                'blender', ...
                File, ...
                'Precision', Prec, ...
                'Origin', Origin, ...
                'ExtraData', ExtraData ...
            );

        end

        function h = plot_model(ob,varargin)
        % Plot the cylinder model as a patch object on the current axis.
        % Cylinder geometry is converted into vertices and faces. Apart 
        % from the QSMBCylindrical object the method does not have other
        % required input parameters, but many optional ones do exist.
        % The optional arguments are given in name-value pairs unless
        % stated otherwise below. Returns the PATCH object.
        %
        %
        % Optional inputs:
        %
        % 'FaceCount'           Set the integer count or faces to compute 
        %                       on each cylinder. Can be either a single
        %                       integer or a two-element vector of 
        %                       integers, in which case the face count is 
        %                       interpolated linearly based on the relative
        %                       radii, i.e., smaller cylinders have less
        %                       faces. Default value is 5.
        %
        % 'Quads'               Boolean. If TRUE use quads instead of
        %                       triangles on cylinder envelope faces.
        %                       By default TRUE. Can be defined without
        %                       value pair, then value is interpreted as 
        %                       TRUE.
        %
        % 'Closed'              Boolean. If TRUE a triangle fan of faces is 
        %                       generated to "close" the cylinder top and
        %                       bottom. By default FALSE.  Can be defined 
        %                       without value pair, then value is 
        %                       interpreted as TRUE.
        %
        % 'TriangleStem'        Boolean. If TRUE triangulated stem data is
        %                       used rather than cylinder geometry when
        %                       plotting the model. By default FALSE.
        %                       Can be defined without value pair, then
        %                       value is interpreted as TRUE. Throws a
        %                       warning if model has no triangulated stem.
        %
        % 'Origin'              Three element vector that is used to move
        %                       the base of the tree model. Empty by
        %                       default. Use QSM.cylinder_start_point(1,:)
        %                       to move base to [0,0,0].
        %
        % 'Filter'              Boolean vertor with one element per
        %                       cylinder or a string descriptor passed to 
        %                       the QSMBCylindrical.get_cylinder_set() 
        %                       method. Only plot cylinders that are set
        %                       to TRUE.
        %
        % 'PlotOptions'         Cell array of parameters to be passed to
        %                       the PATCH command. By default 'EdgeColor'
        %                       is set to 'none' and 'FaceColor' is set to
        %                       'flat'.
        %
        % 'ColorMatrix'         [N x 3] matrix of true colors to be used to
        %                       color cylinders according to discrete
        %                       properties such as branch order or branch
        %                       index. The constant property color_matrix
        %                       contains the default color matrix. A color
        %                       matrix can have any number of rows as it is
        %                       indexed with looping.
        %
        % 'Color'               Either a [N x 1] vector or a [N x 3] matrix
        %                       of color values, with N matching the
        %                       cylinder count. Describing the color values
        %                       for each cylinder. Value can also be a
        %                       single three-element color vector to be
        %                       used for all cylinders. Takes precedence
        %                       over 'ColorSource' when both are present.
        %
        % 'ColorSource'         String descriptor of property to use for
        %                       coloring vertices, faces or cylinders. Use
        %                       as an alternative for the 'Color' 
        %                       attribute. 
        %                       QSMBCylindrical.get_color_distribution()
        %                       method is used to compute the color 
        %                       distribution. Possible values are:
        %
        %       'None'                      No color based on source.
        %
        %       'CylinderLength'            Length of a cylinder.
        %      
        %       'CylinderIsLast'            Highlight last cylinders in 
        %                                   each branch.
        %      
        %       'CylinderAdded'             Highlight cylinders that were
        %                                   added after reconstruction.
        %      
        %       'CylinderRadius'            Cylinder radius.
        %      
        %       'CylinderOrder',
        %       'CylinderBranchOrder',
        %       'BranchOrder',
        %       'Order'                     Branch order of cylinder.
        %      
        %       'CylinderBranchIndex',
        %       'CylinderIndex',
        %       'BranchIndex'               Separate branches with color.
        %      
        %       'BranchLength'              Longer branches have higher 
        %                                   index.
        %      
        %       'BranchAngle'               Branching angle at branch base.
        %      
        %       'BranchHeight'              Mean height of branch.
        %
        %       'VertexHeight'              Height of vertices. This is the
        %                                   default value.
        %
        %       'VertexRadius'              Horizontal distance of vertex
        %                                   from tree base.
        %
        %       'FaceHeight'                Mean height of face vertices.
        %
        %
        % Examples:
        %
        % % Default plotting.
        % QSM.plot_model();
        %
        % % Change face count on single cylinders.
        % QSM.plot_model('FaceCount',[5 10]);
        %
        % % Color by branch order.
        % QSM.plot_model('ColorSource','Order');
        %
        % % Only plot stem.
        % QSM.plot_model('Filter','Stem');
        %
        % % Return PATCH object and manipulate.
        % hQSM = QSM.plot_model();
        % set(hQSM, 'EdgeColor',[1 0 0]);

            % Skip method when empty.
            if ob.block_count == 0
                warning('No cylinders in the model. Exiting.');
                return;
            end

            % Compressed models can not be operated on.
            assert(not(ob.compressed), ...
                ['Model is compressed. ' ...
                'Use the UNCOMPRESS method to uncompress.'] ...
            );

            % Struct to hold and modify method parameters.
            Param = struct();

            % Flags

            % Use quads instead of triangles when true.
            Param.FQuads = true;
            % Close tops and bottoms with triangle fans when true.
            Param.FClosed = false;

            % Number of side faces computed for a single cylinder.
            Param.FaceCount = 5;

            % Logical vector for selecting a subset of cylinders.
            Param.ICyl = [];

            % Color data for cylinders.
            Param.CylColor = [];

            % Vector for moving the tree model start point.
            % By default no tranlation.
            Param.Origin = [0, 0, 0];

            % Default plot options that are passed to the draw commands.
            DefaultPlotOptions = { ...
                'EdgeColor','none', ...
                'FaceColor','flat' ...
            };

            % Plot extra options.
            Param.PlotOptions = {};

            % Default color source. Color by vertex height.
            Param.ColorSource = 'vertexheight';

            % Use default color matrix of object.
            Param.ColorMatrix = ob.color_matrix;

            % Flag: use trianglulated stem vertices and faces
            % and exclude respective cylinders.
            Param.FTriStem = false;

            % List of acceptable arguments for current function.
            Accepted = { ...
                'facecount','closed','quads', 'plotoptions', ...
                'filter','color','origin', 'colorsource', ...
                'colormatrix', 'trianglestem' ...
            };

            % Parse optional input arguments.
            Param = ob.parse_arguments(Param, varargin, Accepted);

            % Convert cylinder parameters to vertices and faces.
            [Vert, Faces, JVertCyl] = ob.compute_geometry( ...
                Param.FaceCount, ...
                Param.ICyl, ...
                Param.FQuads, ...
                Param.FClosed ...
            );

            % Color defined by 'Color' attribute.
            if not(isempty(Param.CylColor))

                % Use color values defined by 'Color' attribute.
                ColorDist = Param.CylColor(JVertCyl,:);

            % If color source is set by 'ColorSource' attribute.
            elseif not(strcmp(Param.ColorSource,'none'))

                % Try to compute color distribution. May throw error if
                % color source string is unknown.
                ColorDist = ob.get_color_distribution(...
                    Vert, ...
                    Faces, ...
                    JVertCyl, ...
                    Param.ColorSource, ...
                    Param.ColorMatrix ...
                );

            % Otherwise, no color values.
            else
                ColorDist = [];
            end

            % If usage of triangulated stem is requested, convert
            % cylinder geometry to triangulated geometry and update
            % possible color data.
            if Param.FTriStem

                % Vertex-cylinder mapping is excluded as output because
                % not used after this point.
                [Vert, Faces, ~, ColorDist] = ob.convert_stem_vert( ...
                    Vert, ...
                    Faces, ...
                    JVertCyl, ...
                    ColorDist ...
                );
            end

            % If color prerent.
            if not(isempty(ColorDist))

                % If color distribution size matches vertex count, 
                % use interpolated vertex-based coloring.
                if size(ColorDist,1) == size(Vert,1)
                    PlotParam = {'FaceColor','interp'};

                % Otherwise color distribution size should match face
                % count and default flat coloring can be used.
                else
                    PlotParam = {};
                end

                % Add color property to patch extra options.
                Param.PlotOptions = cat( ...
                    1, ...
                    PlotParam(:), ...
                    {'FaceVertexCData'; ColorDist}, ...
                    Param.PlotOptions(:) ...
                );
            end

            % If extra plotting options were given, add them after the 
            % default options.
            if not(isempty(Param.PlotOptions))
                Param.PlotOptions = cat( ...
                    1, ...
                    DefaultPlotOptions(:), ...
                    Param.PlotOptions(:) ...
                );
            % Otherwise only use the default options.
            else
                Param.PlotOptions = DefaultPlotOptions;
            end

            % If location is offset, move vertices.
            if any(Param.Origin)
                Vert = bsxfun(@minus, Vert, Param.Origin);
            end

            % Store figure hold status.
            FResetHold = ~ishold;

            % Plot faces with plotting options.
            h = patch( ...
                'Vertices', Vert, ...
                'Faces', Faces, ...
                Param.PlotOptions{:} ...
            );

            % Set axis scaling equal.
            axis equal;
            view(3);

            % Restore the hold status if necessary.
            if FResetHold
                hold off;
            end

        end

        function [Vert, Faces, JVertCyl, JFaceCyl] = compute_geometry( ...
            ob, ...
            NFace, ...
            ICyl, ...
            FQuads, ...
            FClosed ...
        )
        % Compute vertex and face based geometry for selected cylinders
        % of the cylinder model. 
        %
        % Inputs:
        %
        % ob            The QSMBCylindrical object.
        %
        % NFace         Number of side faces on a single cylinder. Can be
        %               either a scalar integer or a pair of integers. In
        %               the latter case face count is scaled between the 
        %               given lower and upper limits based with the linear
        %               scaling parameter interpolated between the minimum
        %               and maximum cylinder radii values of the included
        %               cylinders.
        %
        % ICyl          Logical vector selecting included cylinders. Can be
        %               empty.
        %
        % FQuads        Logical flag for determining face shape. When TRUE
        %               Side faces are quads and triangles otherwise.
        %
        % FClosed       Logical flag for determining whether to close the
        %               cylinder hull with top and bottom faces. Additional
        %               faces are generated when TRUE.
        %
        % Outputs:
        %
        % Vert          Matrix of cylinder vertices. Each row contains the
        %               Three coordinates of a single vertex.
        %
        % Faces         Matrix of face indices. Each row contains the 
        %               indices of vertices forming the face. If quads are
        %               used the matrix has four columns, and if the 
        %               cylinder ends are closed the respective rows 
        %               contain NaNs as the last elements. If quads are not
        %               used, the matrix will have three columns and no
        %               NaNs.
        %
        % JVertCyl      Cylinder index of each vertex. The index values
        %               are in the full set of cylinders, i.e., the closed
        %               interval [0, ob.block_count], regardless of 
        %               optional filtering.
        %
        % JFaceCyl      Cylinder index of each face. The index values
        %               are in the full set of cylinders, i.e., the closed
        %               interval [0, ob.block_count], regardless of 
        %               optional filtering.
        %
        % Examples:
        %
        % % Geometry of all cylinders with five faces each.
        % [Vert, Faces]  = QSM.compute_geometry(5);
        % 
        % % Geometry of all cylinders with between five and ten faces each.
        % [Vert, Faces]  = QSM.compute_geometry([5 10]);
        % 
        % % Geometry of the stem.
        % IStem = QSM.cylinder_branch_order == 0;
        % [Vert, Faces]  = QSM.compute_geometry([5 10], IStem);
        %
        % % Only use triangles.
        % [Vert, Faces]  = QSM.compute_geometry([5 10], [], false);
        %
        % % Close cylinders with triangle fans.
        % [Vert, Faces]  = QSM.compute_geometry([5 10], [], false, true);
        %
        % % Also return vertex cylinder indices.
        % [Vert, Faces, JVertCyl]  = QSM.compute_geometry([5 10]);
        %

            % No closed cylinder ends by default.
            if nargin < 5
                FClosed = false;
            end

            % Use quads by default.
            if nargin < 4
                FQuads = true;
            end

            % By default all cylinders are converted to vertices and faces.
            if nargin < 3 || isempty(ICyl)
                ICyl = true(ob.block_count,1);
            end

            % Check that face count is an integer or a pair of integers.
            assert( ...
                length(NFace(:)) <= 2, ...
                [ ...
                    'Face count must be an integer or a ' ...
                    'two-element vector of integers.' ...
                ] ...
            );

            % Radius extrema.
            rmin = min(ob.cylinder_radius(ICyl));
            rmax = max(ob.cylinder_radius(ICyl));

            % Face count or maximum face count if multiple elements.
            NMax = max(NFace);

            % Number of selected cylinders.
            NCyl = nnz(ICyl);

            % Face count scaling based on individual cylinder radii.
            % Activated if two inequal face counts given and radius
            % values vary.
            if length(NFace) > 1 && NFace(1) ~= NFace(2) && rmin ~= rmax
                
                % Linear scaling parameter for each included cylinder.
                p = (ob.cylinder_radius(ICyl)-rmin)./(rmax-rmin);

                % Linear interpolation of face count values using the
                % scaling parameter and upper and lower face counts.
                FaceCounts = round(p*(NFace(2) - NFace(1)) + NFace(1));

            % Otherwise use the same face count for each cylinder.
            else
                FaceCounts = repmat(NFace(1),ob.block_count,1);
            end

            % Width of vertex matrix.
            NProp = 3;

            % Higher width for quads.
            if FQuads
                NProp = 4;
            end

            % Number of added vertices.
            NVertex = 0;
            % Number of added faces.
            NFace = 0;

            % All vertices. Initialize with the theoretical maximum size
            % with all cylinders having max face count and closed tops
            % and bottoms.
            Vert = zeros((double(FClosed)+NMax)*2*NCyl,3);

            % All faces. Same initialization with the maximum theoretical
            % face count.
            Faces = zeros((1+double(FClosed)*NMax*NCyl),NProp);

            if nargout > 2
                % Cylinder index for each vertex. Same initialization as
                % vertex matrix.
                JVertCyl = zeros((double(FClosed)+NMax)*2*NCyl,1);
            end

            if nargout > 3
                % Cylinder index for each face. Same initialization as
                % face matrix.
                JFaceCyl = zeros((1+double(FClosed)*NMax*NCyl),1);
            end

            % Go through all cylinders.
            for iCyl = 1:NCyl

                % Skip un-selected cylinders.
                if not(ICyl(iCyl))
                    continue;
                end

                % Face count for current cylinder.
                np = FaceCounts(iCyl);

                % Cylinder parameters. Start point, axis direction,
                % length and radius.
                sp = ob.cylinder_start_point(iCyl,:);
                ax = ob.cylinder_axis(iCyl,:);
                h  = ob.cylinder_length(iCyl);
                r  = ob.cylinder_radius(iCyl);
                
                % Generate vertices for a unit cylinder with given face
                % count.
                [x, y, z] = cylinder(1,np);
                
                % Skip replicated point.
                x = x(:,1:end-1);
                y = y(:,1:end-1);
                z = z(:,1:end-1);

                % Normalize axis and convert to row vector.
                ax = ax(:)'/norm(ax);

                % Convert to row vector.
                sp = sp(:)';

                % Adjust radius and height with cylinder parameters.
                x = r*x;
                y = r*y;
                z = h*z;
                
                % Unit vector for adjusting cylinder rotation.
                u = [0 0 1];

                % Find rotation axis.
                raxis = cross(u,ax);
                raxis = raxis/norm(raxis);

                % Rotation angle.
                angle = acos(ax(3)/norm(ax));

                % Size of vertex matrices. Used when reshaping later.
                s = size(x);

                % If rotation required.
                if any(raxis)

                    % Rotate vertex points.
                    X = [x(:) y(:) z(:)]*rotation_matrix(raxis,angle)';

                    % Reshape vectors back to matrices.
                    x = reshape(X(:,1),s);
                    y = reshape(X(:,2),s);
                    z = reshape(X(:,3),s);
                end

                % Translate to cylinder origin.
                x = x' + sp(1);
                y = y' + sp(2);
                z = z' + sp(3);
                
                % Form cylinder vertex matrix.
                CylVert = [x(:) y(:) z(:)];

                % Define faces.

                % Quad faces.
                if FQuads
                    
                    % Init face matrix.
                    CylFaces = zeros(np,NProp);
                   
                    % Go through faces.                   
                    for i = 1:np
                        
                        % Last face.
                        if i == np
                            CylFaces(i,:) = [i 1 np+1 np+i];
                        % Otherwise.
                        else
                            CylFaces(i,:) = [i i+1 np+i+1 np+i];
                        end
                        
                    end
                
                % Triangle faces.
                else
                    % Init face matrix.
                    CylFaces = zeros(2*np,NProp);
                    
                    % Go through first triangles of faces.
                    for i = 1:np
                       
                        % Last face (loop).
                        if i == np
                            CylFaces(i,:) = [i 1 np+i];
                        % Otherwise.
                        else
                            CylFaces(i,:) = [i i+1 np+i];
                        end
                        
                    end
                    
                    % Go through second triangles of faces.
                    for i = 1:np
                       
                        % Last face (loop).
                        if i == np
                            CylFaces(np+i,:) = [np+i 1 np+1];
                        % Otherwise.
                        else
                            CylFaces(np+i,:) = [np+i i+1 np+i+1];
                        end
                        
                    end
                    
                end
                
                % If the cylinder is closed, add top and bottom 
                % triangle fans.
                if FClosed

                    % Add vertex at the start point and at the end point.
                    CylVert = cat(1, CylVert, sp, sp+ax*h);
                    
                    % Top faces.
                    TopFace = zeros(np,3);
                    % Bottom faces.
                    BotFace = zeros(np,3);
                    
                    % One triangle on either end for each side face.
                    for i = 1:np 
                        
                        % Last face (loop).
                        if i == np
                            BotFace(i,:) = [i 2*np+1 1];
                            TopFace(i,:) = [np+i np+1 2*np+2];

                        % Otherwise.
                        else
                            BotFace(i,:) = [i 2*np+1 i+1];
                            TopFace(i,:) = [np+i np+i+1 2*np+2];
                        end
                        
                    end

                    % Combine for possible padding.
                    TBFaces = [BotFace; TopFace];

                    if FQuads
                        % Pad with NaNs if quads are used.
                        TBFaces = cat(2, TBFaces, nan(size(TBFaces,1),1));
                    end

                    % Append top and bottom faces to all cylinder faces.
                    CylFaces = cat(1, CylFaces, TBFaces);
                    
                end

                % Offset face indices with existing vertex count.
                CylFaces = CylFaces + NVertex;

                % Number of new vertices.
                NNewFace = size(CylFaces,1);

                % Number of new vertices.
                NNewVertex = size(CylVert,1);

                % Store new vertices.
                Vert(NVertex+1:NVertex+NNewVertex,:) = CylVert;

                % Store new faces.
                Faces(NFace+1:NFace+NNewFace,:) = CylFaces;

                if nargout > 2
                    % Assign cylinder index for each vertex.
                    JVertCyl(NVertex+1:NVertex+NNewVertex) = iCyl;
                end

                if nargout > 3
                    % Assign cylinder index for each face.
                    JFaceCyl(NFace+1:NFace+NNewFace) = iCyl;
                end


                % Update total counts.
                NVertex = NVertex + NNewVertex;
                NFace = NFace + NNewFace;

            end

            % Trim empty lines from return values.
            Vert = Vert(1:NVertex,:);
            Faces = Faces(1:NFace,:);

            if nargout > 2
                JVertCyl = JVertCyl(1:NVertex);
            end

            if nargout > 3
                JFaceCyl = JFaceCyl(1:NFace);
            end


        end
        
        function ob = compress(ob)
        % Clear out all data from object that can be derived from other
        % variables. Allows the the object size to be compressed without
        % loss of data.

            % Call parent method first.
            ob = compress@QSMB(ob);

            % Cylinder end point.
            ob.cylinder_end_point = zeros(0,3);
            % Cylinder middle point.
            ob.cylinder_mid_point = zeros(0,3);
            % Flag: cylinder is last in branch.
            ob.cylinder_is_last = [];
            % Index of extension.
            ob.cylinder_extension = [];
            % Index inside branch.
            ob.cylinder_index_in_branch = [];
            % Branch order.
            ob.cylinder_branch_order = [];

            % Branch-level data.
            ob.branch_order  = [];
            ob.branch_parent = [];
            ob.branch_volume = [];
            ob.branch_length = [];
            ob.branch_angle  = [];
            ob.branch_height = [];

            % Bounding box.
            ob.tree_limits = [];

        end

        function ob = uncompress(ob)
        % Recomputute data that was cleared with the .COMPRESS() method.

            % Call parent method first.
            ob = uncompress@QSMB(ob);

            % Number of cylinders.
            NBlock = ob.block_count;

            % Initialize properties with correct size.
            ob.cylinder_extension = zeros(NBlock,1);
            ob.cylinder_index_in_branch = zeros(NBlock,1);
            ob.cylinder_branch_order = nan(NBlock,1);

            % Set branch order of cylinders without parents as zeros.
            ob.cylinder_branch_order(ob.cylinder_parent == 0) = 0;

            % Index inside brach for cylinders without parents is one
            % as they are the starting elements.
            ob.cylinder_index_in_branch(ob.cylinder_parent == 0) = 1;

            % Go through cylinders.
            for iBlock = 1:NBlock

                % Branch index of current cylinder.
                jBranch = ob.cylinder_branch_index(iBlock);

                % Find extension, i.e., cylinder that is part of the same
                % branch and has this cylinder as parent. Should be only
                % one.
                jExt = find(...
                    ob.cylinder_parent == iBlock & ...
                    ob.cylinder_branch_index == jBranch, ...
                    1,'first' ...
                );

                % If not empty, extension was found and should be recorded.
                if not(isempty(jExt))
                    ob.cylinder_extension(iBlock) = jExt;
                end

                % Parent index of current cylinder.
                jPar = ob.cylinder_parent(iBlock);

                % If cylinder has parent, compute values. Otherwise data
                % should already be present.
                if jPar > 0

                    % If cylinder parent has a set branch order.
                    if not(isnan(ob.cylinder_branch_order(jPar)))

                        % If cylinder is an extension, i.e.,
                        % part of the same branch.
                        if ob.cylinder_branch_index(iBlock) == ...
                            ob.cylinder_branch_index(jPar)

                            % Set branch order as parent order.
                            ob.cylinder_branch_order(iBlock) = ...
                                ob.cylinder_branch_order(jPar);
                            %-

                        % Cylinder is a branch of the parent.
                        else

                            % Set branch order as parent order + 1.
                            ob.cylinder_branch_order(iBlock) = ...
                                ob.cylinder_branch_order(jPar) + 1;
                            %-

                        end

                    end

                    % If the index in branch is set for the parent.
                    if ob.cylinder_index_in_branch(jPar) > 0

                        % If cylinder is an extension, i.e.,
                        % part of the same branch.
                        if ob.cylinder_branch_index(iBlock) == ...
                            ob.cylinder_branch_index(jPar)

                            % Incease parent index by one.
                            ob.cylinder_index_in_branch(iBlock) = ...
                                ob.cylinder_index_in_branch(jPar) + 1;
                            %-

                        % Cylinder is a branch of the parent.
                        else

                            % Set as first cylinder in branch.
                            ob.cylinder_index_in_branch(iBlock) = 1;

                        end

                    end
                end

            end

            % Compute end point.
            ob.cylinder_end_point = ob.cylinder_start_point ...
                                  + bsxfun(@times,ob.cylinder_axis,...
                                                  ob.cylinder_length);
            %-
            
            % Compute middle point.
            ob.cylinder_mid_point = ob.cylinder_start_point ...
                                  + bsxfun(@times,ob.cylinder_axis,...
                                                  ob.cylinder_length/2);
            %-

            % Initialize last cylinder flag as false.
            ob.cylinder_is_last = false(ob.block_count,1);

            % All branch indices.
            BranchIndex = unique(ob.cylinder_branch_index);

            % Number of branch indices.
            NBranch = length(BranchIndex);

            % Initialize branch-level data.

            ob.branch_count = NBranch;
            ob.branch_order  = zeros(NBranch,1);
            ob.branch_parent = zeros(NBranch,1);
            ob.branch_volume = zeros(NBranch,1);
            ob.branch_length = zeros(NBranch,1);
            ob.branch_angle  = zeros(NBranch,1);
            ob.branch_height = zeros(NBranch,1);

            % Go through branch indices.
            for iBranch = 1:NBranch

                % Logical vector selecting cylinder in current branch.
                IBranch = ob.cylinder_branch_index == BranchIndex(iBranch);

                % Find first two cylinders.
                JFirst = find(IBranch,2,'first');

                % Branch order from first cylinder in branch.
                ob.branch_order(iBranch) = ...
                    ob.cylinder_branch_order(JFirst(1));
                %-

                % Use method for computing branch volume. Multiply by
                % one-thousand to get liters.
                ob.branch_volume(iBranch) = ...
                    1000*ob.compute_property('volume', IBranch);
                %-

                % Use method to compute total length of branch cylinders.
                ob.branch_length(iBranch) = ...
                    ob.compute_property('length', IBranch);
                %-

                % Select child axis for angle computation. If cylinder
                % has at least two cylinders and the first cylinder was 
                % added during reconstruction, use the second cylinder
                % axis.
                if length(JFirst) > 1 && ob.cylinder_added(JFirst(1))
                    ChildAxis = ob.cylinder_axis(JFirst(2),:);

                % Otherwise use the axis of the first cylinder.
                else
                    ChildAxis = ob.cylinder_axis(JFirst(1),:);
                end

                % If the first cylinder has a parent assigned.
                if ob.cylinder_parent(JFirst(1)) > 0

                    % The branch parent is the branch index of the 
                    % branch of the parent of the first cylinder.
                    ob.branch_parent(iBranch) = ...
                        ob.cylinder_branch_index(...
                            ob.cylinder_parent(JFirst(1)...
                        )...
                    );

                    % Branch angle based on the child and parent axis
                    % directions.
                    ob.branch_angle(iBranch) = 180 / pi * acos(...
                        ChildAxis*ob.cylinder_axis(...
                            ob.cylinder_parent(JFirst(1)),:)' ...
                    );
                end

                % Branch height is computed as the mean of the heights
                % of the branch cylinders.
                ob.branch_height(iBranch) = ...
                    mean(ob.cylinder_mid_point(IBranch,3));
                %-

                % Find maximum index inside branch for current branch.
                MaxIndex = max(ob.cylinder_index_in_branch(IBranch));

                % Find last cylinder in branch.
                jLast = IBranch & ob.cylinder_index_in_branch == MaxIndex;

                % Set flag for last cylinder.
                ob.cylinder_is_last(jLast) = true;

            end

            % Find extreme points, add and substract radius to each
            % cylinder start and end point to be sure.
            Ep_pr = bsxfun( ...
                @plus, ...
                ob.cylinder_end_point, ...
                ob.cylinder_radius ...
            );
            
            Ep_mr = bsxfun( ...
                @plus, ...
                ob.cylinder_end_point, ...
                ob.cylinder_radius ...
            );
            
            Sp_pr = bsxfun( ...
                @plus, ...
                ob.cylinder_start_point, ...
                ob.cylinder_radius ...
                );
            
            Sp_mr = bsxfun( ...
                @plus, ...
                ob.cylinder_start_point, ...
                ob.cylinder_radius ...
            );
            
            
            % Upper and lower limits for the tree in the z-axis.
            ob.tree_limits = [ ...
                min([Ep_pr; Ep_mr; Sp_pr; Sp_mr],[],1); ...
                max([Ep_pr; Ep_mr; Sp_pr; Sp_mr],[],1) ...
            ];

            % Ensure correct data types.
            ob = ob.set_datatypes();            

        end

    end

    methods(Access=protected)

        function cc = CylinderToVoxels(ob, iCyl, edge, minp)
        % Find voxel indices of all voxels a single cylinder occupies. The
        % cylinder is defined by an index, and the voxelization with an
        % edge length and a minimum point.

            % Get cylinder properties to create point samples.
            sp = ob.cylinder_start_point(iCyl,:);
            ax = ob.cylinder_axis(iCyl,:);
            h  = ob.cylinder_length(iCyl,:);
            r  = ob.cylinder_radius(iCyl,:);

            % Generate test points in and on cylinder.
            P = QSMBCylindrical.VoxelSamples(edge, sp, ax, h, r);

            % Compute cube coordinates of each sample point.
            cc = floor(bsxfun(@minus,P,minp)./edge) + 1;

            % Only include each voxel once.
            cc = unique(cc,'rows');

        end

        function ob = set_datatypes(ob)
        % Convert the fields to correct data types. 

            ob.cylinder_start_point = single(ob.cylinder_start_point);
            ob.cylinder_axis = single(ob.cylinder_axis);
            ob.cylinder_length = single(ob.cylinder_length);
            ob.cylinder_radius = single(ob.cylinder_radius);
            ob.cylinder_end_point = single(ob.cylinder_end_point);
            ob.cylinder_mid_point = single(ob.cylinder_mid_point);

            % True / false.
            ob.cylinder_added = logical(ob.cylinder_added);
            ob.cylinder_is_last = logical(ob.cylinder_is_last);

            % Maximum 255.
            ob.cylinder_branch_order = uint8(ob.cylinder_branch_order);

            % Maximum 65535.
            ob.cylinder_parent = uint16(ob.cylinder_parent);
            ob.cylinder_extension = uint16(ob.cylinder_extension);
            ob.cylinder_branch_index = uint16(ob.cylinder_branch_index);
            ob.cylinder_index_in_branch = ...
                uint16(ob.cylinder_index_in_branch);
            %-

            % Maximum 255.
            ob.branch_order  = uint8(ob.branch_order);

            % Maximum 65535.
            ob.branch_parent = uint16(ob.branch_parent);
            ob.branch_count = uint16(ob.branch_count);
            ob.block_count = uint16(ob.block_count);

            ob.branch_volume = single(ob.branch_volume);
            ob.branch_length = single(ob.branch_length);
            ob.branch_angle  = single(ob.branch_angle );
            ob.branch_height = single(ob.branch_height);

        end

        function [Vert, Faces, JVertCyl, CDist] = convert_stem_vert( ...
            ob, ...
            Vert, ...
            Faces, ...
            JVertCyl, ...
            CDist ...
        )
        % Replace cylinder representation with the triangulated stem
        % geometry in vertex-face format. The method also updates the 
        % optional color distribtution by assigning the color value of
        % the closest replaced face/vertex onto the replacing geometry.
        % Vertex cylinder index is updated on the precess as well.

            % If no triangle present, skip.
            if not(ob.has_triangle_stem)
                return;
            end

            % No color data is required by default.
            if nargin < 5
                CDist = [];
            end

            % Check whether vertex or face-level color data is present.
            % Assign a flag for later.
            if size(CDist,1) == size(Vert,1)
                FVertLevel = true;
            else
                FVertLevel = false;
            end

            % Indices of cylinders that have been triangulated.
            JTriCyl = ob.triangle_cylinders;

            % Logical indices of vertices of triangulated cylinders.
            ITriVert = ismember(JVertCyl,JTriCyl);

            % Logical indices of faces of triangulated cylinders.
            ITriFace = any(ismember(Faces,find(ITriVert)),2);

            % Vertices or cylinders to be replaced.
            CylVert = Vert(ITriVert,:);

            % Vertices of replacing geometry.
            TriVert = ob.stem_vertices;

            % Faces of replacing geometry.
            StemFaces = ob.stem_faces;

            % Number of new vertices.
            NTriVert = size(TriVert,1);

            % Index of mapping each new vertex onto a cylinder vertex.
            JVertMap = zeros(NTriVert,1);

            % Go through all new vertices.
            for iTriVert = 1:NTriVert

                % Compute distance to all replaced cylinder vertices.
                PDist = bsxfun(@minus,CylVert,TriVert(iTriVert,:));

                % Find minimum distances.
                [~, jMap] = min(sum(PDist.^2,2));

                % Assign mapping.
                JVertMap(iTriVert) = jMap;

            end

            % If color data is present, update.
            if not(isempty(CDist))

                % If vertex-level color data.
                if FVertLevel

                    % Color values of replaced cylinder vertices.
                    TriCDist = CDist(ITriVert,:);
                    % Color values of new vertices, received with the 
                    % vertex mapping vector.
                    TriCDist = TriCDist(JVertMap,:);

                    % Remove replaced cylinder vertex color data.
                    CDist = CDist(not(ITriVert),:);

                % Face-level color data.
                else

                    % Number of new faces.
                    NTriFace = size(StemFaces,1);

                    % Index of mapping each new face onto a cylinder face.
                    JFaceMap = zeros(NTriFace,1);

                    % Go through all the new faces.
                    for iTriFace = 1:NTriFace

                        % "Closest" face is the one that has the most
                        % closest vertices.
                        JFaceMap(iTriFace) = mode( ...
                            JVertMap(StemFaces(iTriFace,:)) ...
                        );

                    end

                    % Colro values of replaced cylinder faces.
                    TriCDist = CDist(ITriFace,:);
                    % Color values of new faces, received with the face
                    % mapping vector.
                    TriCDist = TriCDist(JFaceMap,:);

                    % Remove replaced cylinder face color data.
                    CDist = CDist(not(ITriFace),:);
                end

                % Append new color data to the end of the old color data.
                CDist = cat(1, CDist, TriCDist);
                
            end

            % Offset vector, number of deleted vertices before every 
            % vertex.
            Offset = cumsum(ITriVert);
            
            % If NaNs present, indexing does not work.
            if any(isnan(Faces(:)))
                
                % Add zero offset to end.
                Offset = cat(1,Offset,0);
                
                % Create temp copy.
                TempFaces = Faces;
                
                % Set NaNs to have offset index pointing to zero.
                TempFaces(isnan(Faces)) = length(Offset);
                
                % Offset face vertex indices with number of removals before
                % each vertex.
                Faces = Faces - Offset(TempFaces);
                
            else
                % Offset face vertex indices with number of removals before
                % each vertex.
                Faces = Faces - Offset(Faces);
            end         

            % Remove replaced vertices.
            Vert = Vert(not(ITriVert),:);

            % Remove replaced faces.
            Faces = Faces(not(ITriFace),:);

            % Remove replaced vertex-cylinder mappings.
            JVertCyl = JVertCyl(not(ITriVert));

            % Number of vertices before adding new, used for face
            % index offset.
            NVert = size(Vert,1);

            % Add new vertices.
            Vert = cat(1, Vert, TriVert);
            
            % Check if new face matrix needs to be padded with NaNs.
            if size(Faces,2) > 3
                StemFaces = cat( ...
                    2, ...
                    StemFaces, ...
                    nan(size(StemFaces,1),1) ...
                );
            end
            
            % Offset vertex indices by existing vertex count.
            Faces = cat(1, Faces, StemFaces + NVert);

        end

        % External implementations.
        Param = parse_arguments(ob, Param, Args, Accepted)
        Dist = get_color_distribution(ob, ...
            Vertices, ...
            Faces, ...
            JVertCyl, ...
            ColorSource, ...
            ColorMatrix ...
        )
        
    end
    
    methods(Static)

        function print_number_line(fid, Val, Format, Lead, Trail, Sep)
        % Print a single line of data into a file stream. 
        %
        % Inputs:
        %
        % fid       File stream idenfifier.
        % Val       Vector of numeric values.
        % Format    Print format string.
        % Lead      Leading string to print at line beginning. Empty
        %           by default.
        % Trail     Trailing string to print at the end of the line. Empty
        %           by default.
        % Sep       Separator to print between values. By default a space.
            
            % Default paremeters.
            if nargin < 6
                Sep = ' ';
            end

            if nargin < 5
                Trail = '';
            end

            if nargin < 4
                Lead = '';
            end

            % Ensure vector.
            Val = Val(:);

            % Number of values to print.
            NVal = length(Val);

            % Leading string.
            fprintf(fid, '%s', Lead);

            % Go through values.
            for iVal = 1:NVal

                % Do not print NaNs. Important for printing face indices
                % when both quads and triangles are present.
                if isnan(Val(iVal))
                    break;
                end

                % If multiple formats, use format with respective index.
                if iscell(Format)
                    fprintf(fid, Format{iVal}, Val(iVal));

                % Otherwise use single format for all values.
                else
                    fprintf(fid, Format, Val(iVal));
                end

                % Print separator if not last entry.
                if iVal < NVal && not(isnan(Val(iVal+1)))
                    fprintf(fid, '%s', Sep);
                end
            end

            % Print trailing string if present.
            if not(isempty(Trail))
                fprintf(fid, ' %s\n', Trail);

            % Otherwise end line.
            else
                fprintf(fid, '\n');
            end

        end

        function P = VoxelSamples(edge, sp, ax, h, r)
        % Generate points inside a cylinder with the given parameters:
        % start point <sp>, axis direction <ax>, length <h>, and radius
        % <r>. The <edge> parameter is used to chose the number of samples.

            % Average distance between points.
            PointDist = edge/2;

            % Number vertex rings inside cylinder.
            NRing = ceil(r/PointDist);

            % Radii of each ring. Note that zero is computed separately.
            RingRadius = linspace(0,1,NRing + 1);
            RingRadius = RingRadius(2:end);

            % Number of vertices on ring loop.
            NVertex = max(ceil( pi./asin(PointDist./(2*r*RingRadius)) ),3);

            % Number of vertical layers.
            NLayer = ceil(h/PointDist) + 1;

            % Heights at which ring layers are distributed at.
            z = linspace(0,h,NLayer);

            % Total number of points.
            NPoint = sum(NVertex)*NLayer;

            % Vertices.
            P = zeros(NPoint,3);

            iPoint = 1;

            for iRing = 1:NRing

                NPointRing = NVertex(iRing)*NLayer;

                ri = RingRadius(iRing);

                Ang = linspace(0,2*pi,NVertex(iRing) + 1);
                Ang = Ang(1:end-1);

                [x,y] = pol2cart(Ang,ri*r);

                X = repmat(x,NLayer,1);
                Y = repmat(y,NLayer,1);
                Z = repmat(z(:),1,NVertex(iRing));

                P(iPoint:iPoint+NPointRing-1,:) = [X(:), Y(:), Z(:)];

                iPoint = iPoint + NPointRing;

            end

            % Add points on axis.
            P = vertcat(P,[zeros(NLayer,2), z(:)]);

            % Coordinate system matrix.

            if ax(3) == 1
                R = eye(3);
            elseif ax(3) == -1
                R = [-1 0 0; 0 1 0; 0 0 -1];
            else
                axx = cross(ax,[0 0 1]);
                axx = axx./norm(axx);

                axy = cross(axx,ax);

                R = [axx; axy; ax];
            end


            % Rotate to cylinder coordinate system and translate.
            P = bsxfun(@plus,P*R,sp);
            
            % Debug: plot resulting points.
            if QSMBCylindrical.debug
                plott(P,'b.');
                hold on;
                axis equal;
            end

        end

    end

end

function R = rotation_matrix(u,k)
% Compute rotation around axis <u> by <k> radians.

    R = zeros(3,3);
    c = cos(k); 
    s = sin(k);
    R(1,:) = [u(1)^2+(1-u(1)^2)*c,    ...
              u(1)*u(2)*(1-c)-u(3)*s, ...
              u(1)*u(3)*(1-c)+u(2)*s];
    %-
    R(2,:) = [u(1)*u(2)*(1-c)+u(3)*s, ...
              u(2)^2+(1-u(2)^2)*c,    ...
              u(2)*u(3)*(1-c)-u(1)*s];
    %-
    R(3,:) = [u(1)*u(3)*(1-c)-u(2)*s, ...
              u(2)*u(3)*(1-c)+u(1)*s, ...
              u(3)^2+(1-u(3)^2)*c];
    %-

end

function TwigParam = default_twig_param_dist(sp, ax, h, r, is_last, len)
% Parameters of the twigs. Row = single twig.
% Columns:
%   1: position on axis on the inverval [0,1].
%      Equals 1 if on head.
%   2: position on radial axis on the inverval [0,1]. 
%      Equals 1 if on side, less than 1 when on head.
%   3: rotation around axis. From -Pi to Pi.
%   4: twig elevation. -Pi/2 backwards, Pi/2 towards tip.
%   5: twig azimuth. -Pi/2 left, Pi/2 right.
%   6: twig length.

    % Limits for even distributions.
    Limits = [0 1; 0 1; -pi pi; -pi/2 pi/2; -pi/2  pi/2; len];
    %Limits = [0 1; 0 1; -pi pi;     0 pi/2; -pi/4 -pi/4; len];

    NTwig = size(ax,1);
    NTwig = max(NTwig,size(sp,1));
    
    TwigParam = zeros(NTwig,6);
    
    for iParam = 1:6
        
        TwigParam(:,iParam) = Limits(iParam,1) ...
                            + (Limits(iParam,2) - Limits(iParam,1)) ...
                            * rand(NTwig,1);
        %-
        
    end
    
    % Ratio of area of circle and complete area.
    ratio = r./(r + 2*h);

    % Leafs that are connected to cylinder envelope and not the end 
    % circle.
    onhead = rand(NTwig,1) < ratio(:) & is_last(:);

    % Set axial translation for head-connected leaves to full.
    TwigParam(onhead,1) = 1;
    % Set radial translation for side-connected leaves to full.
    TwigParam(not(onhead),2) = 1;
    
    

end

function plot_axes(config, index, Origin, TwigUp,TwigSide,TwigFront)
% Function to plot coordinate system axis vectors. Used for debugging.

    colors = eye(3);
    E = [TwigUp; TwigSide; TwigFront];

    subplot(config(1),config(2),index);
    
    for i = 1:3
        plott([Origin; Origin+E(i,:)],'-','Color',colors(i,:));
        hold on;
    end
    
    hold off;
    axis equal;
    legend({'up','side','front'});

end