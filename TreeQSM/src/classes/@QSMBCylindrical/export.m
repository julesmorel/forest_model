% Export cylinder geometry (and topology) in various formats.
%
% Inputs:
%
% Format            String descriptor of the export format. The value is 
%                   case-insensitive. Possible values are:
%
% 'OBJ'             Wavefront OBJ format. Cylinders are converted to 
%                   vertices and faces prior to export.
%
% 'PLY'             Polygon file format. Cylinders are converted to 
%                   vertices and faces prior to export.
%
% 'TXT' or 'TEXT'   Proprietary format that contains cylinder, branch and
%                   tree-level geometry and topology data.
%
% 'Blender'         Proprietary format that is designed to for the
%                   qsm-blender-addon. Contains cylinder geometry data.
%
% File              File name to which the exported data is written to.
%
%
% Optional inputs given in name-value pairs. Possible values are:
%
% 'Precision'       Scalar or two-element vector desribing float
%                   printing format. A scalar assigns the number
%                   of included decimals. The first element of a
%                   vector does the same, while the second element
%                   sets the total width of the printed float.
%                   Default value is 4.
%
% 'FaceCount'       Set the integer count or faces to compute 
%                   on each cylinder. Can be either a single
%                   integer or a two-element vector of 
%                   integers, in which case the face count is 
%                   interpolated linearly based on the relative
%                   radii, i.e., smaller cylinders have less
%                   faces. Default value is 5.
%
% 'Quads'           Boolean. If true use quads instead of
%                   triangles on cylinder envelope faces.
%                   By default TRUE. Can be defined without
%                   value pair, then value is interpreted as 
%                   TRUE.
%
% 'Closed'          Boolean. If true a triangle fan of faces is 
%                   generated to "close" the cylinder top and
%                   bottom. By default FALSE.  Can be defined 
%                   without value pair, then value is 
%                   interpreted as TRUE.
%
% 'TriangleStem'    Boolean. If TRUE triangulated stem data is
%                   used rather than cylinder geometry when
%                   exporting the model. By default FALSE.
%                   Can be defined without value pair, then
%                   value is interpreted as TRUE. Throws a
%                   warning if model has no triangulated stem.
%
% 'Origin'          Three element vector that is used to move
%                   the base of the tree model. Empty by
%                   default. Use QSM.cylinder_start_point(1,:)
%                   to move base to [0,0,0].
%
% 'Filter'          Boolean vertor with one element per
%                   cylinder or a string descriptor passed to 
%                   the QSMBCylindrical.get_cylinder_set() 
%                   method. Only export cylinders that are set
%                   to TRUE.
%
% 'Color'           Either a [N x 1] vector or a [N x 3] matrix
%                   of color values, with N matching the
%                   cylinder count. Describes the color values
%                   for each cylinder. Value can also be a
%                   single three-element color vector to be
%                   used for all cylinders. Only applicable for
%                   PLY and Blender file formats.
%
% 'ScaleColor'      Boolean. When TRUE, color values given by the
%                   'Color' attribute are scaled to fill the full
%                   integer interval [0, 255].
%
% 'ColorMatrix'     [N x 3] matrix of true colors to be used to
%                   color cylinders according to discrete
%                   properties such as branch order or branch
%                   index. The constant property color_matrix
%                   contains the default color matrix. A color
%                   matrix can have any number of rows as it is
%                   indexed with looping.
%
% 'ColorSource'     String descriptor of property to use for
%                   coloring vertices, faces or cylinders. Use
%                   as an alternative for the 'Color' attribute.
%                   Only applicable for PLY and Blender file 
%                   formats. Discrete color sources use the colors 
%                   from the bject color_matrix property or an 
%                   override provided by the 'ColorMatrix'
%                   attribute. Continuous color sources use the
%                   current system-wide COLORMAP.
%                   QSMBCylindrical.get_color_distribution()
%                   method is used to compute the color 
%                   distribution. Possible values are:
%
%       'None'                      No color based on source. 
%                                   This is the default value.
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
%       'VertexHeight'              Height of vertices.
%
%       'VertexRadius'              Horizontal distance of vertex
%                                   from tree base.
%
%       'FaceHeight'                Mean height of face vertices.
%
%
% 'Comments'        Boolean. Extra comments are printed in the 
%                   resulting file when TRUE. By default TRUE.  
%                   Can be defined without value pair, then value
%                   is interpreted as TRUE.
%
% 'Header'          Cell-array of string that are printed at the 
%                   top of the export file, one element (row) at 
%                   a time. The user is resposible for including
%                   the format specific comment syntax, when
%                   required.
%
% 'ExtraData'       Numeric array or a cell-array of extra data
%                   to be included in the exported file. Applicable
%                   with the 'TXT' and 'Blender' formats. If the
%                   value is numeric, the row count should match
%                   the object cylinder count, as the values on
%                   each row are printed after the default cylinder
%                   attributes. Can also be given as a cell-array
%                   when using the 'TXT' format. The cell-array 
%                   should have at most three elements, with the
%                   first containing numeric cylinder-level data.
%                   second numeric branch-level data and the third
%                   containing numeric or a cell array of tree-level
%                   data.
%
%
% Examples:
%
% % Basic export in various formats.
% QSM.export('PLY','qms-export.ply');
% QSM.export('OBJ','qms-export.obj');
% QSM.export('TXT','qms-export.txt');
% QSM.export('Blender','qms-export.txt');
%
% % Set precision and closed cylinder ends.
% QSM.export('PLY','qms-export.ply','Precision',[5,10],'Closed');
%
% % Export only stem.
% QSM.export('PLY','qms-export.ply','Filter','Stem');
%
% % Set all cylinder colors as red.
% QSM.export('PLY','qms-export.ply','Color',[255,0,0]);
%
% % Color by source: branch order.
% QSM.export('PLY','qms-export.ply','ColorSource','Order');



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

function export(ob, Format, File, varargin)

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
    % Assign color to vertices or cylinders when true.
    % Activated with the 'Color' argument.
    Param.FColor = false;
    % Scale color values to interval [0, 255].
    Param.FScaleColor = false;
    % Write extra columns to text files when true.
    % Activated with the 'Param.ExtraData' argument.
    Param.FExtraData = false;
    % Only select subset of cylinders.
    Param.FFilter = false;
    % Include comments to output files when true.
    Param.FComment = true;

    % Default values.

    % Cell array of string lines written at appropriate part of
    % output file.
    Param.Header = {};

    % Float precision for printing floating point numbers.
    Param.FloatPrec = 4;

    % Number of side faces computed for a single cylinder.
    Param.FaceCount = 5;

    % Numeric extra data to be output after default paramenters.
    Param.ExtraData = [];

    % Color data for cylinders.
    Param.CylColor = [];
    
    % Logical vector for selecting a subset of cylinders.
    Param.ICyl = [];

    % Vector for moving the tree model start point.
    % By default no tranlation.
    Param.Origin = [0, 0, 0];

    % Use default color matrix of object.
    Param.ColorMatrix = ob.color_matrix;

    % Default color source. Color by vertex height.
    Param.ColorSource = '';

    % Flag: use trianglulated stem vertices and faces
    % and exclude respective cylinders.
    Param.FTriStem = false;

    % List of acceptable arguments for current function.
    Accepted = { ...
        'header','precision','facecount','closed','quads', ...
        'comments','scalecolor','filter','color','origin', ...
        'extradata', 'colorsource', 'colormatrix', 'trianglestem' ...
    };

    % Parse optional input arguments.
    Param = ob.parse_arguments(Param, varargin, Accepted);

    % Color matrix scaling.
    if Param.FScaleColor

        % If filtering is enabled, only select filtered 
        % cylinder color values for extreme value compututation.
        if Param.FFilter
            SelCylColor = Param.CylColor(Param.ICyl);
        else
            SelCylColor = Param.CylColor;
        end

        % Find extreme values.
        cmin = min(SelCylColor(:));
        cmax = max(SelCylColor(:));

        % Project to [0 1].
        cscale = (Param.CylColor - cmin)./(cmax - cmin);

        % Project to integers [0, 255].
        Param.CylColor = floor(cscale*255);
    end

    % Define float print format.

    % If the precision argument has two elements, the first
    % element is the precision and the second is the total 
    % field width.
    if length(Param.FloatPrec) > 1
        FloatFormat = sprintf(...
            '%%%d.%dg', ...
            Param.FloatPrec(2), ...
            Param.FloatPrec(1) ...
        );

    % Otherwise only set the precision.
    else
        FloatFormat = sprintf('%%.%dg', Param.FloatPrec);
    end

    % Define int print format.
    IntFormat = '%d';

    % Number format that can be updated.
    NumFormat = FloatFormat;

    % Detect extra data number formats.
    if Param.FExtraData

        % Cell array can have components for cylinders, branches and
        % tree-level data.
        if iscell(Param.ExtraData)

            % Column format collection. One element per extra data 
            % element.
            ExtraColFormat = cell(size(Param.ExtraData));

            % Each component.
            for iData = 1:numel(Param.ExtraData)

                if iData == 3 && iscell(Param.ExtraData{3})
                    % The component-level format collection is a cell
                    % array with one element per cell.
                    ExtraColFormat{iData} = cell( ...
                        1, ...
                        numel(Param.ExtraData{3}) ...
                    );

                    % Check each element.
                    for iCol = 1:numel(Param.ExtraData{3})

                        % Detect integers.
                        if all(mod(Param.ExtraData{3}{iCol},1) == 0)

                            ExtraColFormat{iData}{iCol} = IntFormat;

                        % Otherwise use float.
                        else

                            ExtraColFormat{iData}{iCol} = FloatFormat;
                        end

                    end
                else
                    % The component-level format collection is a cell
                    % array with one element per column.
                    ExtraColFormat{iData} = cell( ...
                        1, ...
                        size(Param.ExtraData{iData},2) ...
                    );

                    % Check each column.
                    for iCol = 1:size(Param.ExtraData{iData},2)

                        % Detect integers.
                        if all(mod(Param.ExtraData{iData}(:,iCol),1) == 0)

                            ExtraColFormat{iData}{iCol} = IntFormat;

                        % Otherwise use float.
                        else

                            ExtraColFormat{iData}{iCol} = FloatFormat;
                        end

                    end
                end

                
            end

        % Otherwise extra data is a single matrix and contain extra data
        % only for cylinders.
        else

            % Format container can be a single-level cell array with
            % one elemnt per data column.
            ExtraColFormat = cell(1,size(Param.ExtraData,2));

            % Check each column.
            for iCol = 1:size(Param.ExtraData,2)

                % Detect integers.
                if all(mod(Param.ExtraData(:,iCol),1) == 0)

                    ExtraColFormat{iCol} = IntFormat;

                % Otherwise use float.
                else
                    
                    ExtraColFormat{iCol} = FloatFormat;
                end

            end
        end

    end

    % Open file stream.
    fid = fopen(File,'w');

    % Check that file opened properly.
    assert(fid ~= -1, ['Unable to open target file: ''' File '''']);

    % Select type of export.
    switch lower(Format)

        % Polygon File Format.
        case 'ply'

            % Compute cylinder vertices and faces. Filtering can
            % limit inlcuded cylinders. Face count, side face 
            % shape and top and bottom triangle fans can be 
            % configured.
            [Vert, Faces, JVertCyl] = ob.compute_geometry(...
                Param.FaceCount, ...
                Param.ICyl, ...
                Param.FQuads, ...
                Param.FClosed ...
            );

            % If location is offset, move vertices.
            if any(Param.Origin)
                Vert = bsxfun(@minus, Vert, Param.Origin);
            end

            % If manual color is off and some color source is active,
            % compute color based on string descriptor.
            if not(Param.FColor) && ...
                not(isempty(Param.ColorSource)) && ...
                not(strcmp(Param.ColorSource,'none'))

                % Try to compute color distribution. May throw error if
                % color source string is unknown.
                ColorDist = ob.get_color_distribution(...
                    Vert, ...
                    Faces, ...
                    JVertCyl, ...
                    Param.ColorSource, ...
                    Param.ColorMatrix ...
                );

                % Convert indexed to true colors.
                if size(ColorDist,2) == 1

                    % Get current colormap.
                    CMap = colormap();

                    % Number of colors in colormap.
                    NColorMap = size(CMap,1);

                    % Number of colors in color distribution.
                    NColorDist = max(ColorDist) - min(ColorDist);

                    % If color count do not match, scaling is required.
                    if NColorMap ~= NColorDist

                        % Map color distribution to available colors.
                        ColorDist = floor( ...
                            (ColorDist - min(ColorDist))/(NColorDist) ...
                            * (NColorMap-1) ...
                        ) + 1;
                    end

                    % Index color map to receive true colors.
                    ColorDist = CMap(ColorDist,:);
                end

                % Scale to integer values in interval [0, 255].
                ColorDist = round(ColorDist*255);

                % Update cylinder color matrix.
                Param.CylColor = ColorDist;
                % Set color flag on.
                Param.FColor = true;
            end

            % If usage of triangulated stem is requested, convert
            % cylinder geometry to triangulated geometry and update
            % possible color data.
            if Param.FTriStem

                [Vert, Faces, JVertCyl, ColorDist] = ...
                    ob.convert_stem_vert( ...
                        Vert, ...
                        Faces, ...
                        JVertCyl, ...
                        Param.CylColor ...
                );

                Param.CylColor = ColorDist;
            end
        
            % Number of vertices.
            NVertex = size(Vert,1);
            % Number of faces.
            NFace = size(Faces,1);

            % PLY is zero-indexed.
            Faces = Faces - 1;

            % Number of corner vertices on each face. Can vary if
            % using quads and having closed cylinders.
            FaceSizes = sum(~isnan(Faces),2);

            % In PLY face lines start with the face vertex count.
            Faces = cat(2, FaceSizes, Faces);

            % Print default header.
            fprintf(fid, '%s\n','ply');
            fprintf(fid, '%s\n','format ascii 1.0');

            % Print custom header.
            for iLine = 1:numel(Param.Header)
                fprintf(fid,'%s\n',Param.Header{iLine});
            end

            % Print elements and element properties.

            % Vertex count and format.
            fprintf(fid, 'element vertex %d\n',NVertex);
            fprintf(fid, '%s\n','property float x');
            fprintf(fid, '%s\n','property float y');
            fprintf(fid, '%s\n','property float z');

            % If color information should be included.
            if Param.FColor

                % Vertex lines will have three additional elements.
                fprintf(fid, '%s\n','property uchar red');
                fprintf(fid, '%s\n','property uchar green');
                fprintf(fid, '%s\n','property uchar blue');

                % Add color data to vertex matrix.
                if size(Param.CylColor,1) == size(Vert,1)
                    Vert = cat(2, Vert, Param.CylColor);
                else
                    Vert = cat( ...
                        2, ...
                        Vert, ...
                        Param.CylColor(JVertCyl,:) ...
                    );
                end

                % Number format for vertex and color information.
                % Color data is in integer format.
                NumFormat = { ...
                    FloatFormat, ...
                    FloatFormat, ...
                    FloatFormat, ...
                    IntFormat, ...
                    IntFormat, ...
                    IntFormat ...
                };

            end

            % Face count and format.
            fprintf(fid, 'element face %d\n',NFace);
            fprintf(fid, '%s\n','property list uchar int vertex_indices');
            fprintf(fid, '%s\n','end_header');

            % Print vertex data.
            for iVertex = 1:NVertex

                % Trailing line comment.
                if iVertex == 1 && Param.FComment
                    LineComment = '{ start of vertex list }';
                else
                    LineComment = '';
                end

                % Print vertex line.
                QSMBCylindrical.print_number_line(...
                    fid, ...
                    Vert(iVertex,:),  ...
                    NumFormat, ...      % Number format.
                    '', ...             % Line lead.
                    LineComment ...     % Line trail.
                );
            end

            % Print face data.
            for iFace = 1:NFace

                % Trailing line comment.
                if iFace == 1 && Param.FComment
                    LineComment = '{ start of face list }';
                else
                    LineComment = '';
                end

                % Print face line. Numbers are integers.
                QSMBCylindrical.print_number_line(...
                    fid, ...
                    Faces(iFace,:),  ...
                    IntFormat, ...  % Number format.
                    '', ...         % Line lead.
                    LineComment ... % Line trail.
                );
            end

        % Wavefront OBJ file.
        case 'obj'

            % Compute cylinder vertices and faces. Filtering can
            % limit inlcuded cylinders. Face count, side face 
            % shape and top and bottom triangle fans can be 
            % configured.
            [Vert, Faces, JVertCyl] = ob.compute_geometry(...
                Param.FaceCount, ...
                Param.ICyl, ...
                Param.FQuads, ...
                Param.FClosed ...
            );

            % If usage of triangulated stem is requested, convert
            % cylinder geometry to triangulated geometry and update
            % possible color data.
            if Param.FTriStem

                % Vertex-cylinder mapping is excluded as output because
                % not used after this point.
                [Vert, Faces] = ...
                    ob.convert_stem_vert( ...
                        Vert, ...
                        Faces, ...
                        JVertCyl ...
                );
            end
        
            % Number of vertices.
        	NVertex = size(Vert,1);

            % Number of faces.
            NFace = size(Faces,1);

            % If location is offset, move vertices.
            if any(Param.Origin)
                Vert = bsxfun(@minus, Vert, Param.Origin);
            end

            % Print custom header.
            for iLine = 1:numel(Param.Header)
                fprintf(fid,'%s\n',Param.Header{iLine});
            end

            % Print vertex data.
            for iVertex = 1:NVertex

                % Print vertex line. Numbers are floats.
                QSMBCylindrical.print_number_line(...
                    fid, ...
                    Vert(iVertex,:),  ...
                    NumFormat, ...      % Number format.
                    'v ' ...            % Line lead.
                );
            end

            % Print face data.
            for iFace = 1:NFace

                % Print face line. Numbers are integers.
                QSMBCylindrical.print_number_line(...
                    fid, ...
                    Faces(iFace,:),  ...
                    IntFormat, ...  % Number format.
                    'f ' ...        % Line lead.
                );
            end

        % Custom text file format containing cylinder, branch and tree 
        % data.
        case {'txt', 'text'}

            % Print custom header.
            for iLine = 1:numel(Param.Header)
                fprintf(fid,'%s\n',Param.Header{iLine});
            end

            % Section header as comment.
            if Param.FComment
                fprintf(fid,'# Cylinder data\n');
            end

            % Go through cylinders. Print radius, length, start point,
            % axis direction, parent index, extension index, branch index,
            % branch order, index in branch and added flag. Use
            % appropriate number format.
            for iCyl = 1:ob.block_count

                fprintf(...
                    fid, ...
                    FloatFormat, ...
                    ob.cylinder_radius(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' FloatFormat], ...
                    ob.cylinder_length(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' FloatFormat], ...
                    ob.cylinder_start_point(iCyl,:) ...
                );
                fprintf(...
                    fid, ...
                    [' ' FloatFormat], ...
                    ob.cylinder_axis(iCyl,:) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_parent(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_extension(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_branch_index(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_branch_order(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_index_in_branch(iCyl) ...
                );
                fprintf(...
                    fid, ...
                    [' ' IntFormat], ...
                    ob.cylinder_added(iCyl) ...
                );

                % If extra data given, print after standard parameters on
                % the same line.
                if Param.FExtraData

                    % Extra data can be given as a cell array. Then the
                    % first element is expected to be extra data for 
                    % the cylinder lines.
                    if iscell(Param.ExtraData)
                        CylParam.ExtraData = Param.ExtraData{1};
                        ColFormat = ExtraColFormat{1};

                    % Otherwise the entire variable is used as extra
                    % cylinder data.
                    else
                        CylParam.ExtraData = Param.ExtraData;
                        ColFormat = ExtraColFormat;
                    end

                    % Go through the columns and print each additional 
                    % data point.
                    for iCol = 1:size(CylParam.ExtraData,2)

                        fprintf(...
                            fid, ...
                            [' ' ColFormat{iCol}], ...
                            CylParam.ExtraData(iCyl,iCol) ...
                        );

                    end
                end

                fprintf(fid, '\n');

            end

            % Section header as comment.
            if Param.FComment
                fprintf(fid,'# Branch data\n');
            end

            % For each cylinder print order, parent index, volume, length,
            % angle, and height. Use appropriate number format.
            for iBranch = 1:ob.branch_count

                fprintf(fid, ...
                    IntFormat, ...
                    ob.branch_order(iBranch) ...
                );
                fprintf(fid, ...
                    [' ' IntFormat], ...
                    ob.branch_parent(iBranch) ...
                );
                fprintf(fid, ...
                    [' ' FloatFormat], ...
                    ob.branch_volume(iBranch) ...
                );
                fprintf(fid, ...
                    [' ' FloatFormat], ...
                    ob.branch_length(iBranch) ...
                );
                fprintf(fid, ...
                    [' ' FloatFormat], ...
                    ob.branch_angle(iBranch) ...
                );
                fprintf(fid, ...
                    [' ' FloatFormat], ...
                    ob.branch_height(iBranch) ...
                );

                % If extra data given, print one column at a time.
                if Param.FExtraData && ...
                    iscell(Param.ExtraData) && ...
                    numel(Param.ExtraData) > 1

                    % Get branch-level extra data.
                    BranchParam.ExtraData = Param.ExtraData{2};
                    % Column number formats.
                    ColFormat = ExtraColFormat{2};

                    % Print each column.
                    for iCol = 1:size(BranchParam.ExtraData,2)

                        fprintf(...
                            fid, ...
                            [' ' ColFormat{iCol}], ...
                            BranchParam.ExtraData(iBranch,iCol) ...
                        );

                    end
                end

                % End current branch line.
                fprintf(fid, '\n');

            end

            % Print section header as comment.
            if Param.FComment
                fprintf(fid,'# Tree data\n');
            end

            % TotalVolume
            % TrunkVolume
            % BranchVolume
            % TreeHeight
            % TrunkLength
            % BranchLength
            % NumberBranches    Total number of branches
            % MaxBranchOrder 
            % TotalArea 
            % DBHqsm        From the cylinder of the QSM at the right heigth
            % DBHcyl        From the cylinder fitted to the section 1.1-1.5m
            % location      (x,y,z)-coordinates of the base of the tree
            % StemTaper     Stem taper function/curve from the QSM
            % VolumeCylDiam     Distribution of the total volume in diameter classes
            % LengthCylDiam     Distribution of the total length in diameter classes
            % VolumeBranchOrder     Branch volume per branching order
            % LengthBranchOrder     Branch length per branching order
            % NumberBranchOrder     Number of branches per branching order

            % Inline function to handle optional comment string printing.
            % If comments are not included, the function returns an empty
            % string.
            if Param.FComment
                c = @(str) str;
            else
                c = @(str) '';
            end

            % Print tree-level properties.

            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.volume(), ...
                FloatFormat, '', ...
                c('# Total volume') ...
            );

            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.volume('stem'), ...
                FloatFormat, '', ...
                c('# Stem volume') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.volume('branch'), ...
                FloatFormat, '', ...
                c('# Branch volume') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.height(),  ...
                FloatFormat, '', ...
                c('# Tree height') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.length('stem'), ...
                FloatFormat, '', ...
                c('# Stem length') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.length('branch'), ...
                FloatFormat, '', ...
                c('# Branch length') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.branch_count, ...
                IntFormat, '', ...
                c('# Branch count') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                max(ob.branch_order), ...
                IntFormat, '', ...
                c('# Maximum branch order') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.area(), ...
                FloatFormat, '', ...
                c('# Total area') ...
            );
            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.dbh(), ...
                FloatFormat, '', ...
                c('# DBH from QSM') ...
            );

            % Check if custom property stored.
            if ob.has_property('DBH CYL')
                DBH = ob.get_property('DBH CYL');
            else
                DBH = 0;
            end
                
            QSMBCylindrical.print_number_line(...
                fid, ...
                DBH, ...
                FloatFormat, '', ...
                c('# Cylinder fitting DBH') ...
            );

            QSMBCylindrical.print_number_line(...
                fid, ...
                ob.cylinder_start_point(1,:), ...
                FloatFormat, '', ...
                c('# Tree start point') ...
            );

            % Extra tree-level data.
            if Param.FExtraData && ...
                iscell(Param.ExtraData) && ...
                numel(Param.ExtraData) > 2

                % Get extra data.
                TreeParam.ExtraData = Param.ExtraData{3};

                % Column formats.
                ColFormat = ExtraColFormat{3};

                if iscell(TreeParam.ExtraData)

                    % Print each cell.
                    for iCol = 1:numel(TreeParam.ExtraData)

                        % Format for cell content.
                        CellFormat = repmat( ...
                            [' ' ColFormat{iCol}], ...
                            1, ...
                            length(TreeParam.ExtraData{iCol}) ...
                        );

                        % Trim leading space.
                        CellFormat = CellFormat(2:end);

                        % Print scalar or vector.
                        fprintf(...
                            fid, ...
                            [CellFormat '\n'], ...
                            TreeParam.ExtraData{iCol} ...
                        );

                    end
                else

                    % Print each column.
                    for iCol = 1:size(TreeParam.ExtraData,2)

                        fprintf(...
                            fid, ...
                            [ColFormat{iCol} '\n'], ...
                            TreeParam.ExtraData(iCol) ...
                        );

                    end
                end
            end

        % Format for importing data for blender visualizations. The text
        % file will include the minimum amount of cylinder-level data
        % for rendering the model.
        case 'blender'

            % If manual color is off and some color source is active,
            % compute color based on string descriptor.
            if not(Param.FColor) && ...
                not(isempty(Param.ColorSource)) && ...
                not(strcmp(Param.ColorSource,'none'))

                assert(not(any(ismember(Param.ColorSource,{ ...
                    'vertexradius', ...
                    'vertexheight', ...
                    'faceheight' ...
                    }))), ...
                    [ ...
                        'Vertex and face level color sources are not ' ...
                        'available for the ''' ...
                        Format ...
                        '''.' ...
                    ] ...
                );

                % Try to compute color distribution. May throw error if
                % color source string is unknown.
                ColorDist = ob.get_color_distribution(...
                    [], ...
                    [], ...
                    [], ...
                    Param.ColorSource, ...
                    Param.ColorMatrix ...
                );

                % Convert indexed to true colors.
                if size(ColorDist,2) == 1

                    % Get current colormap.
                    CMap = colormap();

                    % Number of colors in colormap.
                    NColorMap = size(CMap,1);

                    % Number of colors in color distribution.
                    NColorDist = max(ColorDist) - min(ColorDist);

                    % If color count do not match, scaling is required.
                    if NColorMap ~= NColorDist

                        % Map color distribution to available colors.
                        ColorDist = floor( ...
                            (ColorDist - min(ColorDist))/(NColorDist) ...
                            * (NColorMap-1) ...
                        ) + 1;
                    end

                    % Index color map to receive true colors.
                    ColorDist = CMap(ColorDist,:);
                end

                % Update cylinder color matrix.
                Param.CylColor = ColorDist;
                % Set color flag on.
                Param.FColor = true;
            end

            % Print custom header.
            for iLine = 1:numel(Param.Header)
                fprintf(fid,'%s\n',Param.Header{iLine});
            end

            % Print properties of each cylinder.
            for iCyl = 1:ob.block_count

                if Param.FFilter && not(Param.ICyl(iCyl))
                    continue;
                end

                % The following properties are exported.
                % - Branch ID
                % - Starting point
                % - Axis direction
                % - Length
                % - Radius
                fprintf(fid,...
                    [IntFormat repmat([' ' FloatFormat],1,8)], ...
                    ob.cylinder_branch_index(iCyl), ...
                    ob.cylinder_start_point(iCyl,:) - Param.Origin, ...
                    ob.cylinder_axis(iCyl,:), ...
                    ob.cylinder_length(iCyl), ...
                    ob.cylinder_radius(iCyl) ...
                );

                % Print color columns if present.
                if Param.FColor

                    for iCol = 1:size(Param.CylColor,2)
                        
                        % Select number format.
                        if all(mod(Param.CylColor(:,iCol),1) == 0)
                            NumFormat = IntFormat;
                        else
                            NumFormat = FloatFormat;
                        end

                        fprintf( ...
                            fid, ...
                            [' ' NumFormat], ...
                            Param.CylColor(iCyl,iCol) ...
                        );
                    end
                end

                % Print extra columns if present.
                if Param.FExtraData

                    % Print each column.
                    for iCol = 1:size(Param.ExtraData,2)
                        
                        % Select number format.
                        if all(mod(Param.ExtraData(:,iCol),1) == 0)
                            NumFormat = IntFormat;
                        else
                            NumFormat = FloatFormat;
                        end

                        fprintf( ...
                            fid, ...
                            [' ' NumFormat], ...
                            Param.ExtraData(iCyl,iCol) ...
                        );
                    end
                end

                % Print line end.
                fprintf(fid,'\n');
            end

        % Unknown export type.
        otherwise
            error(['Unknown argument: ''' Format '''']);
    end

    % Close stream.
    fclose(fid);

end