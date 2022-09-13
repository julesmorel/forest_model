% Compute an indexed color distribution based on various vertex, face,
% cylinder of branch-level attributes. Returns a vector distribution of 
% color indices. The length of the vector matches either the vertex count
% or the face count, depending on the selected color source.
%
% Inputs:
%
% Vertices          Nx3 matrix of vertex (x,y,z)-coordinates. One vertex
%                   per row.
% Faces             MxP matrix of vertex indices. One face per row. Rows
%                   can include NaNs if face vertex count varies.
% JVertCyl          Mapping from vertex index to cylinder index.
% ColorSource       String descriptor of the distribution to compute for
%                   coloring.
%
% Available color sources are:
%
% 'CylinderLength'              Length of a cylinder.
%
% 'CylinderIsLast'              Highlight last cylinders in each branch.
%
% 'CylinderAdded'               Highlight cylinders that were added after
%                               reconstruction.
%
% 'CylinderRadius'              Cylinder radius.
%
% 'CylinderOrder',
% 'CylinderBranchOrder',
% 'BranchOrder',
% 'Order'                       Branch order of cylinder.
%
% 'CylinderBranchIndex',
% 'CylinderIndex',
% 'BranchIndex'                 Separate branches with color.
%
% 'BranchLength'                Longer branches have higher index.
%
% 'BranchAngle'                 Branching angle at branch base.
%
% 'BranchHeight'                Mean height of branch.
%
% 'VertexHeight'                Height of vertices.
%
% 'VertexRadius'                Horizontal distance of vertex from tree 
%                               base.
%
% 'FaceHeight'                  Mean height of face vertices.



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

function Dist = get_color_distribution(ob, ...
	Vertices, ...
	Faces, ...
	JVertCyl, ...
	ColorSource, ...
    ColorMatrix ...
)

    %NVert = size(Vertices,1);
    NFace = size(Faces,1);

    % If called without vertex and face data.
    if isempty(JVertCyl)
        JVertCyl = 1:ob.block_count;
    end

    % Store original colorsource string for error printing.
    OrigColorSource = ColorSource;

    FDiscrete = false;

    % Change color source string for some cases to unify computations.
    switch lower(ColorSource)

        case 'cylinderlength'
            ColorSource = 'cylinder_length';
        case {'cylinderislast', 'cylinderlast'}
            ColorSource = 'cylinder_is_last';
            FDiscrete = true;
        case 'cylinderradius'
            ColorSource = 'cylinder_radius';
        case {'cylinderorder', 'cylinderbranchorder','branchorder','order'}
            ColorSource = 'cylinder_branch_order';
            FDiscrete = true;
        case {'cylinderindex', 'cylinderbranchindex','branchindex'}
            ColorSource = 'cylinder_branch_index';
            FDiscrete = true;
        case 'cylinderadded'
            ColorSource = 'cylinder_added';
            FDiscrete = true;
        case 'branchlength'
            ColorSource = 'branch_length';
        case 'branchangle'
            ColorSource = 'branch_angle';
        case 'branchheight'
            ColorSource = 'branch_height';
    end

    % Distribution computation.
    switch lower(ColorSource)

        % Vertex horizontal distance from tree base.
        case 'vertexradius'
            Dist = bsxfun( ...
                @minus, ...
                Vertices(:,1:2), ...
                ob.cylinder_start_point(1,1:2) ...
            );

            Dist = sqrt(sum(Dist.^2,2));

        % Vertex height.
        case 'vertexheight'
            Dist = Vertices(:,3);

        % Mean height of face vertices.
        case 'faceheight'

            % Get vertex heights.
            VertHeight = Vertices(:,3);

            % Initialize distibution.
            Dist = zeros(NFace,1);

            % Go through face matrix.
            for iFace = 1:NFace

                % Compute the mean value of face vertices with non-nan
                % indices.
                Dist(iFace) = mean( ...
                    VertHeight(Faces( ...
                        iFace, ...
                        not(isnan(Faces(iFace,:))) ...
                    )), ...
                    1 ...
                );

            end

        % Height of mid point of cylinder.
        case 'cylinderheight'
            Dist = ob.cylinder_mid_point(:,3);
            Dist = Dist(JVertCyl);

        % Direct cylinder-level attributes.
        case { ...
            'cylinder_length', ...
            'cylinder_is_last', ...
            'cylinder_radius', ...
            'cylinder_branch_order', ...
            'cylinder_added', ...
            'cylinder_branch_index' ...
        }
            % Get property value.
            Dist = ob.(ColorSource);

            % Move from cylinder-level to vertex-level.
            Dist = Dist(JVertCyl);

        % Direct branch-level attributes.
        case { ...
            'branch_length', ...
            'branch_angle', ...
            'branch_height' ...
        }

            % Get property value.
            Dist = ob.(ColorSource);

            % Move from branch- to cylinder-level.
            Dist = Dist(ob.cylinder_branch_index);
            % Move from cylinder to face-level.
            Dist = Dist(JVertCyl);

        case 'none'
            Dist = [];
            return;

        % Otherwise throw error of unknown source. Use original string.
        otherwise
            error(['Unknown color source: ''' OrigColorSource '''']);
    end

    % Convert distribution to floats for scaling.
    Dist = single(Dist);

    % Extreme values.
    mi = min(Dist);
    ma = max(Dist);

    % Linear scaling to interval [0, 1].
    Dist = (Dist - mi)./(ma - mi);

    % Scale to integer values in interval [0, 255].
    Dist = round(Dist*255);

    if FDiscrete
        JUnique = unique(Dist);

        ColorDist = zeros(size(Dist,1),3);

        NVal = length(JUnique);

        NColor = size(ColorMatrix,1);

        for iVal = 1:NVal
            JColor = mod(iVal - 1,NColor) + 1;

            JVal = Dist == JUnique(iVal);
            NVal = nnz(JVal);

            ColorDist(JVal,:) = repmat(ColorMatrix(JColor,:),NVal,1);
        end
        
        Dist = ColorDist;
    end

end