% QSMB - Quantitative Structure Model Blocks. Abstract class for holding
% quantitative structure information.

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

classdef QSMB

    properties
        
        % Function handle to a function that defines a twig parameter
        % distribution. The inputs are the parameters of the blocks.
        fun_twig_distribution;
        
        % Lower and upper limit for twig length.
        twig_length_limits = [];

    end

    properties (Access=protected)

        % Map structure to hold custom properties.
        custom_properties;

    end

    properties (SetAccess=protected)

        % Number of geometric primitives (blocks).
        block_count = 0;

        % Minimum and maximum point of axis-aligned bounding box.
        tree_limits = [];

        compressed;

    end


    methods
        
        % Constructor.
        function ob = QSMB()
            ob.custom_properties = containers.Map();
            ob.compressed = false;
        end

        function list_properties(ob)
        % List custom properties stored in the object.
        % Shows property names and values.

            % Number of stored properties.
            NProperty = ob.custom_properties.Count;

            % If no properties are set, exit.
            if NProperty == 0
                disp('No custom properties set.');
                return;
            end

            % Maximum length of property names.
            max_name = 0;
            % Maximum length of string presentation of property values.
            max_val = 0;

            % Property names cell array.
            KeySet = keys(ob.custom_properties);
            % Property value cell array.
            ValSet = values(ob.custom_properties);

            % Cell array to store string representations of values.
            StrValSet = cell(NProperty,1);

            % Go through all properties.
            for iProperty = 1:NProperty

                % Update max name length.
                max_name = max(max_name, length(KeySet{iProperty}));

                % Current property value.
                Val = ValSet{iProperty};

                % String representations of value.
                StrVal = '';

                % Numeric vector/matrix values and cell arrays.
                % Print dimensions and type.
                if (isnumeric(Val) && ~isscalar(Val)) || iscell(Val)

                    % Size matrix.
                    s = size(Val);

                    % Go through all dimensions.
                    for iDim = 1:length(s)

                        % Print length in dimensions to string.
                        StrVal = StrVal + sprintf('%d', s(iDim));

                        % Add cross.
                        if iDim < length(s)
                            StrVal = StrVal + 'x';
                        end
                    end

                    % Print property value class.
                    StrVal = cat(2,'[',StrVal, ' ', class(Val),']');

                % Numeric scalar value.
                elseif isnumeric(Val)

                    % Integer.
                    if mod(Val,1) == 0
                        StrVal = sprintf('%d', Val);
                    % Float.
                    else
                        StrVal = sprintf('%f', Val);
                    end

                % String property.
                elseif ischar(Val)

                    StrVal = ['''' Val ''''];

                % Just print property type for other types.
                else

                    StrVal = ['[' class(Val) ']'];

                end

                % Update maximum name length.
                max_val = max(max_val, length(StrVal));

                % Store string value.
                StrValSet{iProperty} = StrVal;

            end

            % Column separator string.
            ColBreak = '   ';
            
            % Column width should be at least the length of 
            % the column title.
            max_name = max(max_name, 4);
            max_val = max(max_name, 5);

            % Number of character to print on limiting lines.
            NDash = max_name+length(ColBreak)+max_val;

            fprintf('\n');
            % Top rule.
            fprintf(repmat('=',1,NDash));
            fprintf('\n');
            % Table header.
            fprintf(['%-' num2str(max_name) 's'], 'Name');
            fprintf(ColBreak);
            fprintf(['%' num2str(max_val) 's'], 'Value');
            fprintf('\n');
            % Middle rule.
            fprintf(repmat('–',1,NDash));
            fprintf('\n');
            for iProperty = 1:NProperty

                % Print property name.
                fprintf(['%-' num2str(max_name) 's'], KeySet{iProperty});
                fprintf(ColBreak);
                % Print property value string representation.
                fprintf(['%' num2str(max_val) 's'], StrValSet{iProperty});
                fprintf('\n');

            end

            % Bottom rule.
            fprintf(repmat('=',1,NDash));
            fprintf('\n');

        end

        function set_property(ob, Name, Value)
        % Add a custom property to the tree model. Existing property 
        % value will be updated.
        % 
        % Inputs:
        %
        % Name      String containing the name of the new property.
        % 
        % Value     Value of the property. Can have any type.

            % Name must be a string.
            assert(ischar(Name), 'Property name must be a string.');

            % Add new value to map.
            ob.custom_properties(Name) = Value;

        end

        function remove_property(ob, Name)
        % Remove a custom property by name. Use .LIST_PROPERTIES
        % to see existing property names. Non-existent property
        % name throws error.
        %
        % Inputs:
        %
        % Name      String containing the name of an existing property.

            % If key not found in map, error.
            assert( ...
                isKey(ob.custom_properties, Name), ...
                ['No property named ''' Name ''' found.'] ...
            );

            % Remove key and value from map.
            remove(ob.custom_properties, Name);

        end

        function Value = get_property(ob, Name)
        % Get the value of a custom property by name. Use .LIST_PROPERTIES
        % to see existing property names. Non-existent property
        % name throws error.

            % If key not found in map, error.
            assert( ...
                isKey(ob.custom_properties, Name), ...
                ['No property named ''' Name ''' found.'] ...
            );

            % Return property value.
            Value = ob.custom_properties(Name);

        end

        function FHasProperty = has_property(ob, Name)
        % Check if custom property exists by name. Return TRUE if
        % value exists and FALSE otherwise.

            % Check if key exists in map.
            FHasProperty = isKey(ob.custom_properties, Name);

        end

        function ob = compress(ob)
        % Compress size by deleting values that can be derived.

            % Update flag property.
            ob.compressed = true;
        end

        function ob = uncompress(ob)
        % Re-compute derived values.

            % Update flag property.
            ob.compressed = false;
        end
        
    end
    
    methods (Abstract)
        
        % Get selected properties of the blocks of the model.
        Props = get_block_properties(ob,PropNames,vargin)
        
        % Generate twigs connecting leaves to the blocks.
        [TwigStart, TwigEnd] = generate_twigs(ob, LeafParent)

        % Detect intersection between a block and a triangle.
        IntersectFlag = block_triangle_intersection(ob, JBlock, Tris)

        % Convert QSM into simple voxelization that shows which voxels are
        % occupied by any part of the blocks.
        BlockVoxelization = toVoxels(ob, edge, minp, maxp)

        % Return the origin point of the tree model.
        Origin = origin(ob)

    end

end