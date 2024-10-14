classdef SymmetryMap
    % SymmetryMap: A class to handle domain symmetries.
    % Properties:
    %   - mapSymmetry: the dictionary of symmetry groups.
    %   - nElements: the total number of elements.
    % Methods:
    %   - symmetry_orthogonal: constructor for orthogonal symmetries.
    %   - symmetry_extraction: from the full space to the reduced symmetric space (product).
    %   - symmetry_reduction: from the full space to the reduced symmetric space (sum).
    %   - symmetry_extension: from the reduced symmetric space to the full space.
    % Static methods:
    %   - symmetry_along_X: compute the symmetry groups for symmetry along the X axis.
    %   - symmetry_along_Y: compute the symmetry groups for symmetry along the Y axis.
    %   - symmetry_along_XY: compute the symmetry groups for symmetry along the X and Y axes.

    properties
        mapSymmetry
        nElements
    end

    methods
        function obj = symmetry_orthogonal(obj, coord, nel, symmetry_type)
            % Compute the symmetry groups for orthogonal symmetry.
            % Inputs:
            %   coord: matrix of size n x 2, where n is the number of elements.
            %          The columns are the x and y coordinates of the element centroids.
            %   nel: vector of size 1 x 2, the number of elements in each direction.
            %   symmetry_type: string, the type of symmetry ('x', 'y', 'xy'). Default is 'xy'.

            % Check symmetry type
            if nargin < 4
                symmetry_type = 'xy';
            end
            switch symmetry_type
                case 'x'
                    symmetry_function = @obj.symmetry_along_X;
                case 'y'
                    symmetry_function = @obj.symmetry_along_Y;
                case 'xy'
                    symmetry_function = @obj.symmetry_along_XY;
                otherwise
                    error("Invalid symmetry type.")
            end

            % Store number of elements
            obj.nElements = prod(nel);

            % Element size
            dEl = (max(coord) + min(coord)) ./ nel;

            % Generate adimensional coordinates
            coordAdim = round(coord ./ dEl  - 0.5);

            % Initialize dictionary
            obj.mapSymmetry = configureDictionary('double', 'cell');

            % Loop over elements
            for ii = 1:size(coord, 1)
                thisCoordAdim = coordAdim(ii, :);

                % Compute the minimum distance to the edge of the domain
                newCoordAdim = symmetry_function(thisCoordAdim, nel);

                % Compute the code according to the lexographical order
                code = 0;
                for jj = length(newCoordAdim):-1:1
                    code = nel(jj) * code + newCoordAdim(jj);
                end

                % Add entry
                if obj.mapSymmetry.isKey(code) == 0
                    obj.mapSymmetry = obj.mapSymmetry.insert(code, {ii});
                else
                    vals = obj.mapSymmetry(code);
                    obj.mapSymmetry(code) = {[vals{:}, ii]};
                end
            end

            % Print info
            fprintf("\nSymmetry groups: %d\n", obj.mapSymmetry.numEntries)
        end

        function obj = symmetry_diagonal(obj, coord, nel)
            % Compute the symmetry groups for diagonal symmetry.
            % Inputs:
            %   coord: matrix of size n x 2, where n is the number of elements.
            %          The columns are the x and y coordinates of the element centroids.
            %   nel: vector of size 1 x 2, the number of elements in each direction.

            % Check if the design domain is a square
            if ~all(nel(1) == nel)
                error("The design domain must be a square.")
            end

            % Store number of elements
            obj.nElements = prod(nel);

            % Element size
            dEl = (max(coord) + min(coord)) ./ nel;

            % Generate adimensional coordinates
            coordAdim = round(coord ./ dEl  - 0.5);

            % Initialize dictionary
            obj.mapSymmetry = configureDictionary('double', 'cell');

            % Loop over elements
            for ii = 1:size(coord, 1)
                thisCoordAdim = coordAdim(ii, :);

                % Switch the coordinates if y is greater than x
                if thisCoordAdim(2) > thisCoordAdim(1)
                    thisCoordAdim = fliplr(thisCoordAdim);
                end

                % Compute the code according to the lexographical order
                code = 0;
                for jj = length(thisCoordAdim):-1:1
                    code = nel(jj) * code + thisCoordAdim(jj);
                end

                % Add entry
                if obj.mapSymmetry.isKey(code) == 0
                    obj.mapSymmetry = obj.mapSymmetry.insert(code, {ii});
                else
                    vals = obj.mapSymmetry(code);
                    obj.mapSymmetry(code) = {[vals{:}, ii]};
                end
            end

            % Print info
            fprintf("\nSymmetry groups: %d\n", obj.mapSymmetry.numEntries)
        end

        function obj = symmetry_octagonal(obj, coord, nel)
            % Compute the symmetry groups for octagonal symmetry. This is the combination of the diagonal and orthogonal xy symmetries.
            % Inputs:
            %   coord: matrix of size n x 2, where n is the number of elements.
            %          The columns are the x and y coordinates of the element centroids.
            %   nel: vector of size 1 x 2, the number of elements in each direction.

            % Check if the design domain is a square
            if ~all(nel(1) == nel)
                error("The design domain must be a square.")
            end

            % Store number of elements
            obj.nElements = prod(nel);

            % Element size
            dEl = (max(coord) + min(coord)) ./ nel;

            % Generate adimensional coordinates
            coordAdim = round(coord ./ dEl  - 0.5);

            % Initialize dictionary
            obj.mapSymmetry = configureDictionary('double', 'cell');

            % Loop over elements
            for ii = 1:size(coord, 1)
                thisCoordAdim = coordAdim(ii, :);

                % Compute the minimum distance to the edge of the domain
                thisCoordAdim = obj.symmetry_along_XY(thisCoordAdim, nel);

                % Switch the coordinates if y is greater than x
                if thisCoordAdim(2) > thisCoordAdim(1)
                    thisCoordAdim = fliplr(thisCoordAdim);
                end

                % Compute the code according to the lexographical order
                code = 0;
                for jj = length(thisCoordAdim):-1:1
                    code = nel(jj) * code + thisCoordAdim(jj);
                end

                % Add entry
                if obj.mapSymmetry.isKey(code) == 0
                    obj.mapSymmetry = obj.mapSymmetry.insert(code, {ii});
                else
                    vals = obj.mapSymmetry(code);
                    obj.mapSymmetry(code) = {[vals{:}, ii]};
                end
            end

            % Print info
            fprintf("\nSymmetry groups: %d\n", obj.mapSymmetry.numEntries)
        end

        function vRed = symmetry_extraction(obj, v)
            % From the full space to the reduced symmetric space (product).
            % Inputs:
            %   v: matrix of size n x m, where n is the number of elements and m is the number of columns.
            % Outputs:
            %   vRed: matrix of size numEntries x m, the reduced symmetric space.

            n_cols = size(v, 2);
            values = obj.mapSymmetry.values;
            vRed = zeros(obj.mapSymmetry.numEntries, n_cols);
            for ii = 1:obj.mapSymmetry.numEntries
                vRed(ii, :) = v(values{ii}(1), :);
            end
        end

        function vRed = symmetry_reduction(obj, v)
            % From the full space to the reduced symmetric space (sum).
            % Inputs:
            %   v: matrix of size n x m, where n is the number of elements and m is the number of columns.
            % Outputs:
            %   vRed: matrix of size numEntries x m, the reduced symmetric space.

            n_cols = size(v, 2);
            values = obj.mapSymmetry.values;
            vRed = zeros(obj.mapSymmetry.numEntries, n_cols);
            for ii = 1:obj.mapSymmetry.numEntries
                vRed(ii, :) = sum(v(values{ii}, :));
            end
        end

        function v = symmetry_extension(obj, vRed)
            % From the reduced symmetric space to the full space.
            % Inputs:
            %   vRed: matrix of size numEntries x m, the reduced symmetric space.
            % Outputs:
            %   v: matrix of size n x m, where n is the number of elements and m is the number of columns.
            
            n_cols = size(vRed, 2);
            values = obj.mapSymmetry.values;
            v = zeros(obj.nElements, n_cols);
            for ii = 1:obj.mapSymmetry.numEntries
                v(values{ii}, :) = vRed(ii, :);
            end
        end
    end

    methods(Static, Access = private)
        function newCoordAdim = symmetry_along_X(thisCoordAdim, nel)
            % Compute the symmetry groups for symmetry along the X axis.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in each direction.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = thisCoordAdim;
            newCoordAdim(2) = min(newCoordAdim(2), nel(2) - 1 - thisCoordAdim(2));
        end

        function newCoordAdim = symmetry_along_Y(thisCoordAdim, nel)
            % Compute the symmetry groups for symmetry along the Y axis.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in each direction.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = thisCoordAdim;
            newCoordAdim(1) = min(newCoordAdim(1), nel(1) - 1 - thisCoordAdim(1));
        end

        function newCoordAdim = symmetry_along_XY(thisCoordAdim, nel)
            % Compute the symmetry groups for symmetry along the X and Y axes.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in each direction.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = min(thisCoordAdim, nel - 1 - thisCoordAdim);
        end
    end
end

