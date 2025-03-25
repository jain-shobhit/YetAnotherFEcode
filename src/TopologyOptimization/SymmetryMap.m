classdef SymmetryMap
    % SymmetryMap: A class to handle domain symmetries.
    % Properties:
    %   mapSymmetry: the dictionary of symmetry groups.
    %   nElements: the total number of elements.
    % Methods:
    %   SymmetryMap: the constructor.
    %   symmetry_extraction: from the full space to the reduced symmetric
    %       space.
    %   symmetry_reduction: from the full space to the reduced symmetric
    %       space (sum).
    %   symmetry_extension: from the reduced symmetric space to the full
    %       space.
    % Static methods:
    %   symmetry_along_X: symmetry along the X axis.
    %   symmetry_along_Y: symmetry along the Y axis.
    %   symmetry_along_XY: symmetry along the X and Y axes.
    %   symmetry_diagonal: symmetry along the diagonal.
    %   symmetry_octagonal: symmetry along the X, Y and diagonal axes.

    properties
        mapSymmetry
        nElements
    end

    methods
        function obj = SymmetryMap(coord, nel, varargin)
            % Compute the symmetry groups. The type of symmetries are:
            %   'x': symmetry along the X axis.
            %   'y': symmetry along the Y axis.
            %   'xy': symmetry along the X and Y axes.
            %   'd': symmetry along the diagonal.
            %   'o': symmetry along the X, Y and diagonal axes.
            % Inputs:
            %   coord: matrix of size n x 2, where n is the number of
            %       elements. The columns are the x and y coordinates of
            %       the element centroids.
            %   nel: vector of size 1 x 2, the number of elements in each
            %       direction.
            %   symmetry_type: string, the type of symmetry (optional,
            %       default is 'xy'). Options are 'x', 'y', 'xy', 'd', 'o'.
            %   timeIt: boolean, whether to time the function (optional,
            %       default is false).
            % Outputs:
            %   obj: the SymmetryMap object.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'symmetry_type', 'xy', ...
                @(x) ischar(x) && ismember(x, {'x', 'y', 'xy', 'd', 'o'}));
            addOptional(p, 'timeIt', false);
            parse(p, varargin{:});
            symmetry_type = p.Results.symmetry_type;
            
            % Start timer
            if p.Results.timeIt
                tStart = tic;
            end

            % Check the symmetry type
            if strcmp(symmetry_type, 'x')
                symmetry_function = @obj.symmetry_along_X;
            elseif strcmp(symmetry_type, 'y')
                symmetry_function = @obj.symmetry_along_Y;
            elseif strcmp(symmetry_type, 'xy')
                symmetry_function = @obj.symmetry_along_XY;
            elseif strcmp(symmetry_type, 'd')
                symmetry_function = @obj.symmetry_diagonal;
            elseif strcmp(symmetry_type, 'o')
                symmetry_function = @obj.symmetry_octagonal;
            else
                error("Invalid symmetry type.")
            end

            % Store number of elements
            obj.nElements = prod(nel);

            % Element size
            dEl = (max(coord) + min(coord)) ./ nel;

            % Generate adimensional coordinates
            coordAdim = round(coord ./ dEl  - 0.5);

            % Initialize symmetry map
            obj.mapSymmetry.keys = [];
            obj.mapSymmetry.values = {};

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
                if ~any(obj.mapSymmetry.keys == code)
                    obj.mapSymmetry.keys = [obj.mapSymmetry.keys; code];
                    obj.mapSymmetry.values{end + 1} = ii;
                else
                    idx = find(obj.mapSymmetry.keys == code);
                    obj.mapSymmetry.values{idx} = [obj.mapSymmetry.values{idx}, ii];
                end
            end

            % Print info
            obj.mapSymmetry.numEntries = length(obj.mapSymmetry.keys);
            fprintf("\nSymmetry groups: %d\n", obj.mapSymmetry.numEntries)

            % Stop timer
            if p.Results.timeIt
                toc(tStart);
            end
        end

        function vRed = symmetry_extraction(obj, v)
            % From the full space to the reduced symmetric space.
            % Inputs:
            %   v: matrix of size n x m, where n is the number of elements
            %       and m is the number of columns.
            % Outputs:
            %   vRed: matrix of size numEntries x m, the reduced symmetric
            %       space.

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
            %   v: matrix of size n x m, where n is the number of elements
            %       and m is the number of columns.
            % Outputs:
            %   vRed: matrix of size numEntries x m, the reduced symmetric
            %       space.

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
            %   vRed: matrix of size numEntries x m, the reduced symmetric
            %       space.
            % Outputs:
            %   v: matrix of size n x m, where n is the number of elements
            %       and m is the number of columns.
            
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
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional
            %       coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            % Outputs:
            %   newCoordAdim: vector of size 1 x nDim, the new adimensional
            %       coordinates of the element.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = thisCoordAdim;
            newCoordAdim(2) = min(newCoordAdim(2), nel(2) - 1 - thisCoordAdim(2));
        end

        function newCoordAdim = symmetry_along_Y(thisCoordAdim, nel)
            % Compute the symmetry groups for symmetry along the Y axis.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional
            %       coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            % Outputs:
            %   newCoordAdim: vector of size 1 x nDim, the new adimensional
            %       coordinates of the element.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = thisCoordAdim;
            newCoordAdim(1) = min(newCoordAdim(1), nel(1) - 1 - thisCoordAdim(1));
        end

        function newCoordAdim = symmetry_along_XY(thisCoordAdim, nel)
            % Compute the symmetry groups for symmetry along the X and Y
            % axes.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional
            %       coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            % Outputs:
            %   newCoordAdim: vector of size 1 x nDim, the new adimensional
            %       coordinates of the element.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = min(thisCoordAdim, nel - 1 - thisCoordAdim);
        end

        function newCoordAdim = symmetry_diagonal(thisCoordAdim, nel)
            % Compute the symmetry groups for diagonal symmetry.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional
            %       coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            % Outputs:
            %   newCoordAdim: vector of size 1 x nDim, the new adimensional
            %       coordinates of the element.

            % Switch the coordinates if y is greater than x
            if thisCoordAdim(2) > thisCoordAdim(1)
                newCoordAdim = fliplr(thisCoordAdim);
            else
                newCoordAdim = thisCoordAdim;
            end
        end

        function newCoordAdim = symmetry_octagonal(thisCoordAdim, nel)
            % Compute the symmetry groups for octagonal symmetry.
            % Inputs:
            %   thisCoordAdim: vector of size 1 x nDim, the adimensional
            %       coordinates of the element.
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            % Outputs:
            %   newCoordAdim: vector of size 1 x nDim, the new adimensional
            %       coordinates of the element.

            % Compute the minimum distance to the edge of the domain
            newCoordAdim = min(thisCoordAdim, nel - 1 - thisCoordAdim);

            % Switch the coordinates if y is greater than x
            if newCoordAdim(2) > newCoordAdim(1)
                newCoordAdim = fliplr(newCoordAdim);
            end
        end
    end
end
