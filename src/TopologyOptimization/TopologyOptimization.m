classdef TopologyOptimization < handle
    % TopologyOptimization: a class that implements the topology
    % optimization algorithm.
    % Properties:
    %   nDim: size of the domain. 2 for 2D and 3 for 3D.
    %   nel: array with the number of elements in each direction.
    %   nElements: total number of elements.
    %   coord: matrix of size nElements x nDim with the coordinates of the
    %       elements centroids.
    %   freeElements: free elements indicator: 1 free, 0 fixed.
    %   mapFea2To: map from the FEA domain to the TO domain.
    %   mapTo2Fea: map from the TO domain to the FEA domain.
    %   d, d_filt, d_proj, d_simp, d_ramp: vector of size nElements x 1
    %       with the densities.
    %   radius: radius of the convolution filter.
    %   H: filter matrix.
    %   beta, eta: parameters for the projection.
    %   dMinSimp, p: parameters for SIMP.
    %   dMinRamp, q: parameters for RAMP.
    % Methods:
    %   TopologyOptimization: constructor.
    %   initialize_filter_2d: initialize the filter for 2D design domains.
    %   initialize_filter_3d: initialize the filter for 3D design domains.
    %   initialize_density: initialize the density field.
    %   set_initial_density_box: set the initial density field inside a
    %       box.
    %   set_density_box: set the density field inside a box and mark the
    %       elements as fixed.
    %   set_initial_density_sphere: set the initial density field inside a
    %       sphere.
    %   set_density_sphere: set the density field inside a sphere and mark
    %       the elements as fixed.
    %   filter: apply the filter.
    %   filter_sensitivity: compute the filter sensitivity.
    %   projection: apply the projection.
    %   projection_sensivitity: compute the projection sensitivity.
    %   simp: apply the SIMP.
    %   simp_sensitivity: compute the SIMP sensitivity.
    %   ramp: apply the RAMP.
    %   ramp_sensitivity: compute the RAMP sensitivity.

    properties
        nDim
        nel
        nElements
        coord
        freeElements
        mapFea2To, mapTo2Fea
        d, d_filt, d_proj, d_simp, d_ramp
        radius, H
        beta, eta
        dMinSimp, p
        dMinRamp, q
    end

    methods
        function obj = TopologyOptimization(nel, coord, radius, beta, eta, dMinSimp, p, varargin)
            % TopologyOptimization: constructor.
            % Inputs:
            %   nel: vector of size 1 x nDim, the number of elements in
            %       each direction.
            %   coord: matrix of size nElements x nDim with the coordinates
            %       of the elements centroids.
            %   radius: radius of the convolution filter.
            %   beta: parameter for the projection.
            %   eta: parameter for the projection.
            %   dMinSimp: parameter for the SIMP.
            %   p: parameter for the SIMP.
            %   dMinRamp: parameter for the RAMP (optional, default 1e-6).
            %   q: parameter for the RAMP (optional, default 4).
            % Outputs:
            %   obj: the TopologyOptimization object.

            % Parse inputs
            parser = inputParser;
            addOptional(parser, 'dMinRamp', 1e-6);
            addOptional(parser, 'q', 4);
            parse(parser, varargin{:});
            dMinRamp = parser.Results.dMinRamp;
            q = parser.Results.q;

            % Store properties
            obj.nDim = length(nel);
            obj.nel = nel(:).';
            obj.nElements = size(coord, 1);
            obj.coord = coord;
            obj.radius = radius;
            obj.beta = beta;
            obj.eta = eta;
            obj.dMinSimp = dMinSimp;
            obj.p = p;
            obj.dMinRamp = dMinRamp;
            obj.q = q;

            % Initialize free elements indicator (all free by default)
            obj.freeElements = ones(obj.nElements, 1);

            % Initialize densities
            obj.d      = zeros(obj.nElements, 1);
            obj.d_filt = zeros(obj.nElements, 1);
            obj.d_proj = zeros(obj.nElements, 1);
            obj.d_simp = zeros(obj.nElements, 1);
            obj.d_ramp = zeros(obj.nElements, 1);

            % Map between FEA and TO domains
            [obj.mapFea2To, obj.mapTo2Fea] = map_fea_2_to(obj.coord);

            % Initialize filter
            if obj.nDim == 2
                obj.initialize_filter_2d;
            elseif obj.nDim == 3
                obj.initialize_filter_3d;
            else
                error('The dimension of the domain should be 2 or 3.');
            end
        end

        function initialize_filter_2d(obj)
            % Initialize the filter for 2D design domains.

            % N-dimensional element map
            mapElements = reshape(obj.mapFea2To, obj.nel);

            % Initialize the triplets of the filter matrix
            nTriplets = prod(obj.nel) * (2*(ceil(obj.radius) - 1) + 1)^2;
            iRow = ones(nTriplets, 1);
            iCol = iRow;
            vals = zeros(nTriplets, 1);

            % Counter
            k = 1;

            % Loop over 2nd direction
            for y = 1:obj.nel(2)
                % Find limits along 2nd direction
                yMin = max(y - obj.radius, 1);
                yMax = min(y + obj.radius, obj.nel(2));
                yVec = yMin:yMax;

                % Loop over 1st direction
                for x = 1:obj.nel(1)
                    % Find limits along 1st direction
                    xMin = max(x - obj.radius, 1);
                    xMax = min(x + obj.radius, obj.nel(1));
                    xVec = xMin:xMax;

                    % Find indices of the neighbor elements
                    neighborIndices = mapElements(xVec, yVec);

                    % Compute relative coordinates
                    [dy, dx] = meshgrid(yVec - y, xVec - x);

                    % Compute relative distances
                    dist = sqrt(dx.^2 + dy.^2);

                    % Compute weights
                    w = obj.radius - dist(:);

                    % Current element index
                    idx = mapElements(x, y);

                    % Extract indices associated to a positive weight
                    positiveIndices = neighborIndices(w > 0);
                    positiveWeights = w(w > 0);
                    positiveWeights = positiveWeights / sum(positiveWeights); % normalize weights

                    % Number of positive weights
                    nPositive = length(positiveIndices);
    
                    % Store the triplets of the filter matrix
                    iRow(k:k+nPositive-1) = idx;
                    iCol(k:k+nPositive-1) = positiveIndices;
                    vals(k:k+nPositive-1) = positiveWeights;
                    k = k + nPositive;
                end
            end

            % Check the allocated memory
            if k-1 > nTriplets
                warning('The number of triplets (%d) for the filter matrix exceeds the allocated memory (%d).', k-1, nTriplets);
            end

            % Assemble the filter matrix
            obj.H = sparse(iRow, iCol, vals, obj.nElements, obj.nElements);
        end

        function initialize_filter_3d(obj)
            % Initialize the filter for 3D design domains.

            % N-dimensional element map
            mapElements = reshape(obj.mapFea2To, obj.nel);

            % Initialize the triplets of the filter matrix
            nTriplets = prod(obj.nel) * (2*(ceil(obj.radius) - 1) + 1)^3;
            iRow = ones(nTriplets, 1);
            iCol = iRow;
            vals = zeros(nTriplets, 1);

            % Counter
            k = 1;

            % Loop over 3rd direction
            for z = 1:obj.nel(3)
                % Find limits along 3rd direction
                zMin = max(z - obj.radius, 1);
                zMax = min(z + obj.radius, obj.nel(3));
                zVec = zMin:zMax;

                % Loop over 2nd direction
                for y = 1:obj.nel(2)
                    % Find limits along 2nd direction
                    yMin = max(y - obj.radius, 1);
                    yMax = min(y + obj.radius, obj.nel(2));
                    yVec = yMin:yMax;

                    % Loop over 1st direction
                    for x = 1:obj.nel(1)
                        % Find limits along 1st direction
                        xMin = max(x - obj.radius, 1);
                        xMax = min(x + obj.radius, obj.nel(1));
                        xVec = xMin:xMax;

                        % Find indices of the neighbor elements
                        neighborIndices = mapElements(xVec, yVec, zVec);

                        % Compute relative coordinates
                        [dz, dy, dx] = meshgrid(zVec - z, yVec - y, xVec - x);

                        % Compute relative distances
                        dist = sqrt(dx.^2 + dy.^2 + dz.^2);

                        % Compute weights
                        w = obj.radius - dist(:);

                        % Physical element index
                        idx = mapElements(x, y, z);

                        % Extract indices associated to a positive weight
                        positiveIndices = neighborIndices(w > 0);
                        positiveWeights = w(w > 0);
                        positiveWeights = positiveWeights / sum(positiveWeights); % normalize weights

                        % Number of positive weights
                        nPositive = length(positiveIndices);
    
                        % Store the triplets of the filter matrix
                        iRow(k:k+nPositive-1) = idx;
                        iCol(k:k+nPositive-1) = positiveIndices;
                        vals(k:k+nPositive-1) = positiveWeights;
                        k = k + nPositive;
                    end
                end
            end

            % Check the allocated memory
            if k-1 > nTriplets
                warning('The number of triplets (%d) for the filter matrix exceeds the allocated memory (%d).', k-1, nTriplets);
            end

            % Assemble the filter matrix
            obj.H = sparse(iRow, iCol, vals, obj.nElements, obj.nElements);
        end

        function initialize_density(obj, dInit)
            % Initialize the density field.
            % Inputs:
            %   dInit: vector of size nElements x 1 or scalar, the initial
            %       density.

            obj.d(:) = dInit;
        end

        function set_initial_density_box(obj, coord, tol, val)
            % Set the initial density field inside a box.
            % Inputs:
            %   coord: vector of size 1 x nDim, the coordinates of the
            %       center of the box.
            %   tol: vector of size 1 x nDim, the half-sizes of the box.
            %   val: scalar, the value of the density.

            dist = obj.coord - coord;
            isInside = all(abs(dist) < tol, 2);
            obj.d(isInside) = val;
        end

        function set_density_box(obj, coord, tol, val)
            % Set the density field inside a box and mark the elements as
            % fixed.
            % Inputs:
            %   coord: vector of size 1 x nDim, the coordinates of the
            %       center of the box.
            %   tol: vector of size 1 x nDim, the half-sizes of the box.
            %   val: scalar, the value of the density.

            dist = obj.coord - coord;
            isInside = all(abs(dist) < tol, 2);
            obj.d(isInside) = val;
            obj.freeElements(isInside) = 0;
        end

        function set_initial_density_sphere(obj, coord, tol, val)
            % Set the initial density field inside a sphere.
            % Inputs:
            %   coord: vector of size 1 x nDim, the coordinates of the
            %       center of the sphere.
            %   tol: scalar, the radius of the sphere.
            %   val: scalar, the value of the density.

            dist = obj.coord - coord;
            isInside = sum(dist.^2, 2) < tol^2;
            obj.d(isInside) = val;
        end

        function set_density_sphere(obj, coord, tol, val)
            % Set the density field inside a sphere and mark the elements
            % as fixed.
            % Inputs:
            %   coord: vector of size 1 x nDim, the coordinates of the
            %       center of the sphere.
            %   tol: scalar, the radius of the sphere.
            %   val: scalar, the value of the density.

            dist = obj.coord - coord;
            isInside = sum(dist.^2, 2) < tol^2;
            obj.d(isInside) = val;
            obj.freeElements(isInside) = 0;
        end

        function filter(obj)
            % Apply the filter.

            obj.d_filt = obj.H * obj.d;
        end

        function filtSens = filter_sensitivity(obj, sens)
            % Compute the filter sensitivity (pdv{J}{d}).
            % Inputs:
            %   sens: vector of size nElements x 1, the sensitivity.
            % Outputs:
            %   filtSens: vector of size nElements x 1, the filter
            %       sensitivity.

            % Apply projection sensitivity
            sens = sens .* obj.projection_sensivitity();

            % Apply filter sensitivity
            filtSens = obj.H.' * sens;
        end

        function projection(obj)
            % Apply the projection.

            obj.d_proj = (tanh(obj.beta*obj.eta) + tanh(obj.beta*(obj.d_filt - obj.eta))) / ...
                (tanh(obj.beta*obj.eta) + tanh(obj.beta*(1 - obj.eta)));
        end

        function projSens = projection_sensivitity(obj)
            % Compute the projection sensitivity (pdv{d_proj}{d_filt}).
            % Outputs:
            %   projSens: vector of size nElements x 1, the projection
            %       sensitivity.

            projSens = obj.beta * (1 - tanh(obj.beta*(obj.d_filt - obj.eta)).^2) / ...
                (tanh(obj.beta*obj.eta) + tanh(obj.beta*(1 - obj.eta)));
        end

        function simp(obj)
            % Apply the SIMP.

            obj.d_simp = obj.dMinSimp + (1 - obj.dMinSimp) * obj.d_proj .^ obj.p;
        end

        function simpSens = simp_sensitivity(obj)
            % Compute the SIMP sensitivity (pdv{d_simp}{d_proj}).
            % Outputs:
            %   simpSens: vector of size nElements x 1, the SIMP
            %       sensitivity.

            simpSens = obj.p * (1 - obj.dMinSimp) * obj.d_proj .^ (obj.p - 1);
        end

        function ramp(obj)
            % Apply the RAMP.

            obj.d_ramp = obj.dMinRamp + (1 - obj.dMinRamp) * obj.d_proj ./ (1 + obj.q * (1 - obj.d_proj));
        end

        function rampSens = ramp_sensitivity(obj)
            % Compute the RAMP sensitivity (pdv{d_ramp}{d_proj}).
            % Outputs:
            %   rampSens: vector of size nElements x 1, the RAMP
            %       sensitivity.

            rampSens = (1 - obj.dMinRamp) * (1 + obj.q) ./ (1 + obj.q * (1 - obj.d_proj)).^2;
        end
    end
end
