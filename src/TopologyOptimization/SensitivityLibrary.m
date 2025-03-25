classdef SensitivityLibrary
    % SensitivityLibrary: A collection of sensitivity analysis methods
    % for topology optimization problems.
    % Methods:
    %   compliance: compliance sensitivity.
    %   compliance_out: output compliance sensitivity.
    %   frequency: frequency sensitivity.

    methods (Static)
        function sens = compliance(myMesh, u, Ke, varargin)
            % Compute the compliance sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected density
            %       (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the compliance sensitivity.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPh', ones(myMesh.nElements, 1));
            parse(p, varargin{:});
            sensPh = p.Results.sensPh;

            % Loop over the elements
            sens = zeros(myMesh.nElements, 1);
            for ii = 1:myMesh.nElements
                iDOFs = myMesh.Elements(ii).Object.iDOFs;
                ue = u(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = -ue.' * dK * ue;
            end
        end

        function sens = compliance_out(myAssembly, Kc, fOutC, u, Ke, varargin)
            % Compute the output compliance sensitivity.
            % Inputs:
            %   myAssembly: the Assembly object.
            %   uc: vector of size nC x 1, the displacement field.
            %   Kc: matrix of size nC x nC, the stiffness matrix (with the
            %       lumped spring).
            %   fOutC: vector of size nC x 1, the output force vector.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected density
            %       (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the compliance sensitivity.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPh', ones(myAssembly.Mesh.nElements, 1));
            parse(p, varargin{:});
            sensPh = p.Results.sensPh;

            % Adjoint variable
            adjC = -Kc \ fOutC;
            adj = myAssembly.unconstrain_vector(adjC);

            % Loop over the elements
            sens = zeros(myAssembly.Mesh.nElements, 1);
            for ii = 1:myAssembly.Mesh.nElements
                iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;
                uE = u(iDOFs, 1);
                adjE = adj(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = adjE.' * dK * uE;
            end
        end

        function sens = frequency(myMesh, omega, u, Ke, Me, varargin)
            % Compute the frequency sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   omega: scalar, the eigenfrequency.
            %   u: vector of size n x 1, the mode shape.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   Me: matrix of size nE x nE, the element mass matrix.
            %   sensPhK: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the stiffness matrix (optional, default is ones(m, 1)).
            %   sensPhM: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the mass matrix (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the frequency sensitivity.
            % Notes:
            %   Make sure the mode shape u is mass-normalized.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPhK', ones(myMesh.nElements, 1));
            addOptional(p, 'sensPhM', ones(myMesh.nElements, 1));
            parse(p, varargin{:});
            sensPhK = p.Results.sensPhK;
            sensPhM = p.Results.sensPhM;

            % Loop over the elements
            sens = zeros(myMesh.nElements, 1);
            for ii = 1:myMesh.nElements
                iDOFs = myMesh.Elements(ii).Object.iDOFs;
                ue = u(iDOFs, 1);
                dK = sensPhK(ii) * Ke;
                dM = sensPhM(ii) * Me;
                sens(ii, :) = ue.' * (dK - omega^2 * dM) * ue;
            end
            sens = sens / (2 * omega) / (2 * pi);
        end
    end
end
