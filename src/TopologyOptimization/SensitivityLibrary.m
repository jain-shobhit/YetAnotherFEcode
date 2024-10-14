classdef SensitivityLibrary
    % SensitivityLibrary: A collection of sensitivity analysis methods
    % for topology optimization problems.
    % Methods:
    %   - compliance: compliance sensitivity
    %   - frequency: frequency sensitivity

    methods (Static)
        function sens = compliance(myMesh, u, Ke, sensPh)
            % Compute the compliance sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh (optional): vector of size m x 1, the sensitivity of the physical density with respect to the projected density.
            % Outputs:
            %   sens: vector of size m x 1, the compliance sensitivity.

            % Check if the sensitivity of the physical density is provided
            if nargin < 4
                sensPh = ones(myMesh.nElements, 1);
            end

            sens = zeros(myMesh.nElements, 1);
            for ii = 1:myMesh.nElements
                iDOFs = myMesh.Elements(ii).Object.iDOFs;
                ue = u(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = -ue.' * dK * ue;
            end
        end

        function sens = frequency(myMesh, omega, u, Ke, Me, sensPhK, sensPhM)
            % Compute the frequency sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   omega: scalar, the eigenfrequency.
            %   u: vector of size n x 1, the mode shape.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   Me: matrix of size nE x nE, the element mass matrix.
            %   sensPhK (optional): vector of size m x 1, the sensitivity of the physical density with respect to the projected density for the stiffness matrix.
            %   sensPhM (optional): vector of size m x 1, the sensitivity of the physical density with respect to the projected density for the mass matrix.
            % Outputs:
            %   sens: vector of size m x 1, the frequency sensitivity.
            % Notes:
            %   - Make sure the mode shape u is mass-normalized.

            % Check if the sensitivity of the physical density is provided
            if nargin < 6
                sensPhK = ones(myMesh.nElements, 1);
                sensPhM = ones(myMesh.nElements, 1);
            end
            
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
