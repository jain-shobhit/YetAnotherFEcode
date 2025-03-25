classdef Hex8Element_EM < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 8          % number of nodes per element
        nDim = 3            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'HEX'
    end
    
    properties
        thickness = 1       % element thickness, by default zero    
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Hex8Element_EM(Material, Ngauss)
            % Self function (constructor)
            if nargin == 1
                Ngauss = 2;
            end
            self.thickness = 1;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;            
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self, X)
            % dH: shape function derivatives in physical coordinates
            % detJ = det(J), J=jacobian
            % G=0: unused, left for compatibility with other methods
            %______________________________________________________________
            g = X(1);
            h = X(2);
            r = X(3);
            xyz = self.nodes;
            % shape functions derivatives (ABAQUS ORDER)
            dHn = [ -((h - 1)*(r - 1))/8    -((g - 1)*(r - 1))/8    -((g - 1)*(h - 1))/8
                     ((h - 1)*(r - 1))/8     ((g + 1)*(r - 1))/8     ((g + 1)*(h - 1))/8
                    -((h + 1)*(r - 1))/8    -((g + 1)*(r - 1))/8    -((g + 1)*(h + 1))/8
                     ((h + 1)*(r - 1))/8     ((g - 1)*(r - 1))/8     ((g - 1)*(h + 1))/8
                     ((h - 1)*(r + 1))/8     ((g - 1)*(r + 1))/8     ((g - 1)*(h - 1))/8
                    -((h - 1)*(r + 1))/8    -((g + 1)*(r + 1))/8    -((g + 1)*(h - 1))/8
                     ((h + 1)*(r + 1))/8     ((g + 1)*(r + 1))/8     ((g + 1)*(h + 1))/8
                    -((h + 1)*(r + 1))/8    -((g - 1)*(r + 1))/8    -((g - 1)*(h + 1))/8].';
            J = dHn*xyz;
            J1 = [0 0 0; 0 0 0; 0 0 0];
            J1(1,1) = J(2,2)*J(3,3) - J(2,3)*J(3,2);
            J1(2,1) = J(2,3)*J(3,1) - J(2,1)*J(3,3);
            J1(3,1) = J(2,1)*J(3,2) - J(2,2)*J(3,1);
            J1(1,2) = J(1,3)*J(3,2) - J(1,2)*J(3,3);
            J1(2,2) = J(1,1)*J(3,3) - J(1,3)*J(3,1);
            J1(3,2) = J(1,2)*J(3,1) - J(1,1)*J(3,2);
            J1(1,3) = J(1,2)*J(2,3) - J(1,3)*J(2,2);
            J1(2,3) = J(1,3)*J(2,1) - J(1,1)*J(2,3);
            J1(3,3) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
            detJ = J(1,1)*J1(1,1) + J(1,2)*J1(2,1) + J(1,3)*J1(3,1);
            J1 = J1/detJ;
            dH = J1*dHn;    % derivatives in physical coordinates,
                            % 3x8 matrix, [dNi_dx; dNi_dy; dNi_dz]
                            % with i = 1...8
            G=0;    % not used, but we keep it for compatibility with other 
                    % methods (e.g. get.area)
        end
        
        % ELECTROSTATICS __________________________________________________

        function Kel = electrostatic_stiffness_matrix(self)
            X = self.quadrature.X;
            W = self.quadrature.W;
            Kel = self.initialization.K;
            for ii = 1:length(W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [~,detJ,dH] = shape_function_derivatives(self, Xi);        
                intKel = dH'*dH*detJ;
                % integration of K and F through Gauss quadrature
                Kel = Kel + we*intKel;
            end
            Kel = Kel * self.Material.permittivity;
        end
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED HEXAHEDRON
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            r = X(3);            
            N = 1/8*[(1-g)*(1-h)*(1-r); (1+g)*(1-h)*(1-r);
                     (1+g)*(1+h)*(1-r); (1-g)*(1+h)*(1-r);
                     (1-g)*(1-h)*(1+r); (1+g)*(1-h)*(1+r);
                     (1+g)*(1+h)*(1+r); (1-g)*(1+h)*(1+r)];
        end
        
        function X = natural_coordinates
            X = [-1 -1 -1
                 +1 -1 -1
                 +1 +1 -1
                 -1 +1 -1
                 -1 -1 +1
                 +1 -1 +1
                 +1 +1 +1
                 -1 +1 +1];
        end
        
        function edges = edge_nodes
            edges = [1 2; 2 3; 3 4; 4 1; ...
                5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
        end

        function [faces, normal_vectors] = face_nodes
            % face nodes are ordered ccw according to the right hand rule,
            % with the thumb oriented as the outward normal to the face
            faces = {[1 2 6 5;
                2 3 7 6;
                3 4 8 7;
                1 5 8 4;
                5 6 7 8;
                1 4 3 2]};
            normal_vectors = [0 -1 0;
                1 0 0;
                0 1 0;
                -1 0 0;
                0 0 1;
                0 0 -1];
        end
    end
        
end % classdef
    
