classdef Tri3Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 3          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'TRI'
    end
    
    properties
        thickness = 0       % element thickness, by default zero    
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Tri3Element(thickness, Material, Ngauss)
            % Self function (constructor)
            if nargin == 3
                warning('note: shape function derivatives are constants, and M is a lumped-parameter mass matrix. No quadrature integration is actually needed. Ngauss is thus set to 1.')
            end
            Ngauss = 1; % note: shape function derivatives are constants,
                        % and M is a lumped-parameter mass matrix. No
                        % quadrature integration is actually needed.
            self.thickness = thickness;
            self.nelDOFs = self.nNodes * self.nDOFPerNode;
            ContinuumElementConstructor(self, Material, Ngauss);
        end
        
        function [G,detJ,dH] = shape_function_derivatives(self, X)
            %______________________________________________________________
            %
            % [G,detJ,dH] = shape_function_derivatives(self, X)
            % G = shape function derivatives in physical coordinates, s.t.
            % th=G*p with th={ux uy vx vy}' (ux=du/dx...)
            % and p={u1,v1,...,u3,v3}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ -1 -1; 1 0; 0 1]';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x6 matrix, [dN_dx; dN_dy]
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end
        
        function Mel = mass_matrix(self)
            % _____________________________________________________________
            %
            % Mel = mass_matrix_global(self,nodes,~);
            % Mel: element-level mass matrix (in global coordinates)
            %______________________________________________________________
            rho = self.Material.DENSITY;
            t = self.thickness;
            m = self.area * rho * t / self.nNodes; % lumped masses are better
            Mel = sparse(eye(self.nNodes*self.nDOFPerNode)*m);
        end
    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h)
            % SHAPE FUNCTIONS FOR A 3-NODED TRIANGULAR
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            N = [1-g-h; g; h];
        end
        
        function X = natural_coordinates
            X = [ ...
                0 0     % node 1
                1 0     % node 2
                0 1];   % node 3
        end
        
        function [edges,normal_vectors] = edge_nodes
            edges = [1 2; 2 3; 3 1]; % counterclockwise
            normal_vectors = [0 -1; sqrt(2)/2 sqrt(2)/2; -1 0];
        end

        function faces = face_nodes
            % face nodes are ordered ccw according to the right hand rule,
            % with the thumb oriented as the outward normal to the face
            faces = {[1 2 3]};
        end

    end
        
end % classdef
    
