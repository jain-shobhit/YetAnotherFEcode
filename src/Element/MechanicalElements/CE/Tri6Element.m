classdef Tri6Element < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 2     % number of DOFs per node
        nNodes = 6          % number of nodes per element
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
        
        function self = Tri6Element(thickness, Material, Ngauss)
            % Self function (constructor)
            if nargin == 2
                Ngauss = 2;
            end
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
            % and p={u1,v1,...,u6,v6}' (nodal displacements).
            % detJ = det(J), J=jacobian
            %______________________________________________________________
            g = X(1);
            h = X(2);
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = [ 4*g + 4*h - 3, 4*g + 4*h - 3
                          4*g - 1,             0
                                0,       4*h - 1
                    4 - 4*h - 8*g,          -4*g
                              4*h,           4*g
                             -4*h, 4 - 8*h - 4*g]';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x6 matrix, [dN_dx; dN_dy]
            G = self.initialization.G;
            G(1:2,1:2:end) = dH;
            G(3:4,2:2:end) = dH;
        end

    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h)
            % SHAPE FUNCTIONS FOR A 6-NODED TRIANGULAR
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.6 Triangular, tetrahedral, and wedge elements
            g = X(1);
            h = X(2);
            N = [...
                    2*(1/2-g-h)*(1-g-h);
                    2*g*(g-1/2);
                    2*h*(h-1/2);
                    4*g*(1-g-h);
                    4*g*h;
                    4*h*(1-g-h)];
        end
        
        function X = natural_coordinates
            X = [ ...
                0 0     % node 1 (corner)
                1 0     % node 2 (corner)
                0 1     % node 3 (corner)
                .5 0    % node 4
                .5 .5   % node 5
                0 .5];  % node 6
        end
        
        function [edges,normal_vectors] = edge_nodes
            edges = [1 4; 4 2; 2 5; 5 3; 3 6; 6 1]; % counterclockwise
            normal_vectors = [0 -1; 0 -1; sqrt(2)/2 sqrt(2)/2; sqrt(2)/2 sqrt(2)/2; -1 0; -1 0];
        end

        function faces = face_nodes
            % face nodes are ordered ccw according to the right hand rule,
            % with the thumb oriented as the outward normal to the face
            faces = {[1 4 2 5 3 6]};
        end

    end
        
end % classdef
    
