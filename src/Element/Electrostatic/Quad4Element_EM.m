classdef Quad4Element_EM < ContinuumElement
    properties
        nodes = []          % global coordinates of element nodes
        nodeIDs = []        % the index location of element nodes
        nDOFPerNode = 1     % number of DOFs per node
        nNodes = 4          % number of nodes per element
        nDim = 2            % number of dimensions in local coordinates
        nelDOFs
        quadrature       	% weights and points for gauss quadrature
        Material           	% Object of class Material
        initialization     	% some 0-matrices to speedup numerical integration
        elType = 'QUAD'
    end
    
    properties
        thickness = 0       % element thickness, by default zero    
    end
    
    methods
        
        % MINIMUM REQUIRED FUNCTIONS ______________________________________
        
        function self = Quad4Element_EM(thickness, Material, Ngauss)
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
            % dH: shape function derivatives in physical coordinates
            % detJ = det(J), J=jacobian
            % G=0: unused, left for compatibility with other methods
            %______________________________________________________________
            r = X(1);
            s = X(2);
            xy = self.nodes;
            % shape function derivatives in natural coordinates
            dHn = 1/4*[ s-1,  r-1; 
                        1-s, -r-1; 
                        s+1,  r+1; 
                       -s-1,  1-r].';
            J = dHn*xy;
            detJ = det(J);
            dH = J\dHn;	% derivatives in physical coordinates,
                      	% 2x4 matrix, [dNi_dx; dNi_dy]
                       	% with i = 1...4
            G=0; % not used, but we keep it for compatibility with other 
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

        function E = electric_field_x(self,x)
            % E = -grad(potential)
            potential = self.extract_element_data(x);
            X = self.natural_coordinates;
            E = zeros(4,1);
            for ii = 1:size(X,1)
                % evaluate on nodes
                [~,~,dH] = shape_function_derivatives(self, X(ii,:));     
                E(ii) = -dH(1,:)*potential;
            end
        end

        function E = electric_field_y(self,x)
            % E = -grad(potential)
            potential = self.extract_element_data(x);
            X = self.natural_coordinates;
            E = zeros(4,1);
            for ii = 1:size(X,1)
                % evaluate on nodes
                [~,~,dH] = shape_function_derivatives(self, X(ii,:));        
                E(ii) = -dH(2,:)*potential;
            end
        end

        function fel_x = electrostatic_force_x(self,x)
            phi = self.extract_element_data(x); % electric potential
            X = self.quadrature.X;
            W = self.quadrature.W;
            fel_x = zeros(4,1);
            for ii = 1:length(W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [~,detJ,dH] = shape_function_derivatives(self, Xi);
                dN = dH';               % dim(dN) = n x 2
                dNx = dN(:,1);          % dim(dNx)= n x 1
                dNy = dN(:,2);          % dim(dNy)= n x 1
                Bt = blkdiag(dN,dN);    % dim(B) = 2n x 4
                M = [phi'*dNx,  -phi'*dNy;
                     0,          2*phi'*dNx;
                     2*phi'*dNy, 0;
                     -phi'*dNx,  phi'*dNy];                 
                intFel = Bt*M*dN'*phi;
                % integration of K and F through Gauss quadrature
                fel_x = fel_x + we*intFel(1:4,1)*detJ;
            end
            fel_x = 1/2 * self.Material.permittivity * ...
                self.thickness * fel_x;
        end

        function fel_y = electrostatic_force_y(self,x)
            phi = self.extract_element_data(x); % electric potential
            X = self.quadrature.X;
            W = self.quadrature.W;
            fel_y = zeros(4,1);
            for ii = 1:length(W)
                Xi = X(:,ii);   % quadrature points
                we = W(ii);     % quadrature weights
                [~,detJ,dH] = shape_function_derivatives(self, Xi);
                dN = dH';               % dim(dN) = n x 2
                dNx = dN(:,1);          % dim(dNx)= n x 1
                dNy = dN(:,2);          % dim(dNy)= n x 1
                Bt = blkdiag(dN,dN);    % dim(B)  = 2n x 4
                M = [phi'*dNx,  -phi'*dNy;
                     0,          2*phi'*dNx;
                     2*phi'*dNy, 0;
                     -phi'*dNx,  phi'*dNy];                
                intFel = Bt*M*dN'*phi;
                % integration of K and F through Gauss quadrature
                fel_y = fel_y + we*intFel(5:8,1)*detJ;
            end
            fel_y = 1/2 * self.Material.permittivity * ...
                self.thickness * fel_y;
        end

    end
    
    methods (Static)
        
        function N = shape_functions(X)
            % N = shape_functions(g,h,r)
            % SHAPE FUNCTIONS FOR A 8-NODED QUADRILATERAL
            % - see Abaqus documentation:
            % Abaqus theory guide > Elements > Continuum elements > ...
            % ... > 3.2.4 Solid isoparametric quadrilaterals and hexahedra
            g = X(1);
            h = X(2);
            N = 1/4*[ +(g - 1)*(h - 1); 
                      -(g + 1)*(h - 1);
                      +(g + 1)*(h + 1); 
                      -(g - 1)*(h + 1)];
        end
        
        function X = natural_coordinates
            X = [-1  -1   % node 1
                  1  -1   % node 2
                  1   1	  % node 3
                 -1   1]; % node 4
        end
        
        function [edges,normal_vectors] = edge_nodes
            edges = [1 2; 2 3; 3 4; 4 1]; % counterclockwise
            normal_vectors = [0 -1; 1 0; 0 1; -1 0];
        end

        function faces = face_nodes
            % face nodes are ordered ccw according to the right hand rule,
            % with the thumb oriented as the outward normal to the face
            faces = {[1 2 3 4]};
        end
    end
        
end % classdef
    
