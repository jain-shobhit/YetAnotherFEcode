classdef DpromAssembly < ReducedAssembly
    properties
        % Already contains properties in the Assembly&ReducedAssembly class        
        U       % Defect basis
    end

    methods
        
        function self = DpromAssembly(Mesh, U, V)
            if nargin < 3
                V = [];
                warning(' Reduction Basis V not defined yet')
            end
            % use ReducedAssembly (superclass) constructor
            self@ReducedAssembly(Mesh, V);
            self.DATA.L1 = L1_matrix(self);
            self.DATA.L2 = L2_matrix(self);
            self.DATA.L3 = L3_matrix(self);
            self.DATA.A2 = A2_fun(self);
            self.DATA.A3 = A3_fun(self);
            self.U = U;                     % set defect basis
        end
        
        function dKdp = stiffness_defect_derivative(self, x, formulation)
            dKdp = self.matrix('stiffness_defect_derivative', x, formulation);
        end
        
        function T = Qtensors(self, formulation, volume, varargin)
            % This function computes the tensors for DpROM.
            
            t0 = tic;            
            nv = size(self.V, 2);
            nu = size(self.U, 2);
            Elements = self.Mesh.Elements;
            data = self.DATA;
            
            if volume == 1
                Nend = nu + 1;
                fname = ['Qten_' upper(formulation) 'd'];
            else
                Nend = 1;
                fname = ['Qten_' upper(formulation) 'n'];
            end
            
            % initialize tensors (9 tensors + 1 for each defect)
            for d = 1 : Nend
                T.Q2n{d}  = zeros(nv,nv);
                T.Q3d{d}  = zeros(nv,nv,nu);
                T.Q4dd{d} = zeros(nv,nv,nu,nu);
                T.Q3n{d}  = zeros(nv,nv,nv);
                T.Q4d{d}  = zeros(nv,nv,nv,nu);
                T.Q5dd{d} = zeros(nv,nv,nv,nu,nu);
                T.Q4n{d}  = zeros(nv,nv,nv,nv);
                T.Q5d{d}  = zeros(nv,nv,nv,nv,nu);
                T.Q6dd{d} = zeros(nv,nv,nv,nv,nu,nu);
            end
            data.Qinit = T;
            
            % Computing element level contributions
            for e = 1 : length(Elements)
                thisElement = Elements(e).Object;
                index = thisElement.iDOFs;          
                Ve = self.V(index, :);
                Ue = self.U(index, :);
                Te = thisElement.DpROM_version(Ve, Ue, data, fname);
                for d = 1 : Nend
                    T.Q2n{d}  = T.Q2n{d}  + Te.Q2n{d};
                    T.Q3d{d}  = T.Q3d{d}  + Te.Q3d{d};
                    T.Q4dd{d} = T.Q4dd{d} + Te.Q4dd{d};
                    T.Q3n{d}  = T.Q3n{d}  + Te.Q3n{d};
                    T.Q4d{d}  = T.Q4d{d}  + Te.Q4d{d};
                    T.Q5dd{d} = T.Q5dd{d} + Te.Q5dd{d};
                    T.Q4n{d}  = T.Q4n{d}  + Te.Q4n{d};
                    T.Q5d{d}  = T.Q5d{d}  + Te.Q5d{d};
                    T.Q6dd{d} = T.Q6dd{d} + Te.Q6dd{d};
                end
            end
            
            % convert into tensors
            for d = 1 : Nend
                T.Q3d{d}  = tensor(T.Q3d{d});
                T.Q4dd{d} = tensor(T.Q4dd{d});
                T.Q3n{d}  = tensor(T.Q3n{d});
                T.Q4d{d}  = tensor(T.Q4d{d});
                T.Q5dd{d} = tensor(T.Q5dd{d});
                T.Q4n{d}  = tensor(T.Q4n{d});
                T.Q5d{d}  = tensor(T.Q5d{d});
                T.Q6dd{d} = tensor(T.Q6dd{d});
            end
            
            T.formulation = formulation;
            T.volume = volume;
            T.time = toc(t0);
        end
        
        function M = ParametricMass(self, varargin)
            % This function computes the reduced order mass matrix for the
            % DpROM, approximating the integral over the defected volume
            % (same procedure used for the stiffness tensors).
            % The defected mass matrix can be obtained as:
            %   Md = M{1} + M{2}*xi(1) + M{3}*xi(2) + ... M{nu+1}*xi(nu)
            
            nv = size(self.V, 2);
            nu = size(self.U, 2);
            Elements = self.Mesh.Elements;
            Nend = nu + 1;
            
            % initialize
            for d = 1 : Nend
                Minit{d}  = zeros(nv,nv); %#ok<AGROW>
            end
            M = Minit;
            
            % Computing element level contributions
            for e = 1 : length(Elements)
                thisElement = Elements(e).Object;
                index = thisElement.iDOFs;          
                Ve = self.V(index, :);
                Ue = self.U(index, :);
                Me = thisElement.('DpROM_mass')(Ve, Ue, Minit);
                for d = 1 : Nend
                    M{d}  = M{d}  + Me{d};
                end
            end
        end
        
        % Ancillary functions _____________________________________________
        
        function set.U(self,U)
            if size(U,1) ~= self.Mesh.nDOFs %#ok<*MCSUP>
                warning(['Reduction basis size incorrect: should have ' ...
                    num2str(self.Mesh.nDOFs) ' rows'])
            end
            self.U = U;
        end
        
        function A2 = A2_fun(self)
            if self.Mesh.nDim == 2
                A2 = @(th) ...
                    -[th(1) th(3) 0     0;
                      0     0     th(2) th(4);
                      th(2) th(4) th(1) th(3)];
            elseif self.Mesh.nDim == 3
                A2 = @(th) -[ ...
                th(1),  th(4),  th(7),      0,      0,      0,      0,      0,      0;
                    0,      0,      0,  th(2),  th(5),  th(8),      0,      0,      0;
                    0,      0,      0,      0,      0,      0,  th(3),  th(6),  th(9);
                th(2),  th(5),  th(8),  th(1),  th(4),  th(7),      0,      0,      0;
                th(3),  th(6),  th(9),      0,      0,      0,  th(1),  th(4),  th(7);
                    0,      0,      0,  th(3),  th(6),  th(9),  th(2),  th(5),  th(8)];
            end
        end
        
        function A3 = A3_fun(self)
            if self.Mesh.nDim == 2
                A3 = @(th) ...
                    -[th(1),      0,      th(3)/2;
                      0,          th(4),  th(2)/2;
                      th(2),      th(3),  (th(1)+th(4))/2];
            elseif self.Mesh.nDim == 3
                A3 = @(th) -[ ...
                th(1) 0     0     th(4)/2         th(7)/2         0;
                0	  th(5) 0     th(2)/2         0               th(8)/2;
                0	  0	    th(9) 0               th(3)/2         th(6)/2;
                th(2) th(4) 0     (th(1)+th(5))/2 th(8)/2         th(7)/2;
                th(3) 0     th(7) th(6)/2         (th(1)+th(9))/2 th(4)/2;
                0     th(6) th(8) th(3)/2         th(2)/2         (th(5)+th(9))/2];
            end
        end
        
        function L1 = L1_matrix(self)
            % L used to compute the quadratic strain matrix: 
            % A = L.th, eps_quad = A*th ("." is the contraction operation)
            if self.Mesh.nDim == 2
                L = zeros([3,4,4]);
                L(1,1,1) = 1; L(3,2,1) = 1; L(3,1,2) = 1; L(2,2,2) = 1;
                L(1,3,3) = 1; L(3,4,3) = 1; L(3,3,4) = 1; L(2,4,4) = 1;
            elseif self.Mesh.nDim == 3
                L = zeros([6,9,9]);
                L(1,1,1)=1; L(4,2,1)=1; L(5,3,1)=1; 
                L(4,1,2)=1; L(2,2,2)=1; L(6,3,2)=1; 
                L(5,1,3)=1; L(6,2,3)=1; L(3,3,3)=1;
                L(1,4,4)=1; L(4,5,4)=1; L(5,6,4)=1; 
                L(4,4,5)=1; L(2,5,5)=1; L(6,6,5)=1; 
                L(5,4,6)=1; L(6,5,6)=1; L(3,6,6)=1;
                L(1,7,7)=1; L(4,8,7)=1; L(5,9,7)=1; 
                L(4,7,8)=1; L(2,8,8)=1; L(6,9,8)=1; 
                L(5,7,9)=1; L(6,8,9)=1; L(3,9,9)=1;
            end
            L1 = L;
        end
        
        function L2 = L2_matrix(self)
            if self.Mesh.nDim == 2
                L2(1,1,1) = 1; L2(3,3,1) = 1; L2(3,1,2) = 1; L2(2,3,2) = 1;
                L2(1,2,3) = 1; L2(3,4,3) = 1; L2(3,2,4) = 1; L2(2,4,4) = 1;
            elseif self.Mesh.nDim == 3
                L2(1,1,1) = 1; L2(4,4,1) = 1; L2(5,7,1) = 1; 
                L2(4,1,2) = 1; L2(2,4,2) = 1; L2(6,7,2) = 1; 
                L2(5,1,3) = 1; L2(6,4,3) = 1; L2(3,7,3) = 1; 
                L2(1,2,4) = 1; L2(4,5,4) = 1; L2(5,8,4) = 1; 
                L2(4,2,5) = 1; L2(2,5,5) = 1; L2(6,8,5) = 1;
                L2(5,2,6) = 1; L2(6,5,6) = 1; L2(3,8,6) = 1; 
                L2(1,3,7) = 1; L2(4,6,7) = 1; L2(5,9,7) = 1; 
                L2(4,3,8) = 1; L2(2,6,8) = 1; L2(6,9,8) = 1; 
                L2(5,3,9) = 1; L2(6,6,9) = 1; L2(3,9,9) = 1;
            end
            L2 = -L2;
        end
        
        function L3 = L3_matrix(self)
            if self.Mesh.nDim == 2
                L3(1,1,1,1) = 2; L3(3,2,1,1) = 1; L3(3,1,2,1) = 1; L3(1,3,3,1) = 2;
                L3(3,4,3,1) = 1; L3(3,3,4,1) = 1; L3(3,1,1,2) = 2; L3(2,2,1,2) = 1;
                L3(2,1,2,2) = 1; L3(3,3,3,2) = 2; L3(2,4,3,2) = 1; L3(2,3,4,2) = 1;
                L3(1,2,1,3) = 1; L3(1,1,2,3) = 1; L3(3,2,2,3) = 2; L3(1,4,3,3) = 1;
                L3(1,3,4,3) = 1; L3(3,4,4,3) = 2; L3(3,2,1,4) = 1; L3(3,1,2,4) = 1;
                L3(2,2,2,4) = 2; L3(3,4,3,4) = 1; L3(3,3,4,4) = 1; L3(2,4,4,4) = 2;
            elseif self.Mesh.nDim == 3
                L3(1,1,1,1) = 2; L3(4,2,1,1) = 1; L3(5,3,1,1) = 1; L3(4,1,2,1) = 1;
                L3(5,1,3,1) = 1; L3(1,4,4,1) = 2; L3(4,5,4,1) = 1; L3(5,6,4,1) = 1;
                L3(4,4,5,1) = 1; L3(5,4,6,1) = 1; L3(1,7,7,1) = 2; L3(4,8,7,1) = 1;
                L3(5,9,7,1) = 1; L3(4,7,8,1) = 1; L3(5,7,9,1) = 1; L3(4,1,1,2) = 2;
                L3(2,2,1,2) = 1; L3(6,3,1,2) = 1; L3(2,1,2,2) = 1; L3(6,1,3,2) = 1;
                L3(4,4,4,2) = 2; L3(2,5,4,2) = 1; L3(6,6,4,2) = 1; L3(2,4,5,2) = 1;
                L3(6,4,6,2) = 1; L3(4,7,7,2) = 2; L3(2,8,7,2) = 1; L3(6,9,7,2) = 1;
                L3(2,7,8,2) = 1; L3(6,7,9,2) = 1; L3(5,1,1,3) = 2; L3(6,2,1,3) = 1;
                L3(3,3,1,3) = 1; L3(6,1,2,3) = 1; L3(3,1,3,3) = 1; L3(5,4,4,3) = 2;
                L3(6,5,4,3) = 1; L3(3,6,4,3) = 1; L3(6,4,5,3) = 1; L3(3,4,6,3) = 1;
                L3(5,7,7,3) = 2; L3(6,8,7,3) = 1; L3(3,9,7,3) = 1; L3(6,7,8,3) = 1;
                L3(3,7,9,3) = 1; L3(1,2,1,4) = 1; L3(1,1,2,4) = 1; L3(4,2,2,4) = 2;
                L3(5,3,2,4) = 1; L3(5,2,3,4) = 1; L3(1,5,4,4) = 1; L3(1,4,5,4) = 1;
                L3(4,5,5,4) = 2; L3(5,6,5,4) = 1; L3(5,5,6,4) = 1; L3(1,8,7,4) = 1;
                L3(1,7,8,4) = 1; L3(4,8,8,4) = 2; L3(5,9,8,4) = 1; L3(5,8,9,4) = 1;
                L3(4,2,1,5) = 1; L3(4,1,2,5) = 1; L3(2,2,2,5) = 2; L3(6,3,2,5) = 1;
                L3(6,2,3,5) = 1; L3(4,5,4,5) = 1; L3(4,4,5,5) = 1; L3(2,5,5,5) = 2;
                L3(6,6,5,5) = 1; L3(6,5,6,5) = 1; L3(4,8,7,5) = 1; L3(4,7,8,5) = 1;
                L3(2,8,8,5) = 2; L3(6,9,8,5) = 1; L3(6,8,9,5) = 1; L3(5,2,1,6) = 1;
                L3(5,1,2,6) = 1; L3(6,2,2,6) = 2; L3(3,3,2,6) = 1; L3(3,2,3,6) = 1;
                L3(5,5,4,6) = 1; L3(5,4,5,6) = 1; L3(6,5,5,6) = 2; L3(3,6,5,6) = 1;
                L3(3,5,6,6) = 1; L3(5,8,7,6) = 1; L3(5,7,8,6) = 1; L3(6,8,8,6) = 2;
                L3(3,9,8,6) = 1; L3(3,8,9,6) = 1; L3(1,3,1,7) = 1; L3(4,3,2,7) = 1;
                L3(1,1,3,7) = 1; L3(4,2,3,7) = 1; L3(5,3,3,7) = 2; L3(1,6,4,7) = 1;
                L3(4,6,5,7) = 1; L3(1,4,6,7) = 1; L3(4,5,6,7) = 1; L3(5,6,6,7) = 2;
                L3(1,9,7,7) = 1; L3(4,9,8,7) = 1; L3(1,7,9,7) = 1; L3(4,8,9,7) = 1;
                L3(5,9,9,7) = 2; L3(4,3,1,8) = 1; L3(2,3,2,8) = 1; L3(4,1,3,8) = 1;
                L3(2,2,3,8) = 1; L3(6,3,3,8) = 2; L3(4,6,4,8) = 1; L3(2,6,5,8) = 1;
                L3(4,4,6,8) = 1; L3(2,5,6,8) = 1; L3(6,6,6,8) = 2; L3(4,9,7,8) = 1;
                L3(2,9,8,8) = 1; L3(4,7,9,8) = 1; L3(2,8,9,8) = 1; L3(6,9,9,8) = 2;
                L3(5,3,1,9) = 1; L3(6,3,2,9) = 1; L3(5,1,3,9) = 1; L3(6,2,3,9) = 1;
                L3(3,3,3,9) = 2; L3(5,6,4,9) = 1; L3(6,6,5,9) = 1; L3(5,4,6,9) = 1;
                L3(6,5,6,9) = 1; L3(3,6,6,9) = 2; L3(5,9,7,9) = 1; L3(6,9,8,9) = 1;
                L3(5,7,9,9) = 1; L3(6,8,9,9) = 1; L3(3,9,9,9) = 2;
            end
            L3 = -L3/2;
        end                
                
    end
end
