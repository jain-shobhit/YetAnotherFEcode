function h = PlotFieldonMesh(Nodes,Elements,c,varargin)
%--------------------------------------------------------------------------
% Synthax:
%   PlotFieldonMesh(Nodes,Elements,c,varargin)
% Description: plot a scalar field on the (undeformed) mesh
% INPUTS:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           c - The scalar field to be plotted
%           -----> c = a column vector in the order of node numbers
% OUTPUT:   h = plotted object handle
%--------------------------------------------------------------------------

meshcolor = parse_inputs(varargin{:});

if iscell(Elements)
    cELEMENTS = Elements;
else
    cELEMENTS{1} = Elements;
end

for ee = 1 : length(cELEMENTS)
    % cycle the PlotMesh function over all element types in the cell array
    Elements = cELEMENTS{ee};

    nnodes = size(Nodes,1);      % number of nodes
    dimension = size(Nodes,2) ;  % Dimension of the mesh
    elementdim = rank(diff(Nodes(Elements(1,:),:))); % Dimension of elements

    nel = size(Elements,1);      % total number of elements
    nnel = size(Elements,2);     % number of nodes per element

    hold on

    if dimension == 3   % For 3D plots

        if elementdim == 3 % solid in 3D when we simply plot the skin elements
            faces = getSkin3D(Elements);
            Elements = faces.';
            nel = size(Elements,1);      % total number of faces
            nnel = size(Elements,2);     % number of nodes per face
        end

        % Preparing for plot with patch (Elements contains the faces in 2D)
        % X,Y,Z are nnel-by-nel matrices containing x,y,z coordinates for
        % each node in the mesh. E.g., X(i,j) contains the x-coordinate of
        % the i-th node of the j-th element in the mesh, i.e., the node
        % with the index Elements(j,i).

        X = Nodes(Elements',1); X = reshape(X, nnel, nel);
        Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
        Z = Nodes(Elements',3); Z = reshape(Z, nnel, nel);

        profile = c(Elements',1);
        profile = reshape(profile,[nnel length(profile)/nnel]);

        view(3); hold on;
        h = patch(X,Y,Z,profile,'EdgeColor',meshcolor,...
            'DisplayName','Field on Mesh');
        rotate3d on;

    elseif dimension == 2           % For 2D plots

        hold on;

        if elementdim == 2 % surface in 2D

            X = Nodes(Elements',1); X = reshape(X,[nnel nel]);
            Y = Nodes(Elements',2); Y = reshape(Y,[nnel nel]);

            profile = c(Elements',1);
            profile = reshape(profile,[nnel length(profile)/nnel]);

            h{1} = patch(X,Y,profile,'EdgeColor',meshcolor);%,'Marker','o','MarkerFaceColor',meshcolor);
        else
            h = plot(Nodes(:,1),Nodes(:,2),'.-','Color', meshcolor, 'Markersize',10);
            c = [];
        end

    end
    axis equal;
    axis off;
    % Colorbar Setting
    if ~isempty(c)
%         SetColorbar
        colorbar
    end

end
end

function [meshcolor] = parse_inputs(varargin)

defaultColor = 'k';
p = inputParser;

addParameter(p,'color',defaultColor)
parse(p,varargin{:});

meshcolor = p.Results.color;

end









