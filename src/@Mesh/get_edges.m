function edges = get_edges(self, S, plotfigure)

nodes = self.nodes;
Ne = self.nElements;
nNodes = self.Elements(1).Object.nNodes;
EDGES = self.Elements(1).Object.edge_nodes; % TO DO: handle multi-elements

edge_nodes  = zeros(size(EDGES,1)*Ne, size(EDGES,2));
edge_normal = zeros(size(EDGES,1)*Ne, 2);
elements = zeros(Ne, nNodes);
kk = 1;
for ii = 1 : Ne % cycle over elements
    ElemIDs = self.Elements(ii).Object.nodeIDs;
    elements(ii,:) = ElemIDs;
    for jj = 1 : size(EDGES,1) % cycle over edges
        edge_nodes(kk,:) = ElemIDs(EDGES(jj,:)); % store node IDs per edge
        
        x = nodes(edge_nodes(kk,:), 1);
        y = nodes(edge_nodes(kk,:), 2);
        dx = x(2) - x(1);
        dy = y(2) - y(1);
        n = [dy, -dx];
        edge_normal(kk,:) = n/norm(n);
        kk = kk+1;
    end
end
edge_element = kron((1 : Ne)', ones(size(EDGES,1),1)); % element IDs to which the edges belong

edge_nodes_sorted = sort(edge_nodes,2); % [1 2] and [2 1] are the same edge
count = count_rows(edge_nodes_sorted); % count each edge appearances

edge_nodes_unique = count==1; % index to edges appearing only once

edges.order = EDGES;
edges.nodeID = edge_nodes(edge_nodes_unique,:);
edges.elementID = edge_element(edge_nodes_unique,:);
edges.normal = edge_normal(edge_nodes_unique, :);

self.DATA.edges = edges;

if plotfigure == 1
    figure
    PlotMesh(nodes, elements(:,1:4), 0);
    
    uedge = edges.nodeID;
    for ii = 1 : size(uedge,1)
        x = nodes(uedge(ii,:), 1);
        y = nodes(uedge(ii,:), 2);
        plot(x, y, 'r-o', 'LineWidth', 2)
        xn = mean(x);
        yn = mean(y);
        quiver(xn,yn,edges.normal(ii,1)*S,edges.normal(ii,2)*S,0,'k','MaxHeadSize',1)
    end
end

end