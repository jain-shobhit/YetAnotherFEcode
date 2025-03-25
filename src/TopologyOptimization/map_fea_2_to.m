function [mapFea2To, mapTo2Fea] = map_fea_2_to(coordFea)
    % Map between the FEA and the TO domains. In the TO domain, the
    % elements are lexicographic sorting according to the coordinates of
    % the elements' centroids.
    % Inputs:
    %   coordFea: matrix of size n x m, where n is the number of elements
    %       and m is the number of geometrical dimensions. Each row
    %       represents the coordinates of the centroid of an element.
    % Outputs:
    %   mapFea2To: the map from the FEA domain to the TO one (n x 1).
    %   mapTo2Fea: the map from the TO domain to the FEA one (n x 1).

    % Dimension of the domain
    nDim = size(coordFea, 2);
    
    % Very important to have a correct significant
    coordFea = round(coordFea, 5, "significant");
        
    % Obtain the map from FEA to TO
    [coordTo, mapFea2To] = sortrows(coordFea, nDim:-1:1);

    % Obtain the map from TO to FEA
    [~, mapTo2Fea] = sort(mapFea2To);
    
    % Check the sorting
    if ~(isequal(coordFea(mapFea2To, :), coordTo))
        error('coordFea(mapFea2To, :) ~= coordTo');
    end
    
    % Check if the sorting is correct or not
    if ~(isequal(coordTo(mapTo2Fea, :), coordFea))
        error('coordTo(mapTo2Fea, :) ~= coordFea');
    end
end