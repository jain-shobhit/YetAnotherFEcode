function [mapFea2To, mapTo2Fea] = map_fea_2_to(xyFea)
% Map the FEA domain to the TO domain and vice versa.
% Inputs:
%   xyFea: matrix of size n x 2, where n is the number of elements.
%          The columns are the x and y coordinates of the element centroids.
% Outputs:
%   mapFea2To: vector of size n x 1, the mapping from the FEA domain to the TO domain.
%   mapTo2Fea: vector of size n x 1, the mapping from the TO domain to the FEA domain.

% Very important to have a correct significant
xyFea = round(xyFea, 5, "significant");
    
% Obtain the "mapFea2To"  
% Sort the element centroids in FEA domain to map those in TO domain
% It is assumed that the order of element numbering in TO domain is followed 
% by X (1st column) -> Y (2nd column)

% Sort the element centroids in FEA domain by X and then by Y
[xyTo_by_x, mapFea2To_by_x] = sortrows(xyFea, 1);
[xyTo_by_xy, mapFea2To_by_xy] = sortrows(xyTo_by_x, 2);

% Create the mapping from FEA domain to TO domain
xyTo = xyTo_by_xy;
mapFea2To = mapFea2To_by_x(mapFea2To_by_xy);

% Check if the sorting is correct or not
if ~(isequal(xyFea(mapFea2To, :), xyTo))
    error('xyFea(mapFea2To,:) ~= xyTo');
end

% Obtain the "mapTo2Fea"
[~, mapTo2Fea] = sort(mapFea2To);

% Check if the sorting is correct or not
if ~(isequal(xyTo(mapTo2Fea, :), xyFea))
    error('xyTo(mapTo2Fea,:) ~= xyFea');
end
