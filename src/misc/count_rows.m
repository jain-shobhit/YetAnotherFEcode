function [c, Mu, ic] = count_rows(M)
% this function counts the unique rows of a matrix. Credits to Star Strider
% https://it.mathworks.com/matlabcentral/answers/282561-count-numbers-of-identical-rows
    [Mu,~,ic] = unique(M,'rows','stable'); % Unique Values By Row
    h = accumarray(ic, 1);                 % Count Occurrences
    c = h(ic);                             % Map Occurrences To ic’ Values
end