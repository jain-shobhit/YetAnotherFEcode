% meshinfo = abqmesh(filename)
%
% Description: parse an Abaqus .inp file to read nodes and elements.
%
% INPUT: filename (string, with or without extension)
% OUTPUT: meshinfo, a structure with the following fields:
%   .nodes:  nodal coordinates
%   .elem:   cell array collecting the element tables (one per type)
%   .eltype: element type (same order as for .elem)
%   .nset:   node sets
% Note: this function uses regular expressions to parse the .inp file. The
% function can thus easily upgraded to ouput more information, if needed.
%
% Author: Jacopo Marconi, PhD, Politecnico di Milano
% Created: 21 March 2022
% Last modified: 22 March 2022

function meshinfo = abqmesh(filename)

if ~strcmp(filename(end-3:end),'.inp')
    filename = [filename '.inp'];
end

fid = fopen(filename);

a = fscanf(fid, '%c');

% GET NODES
nodes = regexp(a,'\<Node\s\s(.*?)\s\s\*','tokens');
nodes = str2num(nodes{1}{1}); %#ok<*ST2NM> 
meshinfo.nodes = nodes(:, 2:end);

% GET ELEMENTS
b = regexp(a,'\<Element, type=([\w]*)\s\s(.*?)\s\s\*','tokens');

% parse elements
for ii = 1 : length(b)
    eltype{ii} = b{ii}{1};
    elem{ii}   = b{ii}{2};

    lines = splitlines(elem{ii});
    for jj = 1 : length(lines)
        c{jj,:} = str2num(lines{jj});
    end
    
    C = [];
    % check first two lines (elements with many nodes go on two lines)
    if length(c{1}) ~= length(c{2})
        % if different, concatenate lines by pairs:
        jj = 1;
        for kk = 1 : length(c)/2
            C(kk,:) = [c{jj} c{jj+1}];
            jj = jj+2;
        end
    else
        % otherwise, convert to matrix:
        for kk = 1 : length(c)
            C(kk,:) = c{kk};
        end
    end
    % remove first column (IDs)
    elem{ii} = C(:,2:end);
end
meshinfo.elem = elem;
meshinfo.eltype = eltype;

% GET NODE SETS
b = regexp(a,'\<Nset, nset=.*?, instance=.*?\s\s(.*?)\s\s\*','tokens');

% parse node sets
for ii = 1 : length(b)
    lines = splitlines(b{ii}{1});
    nset{ii} = [];
    for jj = 1 : length(lines)
        nset{ii} = [nset{ii} str2num(lines{jj})];
    end
end
meshinfo.nset = nset;

fclose(fid);
end
