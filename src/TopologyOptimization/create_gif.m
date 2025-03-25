function create_gif(nel, densHistory, varargin)
% Create gif of density history.
% Inputs:
%   nel: vector of size 2 x 1, the number of elements in each direction.
%   densHistory: matrix of size n x m, the density history.
%   mapFea2To: vector of size n x 1, the mapping from the FEA domain to the
%       TO domain (optional, default is []).
%   fileName: string, the name of the gif file (optional, default is
%       'densHistory').

% Parse the inputs
p = inputParser;
addOptional(p, 'mapFea2To', []);
addOptional(p, 'fileName', 'densHistory');
parse(p, varargin{:});
mapFea2To = p.Results.mapFea2To;
fileName = p.Results.fileName;

% Number of density vectors
nPoints = find(~isnan(densHistory(1, :)), 1, 'last');

% Create figure
fig = figure('Color', 'w');

% Loop
fileName = [fileName, '.gif'];
delay_time = 0.1;
for ii = 1:nPoints
    % Plot
    plot_layout(nel, densHistory(:, ii), mapFea2To);
    drawnow

    % Save gif
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if ii == 1
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf, 'DelayTime', delay_time, 'BackgroundColor', 0);
    else
        imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end
