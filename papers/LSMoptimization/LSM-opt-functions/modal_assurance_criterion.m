function [macIndex, macValue] = modal_assurance_criterion(phi, phiRef, doPrint)
% Calculate the Modal Assurance Criterion (MAC) between a set of modes and a reference mode.
% Inputs:
%   phi: matrix of size n x m, where n is the number of DOFs and m is the number of modes.
%        The columns are the modes.
%   phiRef: vector of size n x 1, the reference mode.
%   doPrint: boolean, print the MAC value and index. Default is false.
% Outputs:
%   macIndex: scalar, the index of the mode with the maximum MAC value.
%   macValue: scalar, the MAC value.

% Check inputs
if nargin == 2
    doPrint = false;
end

% Initialize variables
phiDot = dot(phiRef, phiRef);
macValues = zeros(size(phi, 2), 1);

% Loop over modes
for ii = 1:size(phi, 2)
    thisPhi = phi(:, ii);
    macValues(ii) = dot(thisPhi, phiRef).^2 / (phiDot * dot(thisPhi, thisPhi));
end

% Find max mac
[macValue, macIndex] = max(macValues);

% Print
if doPrint
    fprintf("\nIdx %d has MAC %.4f", macIndex, macValue)
end
