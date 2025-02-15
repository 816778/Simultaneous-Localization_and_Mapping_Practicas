function H = SINGLES (prediction, observations, compatibility)
%-------------------------------------------------------
% University of Zaragoza
% Authors:  J. Neira, J. Tardos
%-------------------------------------------------------
%-------------------------------------------------------
global chi2;
global configuration;

H = zeros(1, observations.m);

% You have observations.m observations, and prediction.n
% predicted features.
%
% For every observation i, check whether it has only one neighbour,
% say feature j, and whether that feature j  has only that one neighbour
% observation i.  If so, H(i) = j.
%
% You will need to check the compatibility.ic matrix
% for this:
%
% compatibility.ic(i,j) = 1 if observation i is a neighbour of
% feature j.

for i = 1:observations.m
    neighbors_i = find(compatibility.ic(i, :) == 1);
    if isscalar(neighbors_i)
        j = neighbors_i(1);
        neighbors_j = find(compatibility.ic(:, j) == 1);
        if isscalar(neighbors_j) && neighbors_j(1) == i
            H(i) = j;
        end
    end
end
configuration.name = 'SINGLES';
