%-------------------------------------------------------
function H = JCBB (prediction, observations, compatibility)
% 
%-------------------------------------------------------
global Best;
global configuration;

Best.H = zeros(1, observations.m);

JCBB_R (prediction, observations, compatibility, [], 1);

H = Best.H;
configuration.name = 'JCBB';

%-------------------------------------------------------
function JCBB_R (prediction, observations, compatibility, H, i)
% 
%-------------------------------------------------------
global Best;
global configuration;

if i > observations.m % leaf node?
    if pairings(H) > pairings(Best.H) % did better?
        Best.H = H;
    end
else
    % complete JCBB here
    for j = 1:prediction.n % for j in {1...n}
        % if individual_compatibility(i,j) and then joint_compatibility(H, i, j)
        if compatibility.ic(i, j) && jointly_compatible(prediction, observations, [H j])
            JCBB_R (prediction, observations, compatibility, [H j], i+1); % pairing (E_i, F_j) accepted
        end
    end
    % if pairing(H) + m - i > pairing(Best)
    if pairings(H) + observations.m - i > pairings(Best.H) % can do better?
        JCBB_R (prediction, observations, compatibility, [H 0] , i+1);  % start node, E_i nor paired
    end
end

%-------------------------------------------------------
% 
%-------------------------------------------------------
function p = pairings(H)

p = length(find(H));