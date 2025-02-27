function [ground, people] = generate_cloister_experiment,
%-------------------------------------------------------
% University of Zaragoza
% Authors:  J. Neira, J. Tardos
%-------------------------------------------------------
% Generate experiment
%-------------------------------------------------------

s = 2;
p = [];
l1 = 10; w = -1.2; d = 0; p = [p genpoints(l1, s, w, d)];
l1 = 10; w = -0.7; d = 0; p = [p genpoints(l1, s, w, d)];
l2 = l1 - s; w = 0.7; d = s/2 ; p = [p genpoints(l2, s, w, d)];
l2 = l1 - s; w = 1.2; d = s/2 ; p = [p genpoints(l2, s, w, d)];
q = tpcomp([0 0 pi/2], p);
q(1,:) = q(1,:) + l1;
r = tpcomp([0 0 pi], [p q]);
r = r + l1;
p = [p q r];

ground.n = size(p, 2);
ground.points = p;

simulate_people(10, 0, l1, 0, l1);

ground.trajectory(1).x = [0 0 0]';
ground.trajectory(1).P = zeros(3, 3);

p = 0;
q = 0;
st = 0.3125;
sa = 0.1;
steps = 4*l1/st + 4/sa;
for step = 1 : steps,
    if (p >= l1)
        q = q + sa*pi/2;
        ground.motion(step).x = [0 0 sa*pi/2]';
        ground.motion(step).P = diag([0.03 0.03 2.5*sa*pi/180].^2);
        if q >= pi/2
           p = 0;
           q = 0;
        end
    else
        p = p + st;
        ground.motion(step).x = [st 0 0]';
        ground.motion(step).P = diag(sqrt(st)*[0.01 0.01 2*pi/180].^2);
    end
    [x, y] = simulate_people;
    people.gx(step,:) = x';
    people.gy(step,:) = y';
end

function points = genpoints (l, s, w, d),
x = [(0 : s : l)] + d;
y = [(w/2)*ones(1, length(x))];

points = [x ; y];
