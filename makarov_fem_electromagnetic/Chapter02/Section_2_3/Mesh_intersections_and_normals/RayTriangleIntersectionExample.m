clear all;
%   SYNTAX
%   RayTriangleIntersectionExample
%   DESCRIPTION
%   This script computes and displays intersection points of a ray with the
%   triangular mesh The mesh must be 2 manifold
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

load VHP_eye_left;

orig0 = [20 -50 3];                             % Origin of the ray
dir0 = [25  -17   39];                          % Direction of the ray
dir0 = dir0/norm(dir0);                         % Normalized direction

vert1 = P(t(:, 1),:);
vert2 = P(t(:, 2),:);
vert3 = P(t(:, 3),:);

% Repmat dimensions uniform
orig  = repmat(orig0, size(vert1, 1),1);
dir   = repmat(dir0, size(vert1, 1),1);

[t1] = RayTriangleIntersection(orig, dir, vert1, vert2, vert3);

intersection = orig + t1*dir0;      % Find intersection points
intersection = unique(intersection,'rows');

%% Plot all patches of the first mesh
color = [.8, .9, 1];
patch('Faces', t, 'Vertices', P, 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'FaceColor', color);
axis 'equal';
xlabel('x, mm'); ylabel('y, mm');  grid on;

hold on;
% origin
text(orig0(1), orig0(2), orig0(3), 'orig');
plot3(orig0(1), orig0(2), orig0(3), 'k.', 'MarkerSize', 15);

% direction
quiver3(orig0(1), orig0(2), orig0(3), dir0(1), dir0(2), dir0(3), 60);

% intersection
plot3(intersection(:,1), intersection(:,2), intersection(:,3), 'r.', 'MarkerSize', 25);
view(-145, 44)
hold off; 

