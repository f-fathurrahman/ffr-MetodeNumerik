clear all;
%   SYNTAX
%   SegmentTriangleIntersectionExample
%   DESCRIPTION
%   This script computes and displays intersection points of a segment with
%   the triangular mesh. The mesh must be 2 manifold
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

load VHP_eye_left;

orig0 = [20 -50 3];    % Origin of the ray
dest = 1.2*[45 -63 45];  % Destination of the ray
dir0 = dest - orig0;     % Direction of the ray
dist0 = sqrt(dir0*dir0');      % Finding distance
dir0 = dir0/dist0;  % Normalized direction

% Vertices
vert1 = P(t(:, 1),:); 
vert2 = P(t(:, 2),:);
vert3 = P(t(:, 3),:);

% Making dimensions uniform
orig  = repmat(orig0, size(vert1, 1),1);
dir   = repmat(dir0, size(vert1, 1),1);
dist  = repmat(dist0, size(vert1, 1),1);

[t1] = SegmentTriangleIntersection(orig, dir, vert1, vert2, vert3, dist);
t1(t1 ==0) = [];                                        % Avoiding origin points when t1 = 0
t1 = unique(t1);
orig  = repmat(orig0, size(t1, 1),1);
dir   = repmat(dir0, size(t1, 1),1);
intersection = orig + t1*dir0; % Finding intersection points

%% Plot all patches of the first mesh
color = [.8, .9, 1];
patch('Faces', t, 'Vertices', P, 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'FaceColor', color);
axis 'equal';
xlabel('x, mm'); ylabel('y, mm');  grid on;

hold on;
% origin and destination
text(orig0(1) + 1, orig0(2) + 1, orig0(3) + 0.5, 'orig');
plot3(orig0(1), orig0(2), orig0(3), 'k.', 'MarkerSize', 25);
plot3(dest(1), dest(2), dest(3), 'k.', 'MarkerSize', 25);

% direction vector
line = [orig0;dest];
plot3(line(:,1), line(:,2), line(:,3), 'LineWidth', 2);

% intersection
plot3(intersection(:,1), intersection(:,2), intersection(:,3), 'r.', 'MarkerSize', 25);
view(-145, 44)

