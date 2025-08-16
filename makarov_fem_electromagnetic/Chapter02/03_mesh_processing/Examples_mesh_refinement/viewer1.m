%   SYNTAX
%   viewer1
%   DESCRIPTION
%   This script displays a mesh from a *.mat P-t file and major mesh
%   parameters: the number of triangles, minimum triangle quality, and
%   minimum edge length
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open'); load(FileName, '-mat');
X = reshape(P(t', 1),[3, size(t, 1)]);
Y = reshape(P(t', 2),[3, size(t, 1)]);
Z = reshape(P(t', 3),[3, size(t, 1)]);
patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);

light('Position', [0 1 0], 'Style', 'local'); axis 'equal';  axis 'tight'; view(0, 0); grid on;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
    
NumberOfTirangles = size(t, 1)
Q = min(simpqual(P, t))
edges = meshconnee(t);
temp = P(edges(:, 1), :) - P(edges(:, 2), :);
minedgelength = min(sqrt(dot(temp, temp, 2)))

    





