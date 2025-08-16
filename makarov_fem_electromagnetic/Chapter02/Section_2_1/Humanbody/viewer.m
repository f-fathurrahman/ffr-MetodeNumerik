%   SYNTAX
%   viewer
%   DESCRIPTION
%   This script displays a mesh from a *.mat P-t file
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open');
load(FileName, '-mat');

X = reshape(P(t', 1),[3, size(t, 1)]);
Y = reshape(P(t', 2),[3, size(t, 1)]);
Z = reshape(P(t', 3),[3, size(t, 1)]);
patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);

axis 'equal'; axis 'tight', set(gca, 'YDir', 'normal');
camproj('orthographic'); light('Position', [1 3 2], 'Style', 'local'); 
material shiny;
view(70, 25); grid on; xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');