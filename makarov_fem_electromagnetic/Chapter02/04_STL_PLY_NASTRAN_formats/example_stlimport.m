clear all;
%   SYNTAX 
%   stlimport
%   DESCRIPTION 
%   A short script to import an stl file into MATLAB
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara

[FileName, PathName] = uigetfile('*.stl','Select the mesh file');
[t, P] = stlread(FileName);

X = reshape(P(t', 1),[3, size(t, 1)]);
Y = reshape(P(t', 2),[3, size(t, 1)]);
Z = reshape(P(t', 3),[3, size(t, 1)]);
patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);
light('Position', [0 1 0], 'Style', 'local'); axis 'equal';  axis 'tight'; view(-114, -78); grid on;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');
   