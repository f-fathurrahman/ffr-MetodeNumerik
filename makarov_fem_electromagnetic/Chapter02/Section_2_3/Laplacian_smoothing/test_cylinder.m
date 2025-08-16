clear all;
%   SYNTAX
%   test_cylinder
%   DESCRIPTION
%   This scripts evaluates function CYLINDER
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

R = 1;      %   Cylinder radius
H = 4;      %   Cylinder height
M = 500;    %   Approximate number of triangles
[P, t] = cylinder(R, H, M);
fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
axis equal; view(160, 60); grid on;
xlabel('x'); ylabel('y');
str.iter = 0; str.quality = min(simpqual(P, t)); str
size(t, 1)