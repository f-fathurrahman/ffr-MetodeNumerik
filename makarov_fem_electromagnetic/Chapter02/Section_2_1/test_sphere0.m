%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear all;
R = 1;              %   sphere radius
Tr = 2000;          %   approximate number of triangles
[P, t] = sphere0(R, Tr);

fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); axis equal; view(40, 30); grid on;
xlabel('x'); ylabel('y')

triangles   = size(t, 1)
[q, R]      = simpqual_size(P, t);
min_quality = min(q)
min_size    = min(R)