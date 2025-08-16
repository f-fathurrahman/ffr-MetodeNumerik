%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear all;
L = 0.3; 
W = 0.1; 
l = L/20;
w = W/10;
indicator = 1;

[P, t] = polygon0(L, W, l, w, indicator);
fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); axis equal; view(0, 90); grid on;
xlabel('x'); ylabel('y'); 

triangles   = size(t, 1)
[q, R]      = simpqual_size(P, t);
min_quality = min(q)
min_size    = min(R)
