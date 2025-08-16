%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear variables; close all;
L = 1; 
W = 1.0; 
H = 0.2;
Nx = 7; 
Ny = 7; 
Nz = 3; 
indicator = 1;

[P, t] = brick0(L, W, H, Nx, Ny, Nz, indicator);

fv.faces = t;
fv.vertices = P;
patch(fv, 'FaceColor', 'y');
axis equal;
view(40, 30);
grid on;
xlabel('x');
ylabel('y');

triangles   = size(t, 1);
[q, R]      = simpqual_size(P, t);
min_quality = min(q);
min_size    = min(R);