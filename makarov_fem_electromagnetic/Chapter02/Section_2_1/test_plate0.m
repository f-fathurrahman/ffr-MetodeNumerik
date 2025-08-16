%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear variables; close all;
L = 0.3; 
W = 0.01; 
Nx = 16; 
Ny = 3;
Tr = 400; 
indicator = 1;
nx = 0; ny = 0; nz = 1;

[P, t] = plate0(L, W, Nx, Ny, indicator, nx, ny, nz);
fv.faces = t;
fv.vertices = P;
patch(fv, 'FaceColor', 'y');
axis equal;
view(0, 90);
grid on;
xlabel('x');
ylabel('y'); 

triangles   = size(t, 1);
[q, R]      = simpqual_size(P, t);
min_quality = min(q);
min_size    = min(R);
