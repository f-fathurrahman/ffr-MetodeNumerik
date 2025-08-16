%   SYNTAX
%   circle0
%   DESCRIPTION
%   This script creates and displays Delaunay triangulation for a circle
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear all;
R = 0.5;     %   Circle radius
N = 10;      %   Number of boundary segments

%   Create boundary nodes    
phi = [0:2*pi/N:2*pi*(N-1)/N];
bx  = R*cos(phi); by = R*sin(phi);
P   = [bx; by]';

%   Create and view the mesh 
DT  = delaunayTriangulation(P)
t   =  DT.ConnectivityList; 

fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
axis equal; view(0, 90); grid on;