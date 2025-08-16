%   SYNTAX
%   constrained
%   DESCRIPTION
%   This script creates and displays constrained Delaunay triangulation for
%   a circle
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear all;

R = 0.5;     %   Circle radius
N = 12;      %   Number of boundary segments

%   Create rectangle nodes inside the circle
P0 = [-R/2 -R/2;...
      -R/2 +R/2;...
      +R/2 +R/2;...
      +R/2 -R/2;...
      +0.0 +0.0];
  
%   Create rectangle edges to be kept 
C = [1 2;...
     2 3;...
     3 4;...
     4 1];
 
%   Create boundary nodes for the circle   
phi = [0:2*pi/N:2*pi*(N-1)/N];
bx  = R*cos(phi); by = R*sin(phi);
P   = [bx; by]';

%   Combine meshes together
P = [P0; P];

%   Create and view the resulting mesh 
DT  = delaunayTriangulation(P, C);
t   =  DT.ConnectivityList;
t(isInterior(DT), :) = [];

fv.faces = t; fv.vertices = P; patch(fv, 'FaceColor', 'y'); 
axis equal; view(0, 90); grid on;