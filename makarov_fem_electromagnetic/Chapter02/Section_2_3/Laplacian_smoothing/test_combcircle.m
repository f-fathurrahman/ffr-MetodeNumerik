%   SYNTAX
%   test_combcircle
%   DESCRIPTION
%   This scripts evaluates function COMBCIRCLE with multiple enclosed
%   domains based on constrained Delaunay triangulation
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

clear all;
R       = 1;                %   Outer (main) circle radius
Tr      = 800;              %   Approximate number of triangles
iter    = 20;               %   Number of iterations

strint.R = [0.2 0.3 0.2];       %   Structure of internal boundaries - electrode radii
strint.x = [-0.5 0.0 0.5];      %   Structure of internal boundaries - x-positions
strint.y = [-0.5 0.0 0.5];      %   Structure of internal boundaries - y-positions
strint.N = [15   20   15];      %   Structure of internal boundaries - number of divisions
strint.c = ['r'  'g'  'b'];     %   Structure of internal boundaries - domain color

[P, t] = combcircle(R, Tr, iter, strint);

size(t, 1)