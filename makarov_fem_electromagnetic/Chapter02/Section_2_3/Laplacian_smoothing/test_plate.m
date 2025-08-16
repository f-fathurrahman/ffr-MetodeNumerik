clear all;
%   SYNTAX
%   test_plate
%   DESCRIPTION
%   This scripts evaluates function PLATE
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

L       = 1;    %   Plate length, m
W       = 1;    %   Plate width, m
Tr      = 50;   %   Approximate number of triangles 
iter    = 25;   %   Number of iterations
[P, t] = plate(L, W, Tr, iter);
size(t, 1)

