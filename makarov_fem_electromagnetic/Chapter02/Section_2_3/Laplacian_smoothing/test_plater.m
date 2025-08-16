clear all;
%   SYNTAX
%   test_plater
%   DESCRIPTION
%   This scripts evaluates function PLATER
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

L       = 1;    %   Plate length, m
W       = 1;    %   Plate width, m
Tr      = 200;  %   Approximate number of triangles 
iter    = 20;   %   Number of iterations
[P, t] = plater(L, W, Tr, iter);
triangles = size(t, 1)

