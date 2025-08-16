clear all;
%   SYNTAX 
%   stlexport
%   DESCRIPTION 
%   A short script to create an stl file from MATLAB arrays P, t
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

[FileName, PathName] = uigetfile('*.mat','Select the mesh file');
load(FileName);
%   write stl file
fv.vertices = P;
fv.faces = t;
NewName =  strcat(FileName(1:end-4), '.stl');
stlwrite(NewName, fv);
