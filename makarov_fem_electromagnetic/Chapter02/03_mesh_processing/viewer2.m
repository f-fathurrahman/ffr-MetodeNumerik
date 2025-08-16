clear all
%   SYNTAX
%   viewer2
%   DESCRIPTION
%   This script displays multiple meshes from *.mat P-t files and major mesh
%   parameters: number of triangles, minimum triangle quality, and
%   minimum edge length
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.

FileName = uigetfile('*.mat','Select the tissue mesh file to open', 'MultiSelect', 'on');
if iscell(FileName)
    for m = 1:length(FileName)
        S = load(FileName{m}, '-mat');
        P = S.P; t = S.t; 
        graphics(P, t);
    end
else
    S = load(FileName, '-mat');
    P = S.P; t = S.t; 
    graphics(P, t);
end


  





