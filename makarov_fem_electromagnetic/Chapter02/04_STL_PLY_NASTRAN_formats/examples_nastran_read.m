clear all; close all;
%   SYNTAX 
%   examples_nastran_read
%   DESCRIPTION 
%   The NASTRAN format is the oldest ASCII format for saving triangular
%   meshes. The present script converts a NASTRAN file to the P, t mesh
%   data in MATLAB
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

FileName    = uigetfile('*.nas', 'Select the mesh file to open');
fileID      = fopen(FileName, 'r');
A           = fscanf(fileID, '%49c', [49, inf]); B = A';
%   Extract P and t arrays
P = []; t = [];
for m = 1:size(B, 1)
    Nodes = strfind(B(m, :), 'GRID');
    if ~isempty(Nodes) 
        P = [P; [str2num(B(m, 25:32)) str2num(B(m, 33:40)) str2num(B(m, 41:48))]];
    end
    Tri  = strfind(B(m, :), 'CTRIA3');
    if ~isempty(Tri)
        t = [t; str2num(B(m, 25:end))];
    end
end

fclose(fileID);
save(strcat(FileName(1:end-4), '.mat'), 'P', 't');
