function edges = meshconnee(t) 
%   SYNTAX 
%   edges = meshconnee(t) 
%   DESCRIPTION 
%   This function establishes all edges of a triangular mesh given its
%   t-array. The result will be sorted. Original code - P.-O. Persson
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    edges      = [t(:,[1, 2]); t(:,[1, 3]); t(:,[2, 3])];         %  All edges duplicated
    edges      = unique(sort(edges, 2), 'rows');                  %  Unique edges as node pairs (sorted)
end
   
    
