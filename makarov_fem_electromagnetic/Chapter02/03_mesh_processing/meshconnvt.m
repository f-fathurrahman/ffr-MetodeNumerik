function [si] = meshconnvt(t)
%   SYNTAX 
%   [si] = meshconnvt(t)
%   DESCRIPTION 
%   This function returns triangles attached to every vertex (a cell array si) 
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

   N  = max(max(t));
    si = cell(N, 1);
    for m = 1:N
        temp  = (t(:,1)==m)|(t(:,2)==m)|(t(:,3)==m);
        si{m} = find(temp>0);
    end 
end

