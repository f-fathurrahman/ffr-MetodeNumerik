function se = meshconnve(t, edges)
%   SYNTAX 
%   se = meshconnve(t, edges)
%   DESCRIPTION 
%   This function returns edges (indexes into array edges) attached to
%   every vertex (a cell array se)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.
%   Find edges (indexes) attached to every triangle (Nx3 array)

    N  = max(max(t));
    se = cell(N, 1);
    for m = 1:N
        temp  = (edges(:, 1)==m)|(edges(:, 2)==m);
        se{m} = find(temp>0);
    end
end

