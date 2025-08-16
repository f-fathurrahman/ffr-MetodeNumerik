function sj = meshconnvv(t) 
%   SYNTAX 
%   sj = meshconnvv(t) 
%   DESCRIPTION 
%   Given t-array, this function returns vertices directly connected to
%   every vertex - a cell array sj
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    N = max(max(t));    %   Number of nodes in the mesh
    si = cell(N, 1);
    sj = cell(N, 1);
    for m = 1:N    
        temp  = (t1(:,1)==m)|(t1(:,2)==m)|(t1(:,3)==m);
        si{m} = find(temp>0);
        sj{m} = unique(reshape(t(si{m}, :), 1, 3*length(si{m})));
    end 
end
   
    
