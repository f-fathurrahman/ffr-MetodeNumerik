function [nodes1, nodes2] = checknodes(P1, P2);
%   SYNTAX 
%   [nodes1, nodes2] = checknodes(P1, P2); 
%   DESCRIPTION 
%   This function returns all common nodes (indexes) of node sets P1 and P2.
%   Both sets should be unique. Tolerance is 1024*eps. Arguments nodes1 and
%   nodes2 are indexes into common nodes.
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

N1 = size(P1, 1);
N2 = size(P2, 1);
nodes1 = [];
nodes2 = []; 

for m = 1:N1
    temp = P2  - repmat(P1(m, :), N2, 1);
    dist = sqrt(dot(temp, temp, 2));
    [Min, indexmin] = min(dist);
    [Max, indexmax] = max(dist);
    if Min/Max<1024*eps
        nodes1 = [nodes1 m];
        nodes2 = [nodes2 indexmin];
    end
end
  