function [P, t] = meshcollapse(P, t, edge, AttachedTriangles, No)
%   SYNTAX 
%   [P, t] = meshcollapse(P, t, edge, AttachedTriangles, No);
%   DESCRIPTION 
%   This function performs single elementary edge collapse (toward either
%   node or toward edge's center)
%   Inputs: array of nodes P, triangulation t, mesh edge to be removed - 
%   edge(1, 2), a pair of attached triangles - AttachedTriangles(1, 2)
%   No = 1 -> collapses two nodes toward the first edge node
%   No = 2 -> collapses two nodes toward the second edge node
%   No = 3 -> collapses two nodes toward the edge center
%   Output: reduced mesh P, t
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    t(AttachedTriangles, :) = [];
    if No==1    
        t(t==edge(2)) = edge(1);  
        t(t>edge(2))  = t(t>edge(2))-1; 
        P(edge(2), :) = [];
    end
    if No==2    
        t(t==edge(1)) = edge(2);  
        t(t>edge(1))  = t(t>edge(1))-1; 
        P(edge(1), :) = [];
    end
    if No==3    
        P(edge(1), :) = (P(edge(1), :) + P(edge(2), :))/2; 
        t(t==edge(2)) = edge(1);  
        t(t>edge(2))  = t(t>edge(2))-1; 
        P(edge(2), :) = []; 
    end
end 
