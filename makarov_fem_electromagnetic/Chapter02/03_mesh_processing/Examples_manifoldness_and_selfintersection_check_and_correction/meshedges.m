function [NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t)
%   SYNTAX 
%   [NonManifoldAttached, edges, edgesb, edgesnm] = meshedges(P, t)
%   DESCRIPTION 
%   This function outputs all edges (array edges), boundary edges (array
%   edgesb), and non-manifold edges (array edgesnm) of a mesh P, t. It also
%   output triangles attached to the boundary edges (array
%   NonManifoldAttached). This function is mostly redundant, it may be
%   replaced by meshconnee, meshconnet

%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    count = 0;
    NonManifoldAttached = [];
    edges = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];         % Interior edges duplicated
    edges = unique(sort(edges, 2), 'rows');               % Edges as node pairs 
    edgesb = [];
    edgesnm = [];
    EDGES = size(edges, 1);
    for m = 1:EDGES
        ind1 = (t(:, 1)==edges(m, 1));    
        ind2 = (t(:, 2)==edges(m, 1));
        ind3 = (t(:, 3)==edges(m, 1));    
        IND1 = ind1|ind2|ind3;          % Index into trianges that include the first edge vertex
        ind1 = (t(:, 1)==edges(m, 2));    
        ind2 = (t(:, 2)==edges(m, 2));
        ind3 = (t(:, 3)==edges(m, 2));    
        IND2 = ind1|ind2|ind3;         % Index into trianges that include the second edge vertex
        IND  = find(IND1&IND2);        % Index into triangles that include the given edge   
        if size(IND, 1)==1             % Boundary triangle found 
           edgesb = [edgesb; edges(m, :)];
           NonManifoldAttached = [NonManifoldAttached IND];
        end
        if size(IND, 1)>2              % Non-manifold edge found with three or more attached trianges
          count = count + 1;
          warning('Non-manifold edges found');
          edgesnm = [edgesnm; edges(m, :)];
        end
    end
    NonManifoldAttached = NonManifoldAttached';
end
