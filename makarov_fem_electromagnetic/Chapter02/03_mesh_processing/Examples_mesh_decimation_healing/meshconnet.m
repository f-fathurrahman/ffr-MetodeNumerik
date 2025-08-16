function AttachedTriangles = meshconnet(t, edges, flag)
%   SYNTAX 
%   AttachedTriangles = meshconnet(t, edges, flag)
%   DESCRIPTION 
%   This function finds triangles attached to every edge of a triangular mesh
%   Inputs: array of nodes P, triangulation t, flag
%   Output: numerical (manifold meshes, flag='manifold') or cell
%   (presumably non-manifold meshes) array of triangles attached to each
%   edge. Array size is (N, 2) for manifolds
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    EDGES = size(edges, 1);  temp = cell(EDGES, 1); 
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
        temp{m}    = IND;
    end
    if strcmp(flag, 'manifold')
        AttachedTriangles = cell2mat(temp')';
    else
       AttachedTriangles = temp; 
    end
end