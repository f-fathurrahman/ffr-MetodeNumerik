function nonmanifold = checkmanifoldvertices(t)
%   SYNTAX
%   nonmanifold = checkmanifoldvertices(t);
%   DESCRIPTION
%   This function checks non-manifold nodes (for the meshes with manifold edges)
%   Returns parameter nonmanifold - indexes into non-manifold nodes
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    N = max(max(t));
    edges               = meshconnee(t);
    AttachedTriangles   = meshconnet(t, edges, 'nonmanifold');
    %   Find triangles attached to every vertex (cell array)
    si                  = meshconnvt(t);
    %   Find edges attached to every vertex (cell array)
    [se] = meshconnve(t, edges);
    %   create pairs of attached triangles and collect non-manifold nodes
    nonmanifold = [];
    for m = 1:N
        E = length(se{m});
        edgeslocal = se{m};
        AttTri = zeros(E, 2);
        for n = 1:E
            AttTri(n, :) = AttachedTriangles{edgeslocal(n)};
        end
        A = AttTri(1, 1); 
        B = AttTri(1, 2);
        AttTri(1, :) = [];
        while 1         
            temp1 = AttTri(:, 1)==A|AttTri(:, 2)==A;
            A = setdiff(AttTri(temp1, :), A); 
            temp2 = AttTri(:, 1)==B|AttTri(:, 2)==B;
            B = setdiff(AttTri(temp2, :), B); 
            if isempty(A)|isempty(B)
                nonmanifold = [nonmanifold m]
                break;
            end            
            out = unique([find(temp1) find(temp2)]);
            AttTri(out, :) = [];
            if isempty(AttTri)
                break;
            end
        end
    end
end
    
    
    
    
   