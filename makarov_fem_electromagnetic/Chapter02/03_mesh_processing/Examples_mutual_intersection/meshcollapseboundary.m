function [Pref1, tref1, Pref2, tref2, steps] = meshcollapseboundary(Pref1, tref1, Pref2, tref2, minlength);
%   SYNTAX 
%   [Pref1, tref1, Pref2, tref2, steps] = meshcollapseboundary(Pref1, tref1, Pref2, tref2, minlength);
%   DESCRIPTION 
%   This function performs edge collapse at the boundary between two
%   intersecting meshes until the minimum edge length at the boundary
%   approaches minlength Outputs two reduced meshes. No extra operations
%   implied
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    steps = 0;
    while 1
        %   Find all common nodes of node sets P1 and P2
        [nodes1, nodes2] = checknodes(Pref1, Pref2);    
        %   Find all edges and attached triangles of mesh #1 
        edges1              = meshconnee(tref1);
        AttachedTriangles1  = meshconnet(tref1, edges1, 'manifold');                     
        %   Find all edges and attached triangles of mesh #2 
        edges2              = meshconnee(tref2);
        AttachedTriangles2  = meshconnet(tref2, edges2, 'manifold');       
        %   Find all boundary edges using mesh #1 (mesh #2 will have the same)    
        iedges = zeros(size(edges1)); 
        for n = 1:size(nodes1, 2)
            index = find(edges1(:, 1)==nodes1(n));
            iedges(index, 1) = n;
            index = find(edges1(:, 2)==nodes1(n));
            iedges(index, 2) = n;
        end
        index1 = find(iedges(:, 1)>0&iedges(:, 2)>0); 
        %   Find lengths of such edges and the minimum edge
        temp    = Pref1(edges1(index1, 1), :) - Pref1(edges1(index1, 2), :);    
        [length1, ind1] = min(sqrt(dot(temp, temp, 2))); 
        %   Check if we need to stop 
        if  length1 > minlength 
            break; 
        end;   
        %   Find the corresponding boundary edge of mesh #2
        edge1   = index1(ind1);               %   global number of min edge1
        common  = iedges(edge1, :);           %   indexes into common nodes (both sets)
        temp    = sort(nodes2(common));
        edge2   = find(edges2(:, 1)==temp(1)&(edges2(:, 2)==temp(2)));
        %   Perform simultaneous edge collapse toward edge center        
        [Pref1, tref1] = meshcollapse(Pref1, tref1, edges1(edge1, :), AttachedTriangles1(edge1, :), 3);
        [Pref2, tref2] = meshcollapse(Pref2, tref2, edges2(edge2, :), AttachedTriangles2(edge2, :), 3); 
        steps = steps + 1;
    end
end

