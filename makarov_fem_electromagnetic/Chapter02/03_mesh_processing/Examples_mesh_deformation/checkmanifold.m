function [t, flag] = checkmanifold(t);
%   SYNTAX
%   [t, flag] = checkmanifold(t);
%   DESCRIPTION
%   This function checks mesh holes and non-manifold edges and attempts to
%   correct some non-manifold errors
%   Returns flag (0 - good; 1 - bad) and a corrected array of vertices t
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    edges             = meshconnee(t);
    AttachedTriangles = meshconnet(t, edges, 'nonmanifold');
    flag = 0;
    NonManifoldAttached = [];
    for m = 1:size(edges, 1)
        if length(AttachedTriangles{m})~=2        
            NonManifoldAttached = [NonManifoldAttached; AttachedTriangles{m}];
        end
    end        
    if ~isempty(NonManifoldAttached)       
        %   attempt to correct some errors (eliminate repeated triangles from the list)
        temp = sort(NonManifoldAttached);
        remove = [];
        for m = 2:length(temp)
            if temp(m)==temp(m-1)
                remove = [remove temp(m)];
            end
        end
        t(remove, :) = [];
        edges             = meshconnee(t);
        AttachedTriangles = meshconnet(t, edges, 'nonmanifold');        
        NonManifoldAttached = [];
        for m = 1:size(edges, 1)
            if length(AttachedTriangles{m})~=2     
                NonManifoldAttached = [NonManifoldAttached; AttachedTriangles{m}];
            end
        end   
        if ~isempty(NonManifoldAttached) % if is still bad
            flag = 1;
            %             patch('Faces', t, 'Vertices', P, 'FaceColor', 'y', 'EdgeColor', 'k', 'FaceAlpha', 0.5);  
            %             patch('Faces', t(NonManifoldAttached, :), 'Vertices', P, 'FaceColor', 'r', 'EdgeColor', 'k', 'FaceAlpha', 1.0);   
            %             grid on; axis equal; axis tight;
            %             error('manifoldness check failed')
        else
            display('Non-manifold edges are not found. Boundary edges are not found')
        end
    else       
        display('Non-manifold edges are not found. Boundary edges are not found')
    end    
end