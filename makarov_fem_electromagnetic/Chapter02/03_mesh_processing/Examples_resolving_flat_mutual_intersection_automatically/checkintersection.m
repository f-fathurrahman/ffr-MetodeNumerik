function flag = checkintersection(P, t);
%   SYNTAX
%   flag = checkintersection(P, t); 
%   DESCRIPTION
%   This function checks self-intersecting meshes and return flag (0 -
%   good; 1 - bad)
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    edges = meshconnee(t); 
    [si11, si12]    = meshedgeintersect(P, t, P, edges);
    intersections = [];
    for m = 1:size(si12, 1)
        intersections = [intersections; si12{m}];
    end
    if ~isempty(intersections)       
        flag = 1;
    %         patch('Faces', t, 'Vertices', P, 'FaceColor', 'y', 'EdgeColor', 'k', 'FaceAlpha', 0.5);  
    %         hold on;        
    %         plot3(intersections(:, 1), intersections(:, 2), intersections(:, 3), 'r.', 'MarkerSize', 25);
    %         grid on; axis equal; axis tight;
    %         error('intersection check failed')
    else
        flag = 0; 
        display('The mesh has no self-intersections')
    end    
end