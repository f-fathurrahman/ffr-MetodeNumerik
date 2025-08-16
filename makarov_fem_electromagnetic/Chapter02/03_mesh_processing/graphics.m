function [] = graphics(P, t)
%   SYNTAX
%   graphics(P, t)
%   DESCRIPTION
%   This function is used by script viewer2.m to display multiple meshes
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2105, 1st ed.
    X = reshape(P(t', 1),[3, size(t, 1)]);
    Y = reshape(P(t', 2),[3, size(t, 1)]);
    Z = reshape(P(t', 3),[3, size(t, 1)]);
    patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);
    
    axis 'equal'; axis 'tight', set(gca, 'YDir', 'normal');
%     camproj('orthographic'); light('Position', [1 3 2], 'Style', 'local'); 
%     material shiny;
    view(70, 25); 


    NumberOfTirangles = size(t, 1)
    Q = min(simpqual(P, t))
    edges = meshconnee(t);
    temp = P(edges(:, 1), :) - P(edges(:, 2), :);
    minedgelength = min(sqrt(dot(temp, temp, 2)))
end
    





