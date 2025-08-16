function t = combdomain(P, t, strint)
%   SYNTAX
%   t = combdomain(P, t, strint)
%   DESCRIPTION
%   This function assigns domain number for every triangular circular subdomain
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    TriCenter = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) +  P(t(:, 3), :));

    t(:, 4) = 0;
    for m = 1:size(strint.R, 2)
        temp        = (TriCenter(:, 1) - strint.x(m)).^2 + ...
                      (TriCenter(:, 2) - strint.y(m)).^2 ;
        t(temp<strint.R(m)^2, 4) = m; 
    end
end

