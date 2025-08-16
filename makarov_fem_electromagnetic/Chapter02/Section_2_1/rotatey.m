function Pnew = rotatey(P, yangle)
%   SYNTAX
%   P = rotatey(P, yangle) where the angle is in degrees
%   DESCRIPTION
%   This function rotates the mesh about the y-axis
%   The local coordinate system may be used with the origin at the center of gravity
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    angley = yangle/180*pi;
    %   Rotation about the y-axis
    LCX  = mean(P(:, 1)); LCX = 0;
    LCZ  = mean(P(:, 3)); LCZ = 0;   
    Pnew(:, 1) = +(P(:, 1) - LCX)*cos(angley) + (P(:, 3) - LCZ)*sin(angley);
    Pnew(:, 2) = +P(:, 2);
    Pnew(:, 3) = -(P(:, 1) - LCX)*sin(angley) + (P(:, 3) - LCZ)*cos(angley);
    Pnew(:, 1) = Pnew(:, 1) + LCX;
    Pnew(:, 3) = Pnew(:, 3) + LCZ;
end

