function Pnew = rotatex(P, xangle)
%   SYNTAX
%   P = rotatex(P, xangle) where the angle is in degrees
%   DESCRIPTION
%   This function rotates the mesh about the x-axis
%   The local coordinate system may be used with the origin at the center of gravity
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    anglex = xangle/180*pi;
    %   Rotation about the x-axis
    LCY  = mean(P(:, 2)); LCY = 0;
    LCZ  = mean(P(:, 3)); LCZ = 0;
    Pnew(:, 1) = +P(:, 1);
    Pnew(:, 2) = +(P(:, 2) - LCY)*cos(anglex) - (P(:, 3) - LCZ)*sin(anglex);
    Pnew(:, 3) = +(P(:, 2) - LCY)*sin(anglex) + (P(:, 3) - LCZ)*cos(anglex);
    Pnew(:, 2) = Pnew(:, 2) + LCY;
    Pnew(:, 3) = Pnew(:, 3) + LCZ;
end

