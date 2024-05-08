function [CL]=CL(V);

% ***** LIFT COEFFICIENT OF THE AIRFOIL EXPRESSED AS A FUNCTION OF VERTICAL
% VERTICAL VELOCITY OF THE WING *****

global ALPHA0 BETA PI U

ALPHA = ALPHA0 - atan(V/U);
CL = 0.0;
if (abs(ALPHA) <= PI/10.)
    CL=2.*PI*ALPHA;
end

    