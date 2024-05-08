function [FY] = FY(X, Y, U, V, T);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.5.2) *****
% ( IN THE PRESENT PROBLEM THE FLUID IS STATIONARY )

global A B C D NU PI

UF = 0.0;
VF = 0.0;
WR = sqrt( (UF-U)^2 + (VF-V)^2 );
FY = ( -B + C*CD(WR)*(VF-V)*WR )/A;
