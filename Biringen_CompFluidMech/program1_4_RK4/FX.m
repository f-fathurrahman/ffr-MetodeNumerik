function [FX] = FX(X, Y, U, V, T);

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.5.1) *****
% ( IN THE PRESENT PROBLEM THE FLUID IS STATIONARY )

global A B C D NU PI

UF = 0.0;
VF = 0.0;
WR = sqrt( (UF-U)^2 + (VF-V)^2 );
FX = C*CD(WR)*(UF-U)*WR/A;
