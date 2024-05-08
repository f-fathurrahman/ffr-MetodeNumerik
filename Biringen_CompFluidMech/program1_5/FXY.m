function [RES] = FXY(T,VEC);

global A B C D NU PI UF VF

X = VEC(1);
U = VEC(2);
Y = VEC(3);
V = VEC(4);

WR = sqrt( (UF-U)^2 + (VF-V)^2 );

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.5.1) *****
% ( IN THE PRESENT PROBLEM THE FLUID IS STATIONARY )

FX = C*CD(WR)*(UF-U)*WR/A;

% ***** IT REPRESENTS THE RIGHT-HAND SIDE OF EQUATION (1.5.2) *****
% ( IN THE PRESENT PROBLEM THE FLUID IS STATIONARY )

FY = ( -B + C*CD(WR)*(VF-V)*WR )/A;

RES = [U; FX; V; FY];
