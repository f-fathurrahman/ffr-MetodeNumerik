function [X,Y,U,V,T] = KUTTA(X,Y,U,V,T,H,F1,F2);

% ***** FOURTH-ORDER RUNGE-KUTTA FORMULAE (1.4.3) AND (1.1.13) FOR SOLVING
% TWO SIMULATIONS SECOND-ORDER DIFFERENTIAL EQUATIONS WRITTEN AS 
% DX/DT = U, DU/DT = F1(X,Y,U,V,T),
% DY/DT = V, DV/DT = F2(X,Y,U,V,T).
% H IS STEP SIZE IN T. F1, F2 ARE DEFINED IN OTHER SUBPROGRAMS.

D1X = H * U;
D1Y = H * V;
D1U = H * feval(F1, X, Y, U, V, T);
D1V = H * feval(F2, X, Y, U, V, T);

D2X = H * (U + D1U/2.);
D2Y = H * (V + D1V/2.);
D2U = H * feval(F1, X + D1X/2., Y + D1Y/2., U + D1U/2., V + D1V/2., T + H/2.);
D2V = H * feval(F2, X + D1X/2., Y + D1Y/2., U + D1U/2., V + D1V/2., T + H/2.);

D3X = H * (U + D2U/2.);
D3Y = H * (V + D2V/2.);
D3U = H * feval(F1, X + D2X/2., Y + D2Y/2., U + D2U/2., V + D2V/2., T + H/2.);
D3V = H * feval(F2, X + D2X/2., Y + D2Y/2., U + D2U/2., V + D2V/2., T + H/2.);

D4X = H * (U + D3U);
D4Y = H * (V + D3V);
D4U = H * feval(F1, X + D3X, Y + D3Y, U + D3U, V + D3V, T + H);
D4V = H * feval(F2, X + D3X, Y + D3Y, U + D3U, V + D3V, T + H);

T = T + H;
X = X + (D1X + 2.*D2X + 2.*D3X + D4X) / 6.;
Y = Y + (D1Y + 2.*D2Y + 2.*D3Y + D4Y) / 6.;
U = U + (D1U + 2.*D2U + 2.*D3U + D4U) / 6.;
V = V + (D1V + 2.*D2V + 2.*D3V + D4V) / 6.;

