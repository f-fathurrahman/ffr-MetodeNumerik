function [F] = F(X,Y);

% ***** IT REPRESENTS RIGHT-HAND SIDE OF EQUATION (2.4.8) WITH THE DOUBLETS
% SHIFTED TO THE RIGHT THROUGH A SHORT DISTANCE *****

C = [0.15, 0.3, 0.2, 0.1, 0.05];
D = [-1.0, -0.5, 0., 0.5, 1.0];

F = Y;
for I=1:5
    F = F - C(I)*Y/((X-D(I)-1.0E-6)^2+Y^2); 
end