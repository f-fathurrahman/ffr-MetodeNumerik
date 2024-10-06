function [x,y] = rotation(xRel,yRel,xCOM,yCOM,theta)
% PURPOSE: Rotate a polygonal particle for a given 
%  angle counterclockwise
% OUTPUT: x,y The new absolute position of the corners
% USAGE: after the center of mass of a particle has
%  been shifted into the origin
% CAVEAT: The unit of "theta" is in radian
% LITERATURE: p. 252 
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

J = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)];

A = J * [xRel';yRel'];
xRel = A(1,:)';
yRel = A(2,:)';

% Absolute position 
x = xCOM - xRel;
y = yCOM - yRel;

return