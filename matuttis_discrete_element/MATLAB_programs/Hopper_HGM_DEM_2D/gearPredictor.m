function [PX,VX,PY,VY,PF,VF]=gearPredictor(dt,freeParticle,X,Y,F)
% change order of PX, VX, PY, VY, PF, VF
% gearPredictor.m
% PURPOSE: Compute the six-point (fifth order) or
%  three-point (second order) predictor after Gear
%  ("Backward-difference formula")
%  see also gearCorrector
% INPUT: Timestep dt and corrected values X,Y,F 
%  from the previous timestep in SI-units
% OUTPUT: Predicted positions PX,PY,PF and velocities VX,VY,VF
%  in SI-units
%  XP YP FP (position/orientations and their higher 
%  time derivatives)in scaled units, via global variables
% USAGE: after the force calculation;
%  comment out either the fith order or the second 
%  order coefficients to change the order
% ALGORITHM: For the six-point form (positions and its
%  five subsequent time-derivatives)
%        | xp|   | 1 1 1 1 1 1 |   | x|
%        |vxp|   | 0 1 2 3 4 5 |   |vx|
%        |axp| = | 0 0 1 3 6 10| * |ax|
%        |bxp|   | 0 0 0 1 4 10|   |bx|.
%        |cxp|   | 0 0 0 0 1 5 |   |cx|
%        |dxp|   | 0 0 0 0 0 1 |   |dx|
%  Internally, the predictor uses dimensionless units, 
%  see p. 96, but the output is in SI-units
% CAVEAT: Higher order does not necessarily mean
%  higher accuracy, but higher stability. If high
%  accuracy for the positions is desired (to simulate
%  static heaps), the second order method is preferable.
%  p. ????
%  There is no stable Backward-difference formula
%  of higher than 5; All methods of higher order are
%  unconditionally unstable
% REVISION HISTORY:S. H. Ng, HG Matuttis 15-May-2014
% LITERATURE: 
%  advantages over other methods: p. 226
%  structure/formule: sec. 2.6
%  use of second or fifth order: sec. 7.6

global XP YP FP

% Set the higher order time derivatives of 
% particles which don't move according to 
% Newton's equation of motion  to zero
X(~freeParticle,2:6) = 0.0;
Y(~freeParticle,2:6) = 0.0;
F(~freeParticle,2:6) = 0.0;

% Prefactors to obtain the dimensionless units from
% the SI-units
idt = 1/dt;
pdt2 = 0.5 * dt * dt;  
pdt3 = pdt2 / 3.0 * dt;
pdt4 = pdt3 / 4.0 * dt;
pdt5 = pdt4 / 5.0 * dt;

% Uncomment this and comment out the next
% assignment of A if you want to work with 
% the 5th order method;  don't forget to
% change the corrector in gearCorrector.m, too:
% A fifth order corrector with a second order corrector
% would be of second order.
% A = [1 1 1 1 1 1;
%      0 1 2 3 4 5;
%      0 0 1 3 6 10;
%      0 0 0 1 4 10;
%      0 0 0 0 1 5;
%      0 0 0 0 0 1]; % 5th-order

A = [1 1 1 0 0 0;
     0 1 2 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
     0 0 0 0 0 0]; % 2nd-order

X = X';
Y = Y';
F = F';

% Successive scaled time derivatives (dimensionless units).
X(2,:) = X(2,:)*dt;
X(3,:) = X(3,:)*pdt2;
X(4,:) = X(4,:)*pdt3;
X(5,:) = X(5,:)*pdt4;
X(6,:) = X(6,:)*pdt5;

Y(2,:) = Y(2,:)*dt;
Y(3,:) = Y(3,:)*pdt2;
Y(4,:) = Y(4,:)*pdt3;
Y(5,:) = Y(5,:)*pdt4;
Y(6,:) = Y(6,:)*pdt5;

F(2,:) = F(2,:)*dt;
F(3,:) = F(3,:)*pdt2;
F(4,:) = F(4,:)*pdt3;
F(5,:) = F(5,:)*pdt4;
F(6,:) = F(6,:)*pdt5;

XP = A*X;
YP = A*Y;
FP = A*F;

% "Un-scaled" time derivatives (in SI-units).
PX = XP(1,:)';
VX = XP(2,:)'*idt;
PY = YP(1,:)';
VY = YP(2,:)'*idt;
PF = FP(1,:)';
VF = FP(2,:)'*idt;

return

