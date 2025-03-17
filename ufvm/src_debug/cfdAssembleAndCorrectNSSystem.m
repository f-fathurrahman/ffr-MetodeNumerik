function cfdAssembleAndCorrectNSSystem
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2018
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%
%--------------------------------------------------------------------------

% Momentum
cfdAssembleAndCorrectMomentumEquation;

% Continuity
cfdAssembleAndCorrectContinuityEquation;
