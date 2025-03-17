function cfdReadSystem
%--------------------------------------------------------------------------
%
%  Written by the CFD Group @ AUB, Fall 2018
%  Contact us at: cfd@aub.edu.lb
%==========================================================================
% Routine Description:
%   This function read fvSchemes, fvSolution and controlDict files
%--------------------------------------------------------------------------
%

fprintf('\nReading System Dictionaries ...\n');

cfdReadControlDictFile;
cfdReadFvSchemesFile;
cfdReadFvSolutionFile;



